#! /usr/bin/env python

'''
Organize science files to be used by mkrefs.py to create reference files.
Convert a list of input files and required pipeline steps into
a list of strun commands that can be farmed out to condor.

Inputs:

inputs - 
         Astropy table of filenames, reference file types, and required
         calibration pipeline steps. 

input table example:

index    file_type    filename              req_steps   
1         gain      something_uncal.fits     bpm,sat
2         gain      another_uncal.fits       bpm,sat
3        readnoise  something_uncal.fits     bpm,sat,superbias

search_dir - 
             Directory to search for science files (from the 'filenames'
             column of the input table. Search looks into any 
             subdirectories as well, to allow the user flexibility when
             organizing input files 

output_dir - 
             Directory in which the pipeline-processed science files will
             be placed. This code does not actually run the pipeline, but
             only generates the commands needed to run the pipeline. The 
             output_dir will be added to these commands.

use_only_given - 
                 Use exactly the file listed in the 'filenames'
                 column of the input table. If this is False,
                 the code will search the search_dir to see 
                 what versions of the input file are present,
                 and choose the most appropriate version for the
                 requested pipeline processing state

overwrite_existing_files - 
                           If True, any existing files in the output
                           directory that have names matching the output
                           names will be over-written. This is essentially
                           a way to force the calibration pipeline to be 
                           re-run on a given input file
                
Outputs:

proc_table - 
             Copy of the input table with several columns added. 

             real_input_file:
             This added column lists the name of the file to be run through 
             the pipeline. This may differ from the filename originally given 
             as input, in the case where a version of the input file that has 
             already been partially processed is found and will be used to 
             save time. 

             output_name:
             This column lists the name of the output file that will be 
             produced by the pipeline

             steps_to_run:
             The column lists the pipeline steps that need to be run in order
             to transform the real_input_file into the output_name file.


strun - 
        This is a list of calls to the calibration pipeline. In order to use
        condor, we must use the command-line, strun, calls to the pipeline.


HISTORY:

Version 0.1 - Created by Hilbert, 11 Dec 2017.
'''
import os
import sys
import copy
from collections import OrderedDict
from glob import glob
import numpy as np
from astropy.io import ascii
from astropy.table import Column, Table
from jwst import datamodels

class CalibPrep:
    def __init__(self):
        self.__version__ = 0.1
        self.verbose = True
        self.inputs = None
        self.search_dir = ''
        self.use_only_given = False
        self.overwrite_existing_files = True
        self.output_dir = ''

        pipe_steps = [('dq','dq_init'),('sat','saturation'),('super','superbias'),
                      ('ref','refpix'),('ipc','ipc'),('lin','linearity'),
                      ('persistence','persistence'),('dark','dark_current'),
                      ('jump','jump'),('rampfit','rate')]
        self.pipe_step_dict = OrderedDict(pipe_steps)

        
    def choose_file(self,filelist,req,current):
        '''Given a list of files, their calibration state
        and the required calibration state, choose the best
        file to start with'''
        usefile = None
        first_false = -1
        for file in filelist:
            comp = []
            for key in req:
                comp.append(req[key] == current[file][key])
            comp = np.array(comp)
            f = np.where(np.array(comp) == False)[0]
            if ((len(f) > 0) and (f[0] > first_false)):
                first_false = f[0]
                usefile = file
            elif len(f) == 0:
                return file
        return usefile

        
    def completed_steps(self,file,checkheader=True):
        '''identify and return the pipeline steps completed
        for the input file'''
        done = copy.deepcopy(self.pipe_step_dict)
        data = datamodels.open(file)
        status = data.meta.cal_step._instance #returns dict of completed steps
        finsteps = list(status.keys())
        for key in self.pipe_step_dict:
            if self.pipe_step_dict[key] in finsteps:
                done[key] = True
            else:
                done[key] = False
        return done


    def create_output(self,base,req):
        '''Create the output name of the pipeline-processed
        input file given the list of required steps'''
        # base (e.g.) something_dq_init
        # Find the latest pipeline step that needs to be done
        # Keep track of any intemediate skipped steps (we
        # probably want to differentiate between files that have
        # had all pipeline steps done, compared to those which
        # share the latest completed step but may have skipped
        # intervening steps
        step = None
        suffix = 'uncal'
        skip = list(self.pipe_step_dict.values())
        baseend = len(base)
        for key in self.pipe_step_dict:
            if req[key]:
                step = key
                suffix = self.pipe_step_dict[key]
                skip.remove(self.pipe_step_dict[key])

                if self.pipe_step_dict[key] in base:
                    idx = base.index(self.pipe_step_dict[key])
                    if ((idx < baseend) and (idx != -1)):
                        baseend = copy.deepcopy(idx) - 1

        # Remove the entries in skip that are after the last
        # required pipeline step
        stepvals = np.array(list(self.pipe_step_dict.values()))
        lastmatch = np.where(stepvals == suffix)[0][0]
        for sval in stepvals[lastmatch+1:]:
            skip.remove(sval)
        
        true_base = base[0:baseend]
        true_base = true_base.replace('_uncal','')
        if self.verbose:
            print("True base name is {}".format(true_base))
            print("Suffix is {}".format(suffix))
            print("Skip vals are {}".format(skip))
        ofile = true_base + '_' + suffix
        if len(skip) > 0:
            for val in skip:
                ofile = ofile + '_' + val
        ofile = ofile + '.fits'
        return os.path.join(self.output_dir,ofile), true_base
        

    def file_search(self,base,dir):
        '''Search for all versions of a particular
        file in the given directory'''
        #files = glob(os.path.join(dir,base+'*fits'))

        # This version searches all subdirectories also...
        files = []
        for dirpath,dirnames,fnames in os.walk(dir,topdown=True):
            print(base)
            print(fnames)
            mch = [f for f in fnames if base in f]
            # or could try:
            # fnmatch.filter(fnames,base)
            for m in mch:
                files.append(os.path.join(dirpath,m))
        return files
                

    def find_repeats(self,tab):
        '''find repeated filenames in the input
        table
        CURRENTLY UNUSED'''
        names = []
        indices = []
        for line in tab:
            name = line['filename']
            index = tab['filename'].data == name
            names.append(name)
            indices.append(index)
        # Remove repeated entries
        final_names,findex = np.unique(names,return_index=True)
        final_indices = indices[findex]
        return final_names,final_indices

    
    def get_file_basename(self,file):
        '''Determine a given file's basename.
        Currently this is the full name minus
        the .fits at the end'''
        df = file.rfind('.fits')
        if df != -1:
            base = file[0:df]
        else:
            print(("WARNING: {} is not a valid fits file name."
                   .format(file)))
            sys.exit()
        return base
        
        
    def output_exist_check(self,filename):
        '''Check to see if the given file already
        exists. Remove if the user allows it.'''
        if os.path.isfile(filename):
            if self.overwrite_existing_files:
                try:
                    os.remove(filename)
                except:
                    print(("WARNING: unable to remove existing "
                           "filename {}".format(filename)))
                    sys.exit()
                if os.path.isfile(filename):
                    print(("WARNING: Removal of existing file "
                           "{} failed. Check permissions."
                           .format(filename)))
                    sys.exit()
            elif ((self.overwrite_existing_files == False) and
                  (self.output_dir == os.path.dirname(filename))):
                print(("WARNING: Output file {} already exists, "
                       "and overwrite_existing_files set to False."
                       .format(filename)))
                sys.exit()


    def prepare(self):
        '''Main function'''

        #find_repeats should located repeated entries in the
        #input table. but where should we do that, and how do
        #we proceed once we have the indices of the repeats?
        
        # Column of output names to add to table
        outfiles = []
        self.strun = []
        realinput = []
        all_to_run = []
        for line in self.inputs:
            file = line['filename']
            req_steps = self.step_dict(line['req_steps'])

            if self.verbose:
                print("")
                print("File: {}".format(file))
                print("Input required steps: {}".format(line['req_steps']))
                print("Required steps: {}".format(req_steps))
            
            # In order to search for other
            # versions of the file we need
            # to know the basename
            basename = self.get_file_basename(file)

            if self.verbose:
                print("Basename: {}".format(basename))
            
            # Create the output filename based on
            # the required pipeline steps
            outname, true_base = self.create_output(basename,req_steps)
            outfiles.append(outname)

            if self.verbose:
                print("Output name: {}".format(outname))

            # Check to see if a file with the same name as the
            # output filename already exists in the output directory.
            # If so and the user allows it, remove the file. If the user
            # has not allowed the removal of files, throw an error.
            # Similarly, if permissions prevent you from successfully
            # removing the file, throw an error
            self.output_exist_check(os.path.join(self.output_dir,outname))
            
            # Search
            if not self.use_only_given:
                files = self.file_search(true_base,self.search_dir)
            else:
                files = [os.path.join(self.search_dir,file)]

            if self.verbose:
                print("Found files: {}".format(files))
                
            # Determine the completed calibration steps
            # for each file
            current_state = {}
            for f in files:
                state = self.completed_steps(f)
                current_state[f] = state

                if self.verbose:
                    print("    file {}, current state: {}".format(f,current_state[f]))
                
            # Select the file to use
            if self.use_only_given:
                input_file = file
            else:
                input_file = self.choose_file(files,req_steps,current_state)

            realinput.append(input_file)
            if self.verbose:
                print("File to use: {}".format(input_file))
                
            # Create list of pipeline steps that must be
            # run on the file
            to_run = self.steps_to_run(input_file,req_steps,
                                       current_state[input_file])

            # Add column to input table
            tr = [key for key in to_run if to_run[key]]
            trstr = ''
            for s in tr:
                trstr = trstr + ',' + s
            if len(trstr) > 0:
                trstr = trstr[1:]
            else:
                trstr = 'None'
            all_to_run.append(trstr)
            
            if self.verbose:
                print("Steps that need to be run: {}".format(to_run))

            # Create strun command
            #command = self.strun_command(input_file,to_run,outname)
            #print('will need options for overriding reffiles here as well')
            #self.strun.append(command)

            #if self.verbose:
            #    print("Necessary strun command:")
            #    print("{}".format(command))
            
        # Add the output filename column to the input table
        realcol = Column(data=realinput,name='real_input_file')
        self.inputs.add_column(realcol)
        outcol = Column(data=outfiles,name='output_name')
        self.inputs.add_column(outcol)
        toruncol = Column(all_to_run,name='steps_to_run')
        self.inputs.add_column(toruncol)
        
        self.proc_table = copy.deepcopy(self.inputs)

        if self.verbose:
            print("Input table updated.")
            print(self.proc_table)
            ascii.write(self.proc_table,'test_out.txt',overwrite=True)

        # Turn the table into a series of strun commands
        self.strun = self.strun_command(realcol,toruncol,outcol)#,reffiles=??)

        if self.verbose:
            print(self.strun)
            
        c = Column(self.strun,name='strun_command')
        c_tab = Table()
        c_tab.add_column(self.inputs['index'])
        c_tab.add_column(c)

        if self.verbose:
            ascii.write(c_tab,'test_strun_commands.txt',overwrite=True)
        
        # Done
        # calling function can now use self.proc_table and
        # self.strun to access results


    def step_dict(self,stepstr):
        '''Translate the input list of required pipeline steps
        into a dictionary with boolean entries'''
        req = copy.deepcopy(self.pipe_step_dict)
        for key in req:
            req[key] = False
        stepslist = stepstr.split(',')
        for ele in stepslist:
            if ele not in list(req.keys()):
                print(("WARNING: unrecognized pipeline step: {}"
                       .format(ele)))
                sys.exit()
            try:
                req[ele] = True
            except:
                print(("WARNING: unrecognized required "
                       "pipeline step: {}".format(ele)))
                sys.exit()
        return req

    
    def steps_to_run(self,infile,req,current):
        '''Return a list of the pipeline steps that need
        to be run to bring the input file up to the state
        described by req. Return a dictionary of boolean
        values just like the req and current inputs'''
        torun = copy.deepcopy(current)
        for key in req:
            if req[key] == current[key]:
                torun[key] = False
            elif ((req[key] == True) & (current[key] == False)):
                torun[key] = True
            elif ((req[key] == False) & (current[key] == True)):
                print(("WARNING: Input file {} has had {} step run, "
                       "but the requirements say that it should not "
                       "be. Need a new input file.".format(infile,key)))
        return torun


    def strun_command(self,input,steps_to_run,outfile,overrides=[]):
        '''Create the necessary strun command to run the
        appropriate pipeline steps'''

        # Quick fix for ramp fitting, which uses a different
        # output suffix than step name, unlike the other steps
        step_names = copy.deepcopy(self.pipe_step_dict)
        step_names['rampfit'] = 'ramp_fitting'
        
        cmds = []
        
        # Determine appropriate reference file overrides
        #something

        initial = 'strun calwebb_detector1.cfg '
        for infile, steps, outfile in zip(input,steps_to_run,outfile):
            with_file = initial + infile
            
            if steps == 'None':
                cmd = 'None'
            else:
                skip_text = ''
                out_text = ''
                for key in step_names:

                    # Add skip statements for the steps to be skipped
                    if key not in steps:
                        skip_text += ' --steps.{}.skip=True'.format(step_names[key])

                    # Add override reference files

                # Add output filename
                finstep = steps.split(',')[-1]
                final_step = step_names[finstep]
                out_text += (' --steps.{}.output_file={}'
                             .format(final_step,outfile))
                
                # Put the whole command together
                cmd = with_file + skip_text + out_text #+override_text
            cmds.append(cmd)
        return cmds
        

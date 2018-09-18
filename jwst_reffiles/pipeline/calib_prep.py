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

index    reftype     fitsfile             ssbsteps
1         gain      something_uncal.fits     bpm,sat
2         gain      another_uncal.fits       bpm,sat
3        readnoise  something_uncal.fits     bpm,sat,superbias


colnames for input table, as defined in mkrefs.py
['cmdID', 'reftype', 'imlabel', 'imtype', 'detector', 'ssbsteps', 'imindex', 'imID', 'MJD', 'fitsfile']



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

import copy
from collections import OrderedDict
from glob import glob
import itertools
import os
# import re
import sys
import time

from astropy.io import ascii, fits
from astropy.table import Column, Table
import numpy as np
from jwst import datamodels

from jwst_reffiles.utils.definitions import PIPE_STEPS, PIPE_KEYWORDS


class CalibPrep:
    def __init__(self):
        self.__version__ = 0.1
        self.verbose = True
        self.inputs = None
        self.search_dir = ''
        self.use_only_given = False
        self.overwrite_existing_files = True
        self.output_dir = ''

        #pipe_steps = [('dq', 'dq_init'), ('sat', 'saturation'), ('super', 'superbias'),
        #              ('ref', 'refpix'), ('ipc', 'ipc'), ('lin', 'linearity'),
        #              ('persistence', 'persistence'), ('dark', 'dark_current'),
        #              ('jump', 'jump'), ('rampfit', 'rate')]
        self.pipe_step_dict = OrderedDict(PIPE_STEPS)

    def build_search_generator(self):
        '''Create a generator to use for searching for files

        Returns
        -------
        generator : generator
            Generator from os.walk
        '''
        if isinstance(self.search_dir, str):
            generator = [os.walk(self.search_dir, topdown=True)]
        elif isinstance(self.search_dir, list):
            # If the search directory is a list, then construct a generator
            # for each directory, and chain them together into a single generator
            generators = []
            for searchdir in self.search_dir:
                generators.append(os.walk(searchdir, topdown=True))

            # Chain together
            generator = itertools.chain()
            for gen in generators:
                generator = itertools.chain(generator, gen)
        return generator

    def choose_file(self, filelist, req, current):
        '''Given a list of files, their calibration state
        and the required calibration state, choose the best file
        with which to proceed

        Paramters
        ---------
        filelist : list
            List of fits filenames

        req : dict
            Dictionary of pipeline steps containing boolean entries to describe
            the required calibration state of the file

        current : dict
            Nested dictionary: for each filename (key), the value is a dictionary
            of pipeline steps containing boolean entires to describe the current
            calibration state of the file

        Returns
        -------
        usefile : str
            Name of the fits file to begin further processing with
        '''
        usefile = None
        first_false = -1
        for file in filelist:
            comp = []
            for key in req:
                comp.append(req[key] == current[file][key])
            comp = np.array(comp)
            f = np.where(np.array(comp) is False)[0]
            if ((len(f) > 0) and (f[0] > first_false)):
                first_false = f[0]
                usefile = file
            elif len(f) == 0:
                return file
        return usefile

    def completed_steps(self, input_file):
        '''Identify and return the pipeline steps completed
        for the input file

        Parameters
        ----------
        input_file : str
            Name of input fits file

        Returns
        -------
        completed : dict
        '''
        completed = copy.deepcopy(self.pipe_step_dict)
        for key in completed:
            completed[key] = False

        #starttime = time.time()
        #data = datamodels.open(file)
        #datamodels_time = time.time() - starttime
        #print("datamodels_time: {}".format(datamodels_time))
        #status = data.meta.cal_step._instance  # returns dict of completed steps
        #finsteps = list(status.keys())
        #for key in self.pipe_step_dict:
        #    if self.pipe_step_dict[key] in finsteps:
        #        done[key] = True
        #    else:
        #        done[key] = False
        header = fits.getheader(input_file)
        for key in header.keys():
            if key in PIPE_KEYWORDS.keys():
                value = header.get(key)
                if value == 'COMPLETE':
                    completed[PIPE_KEYWORDS[key]] = True
        return completed

    def create_output(self, base, req):
        '''Create the output name of the pipeline-processed
        input file given the list of required steps

        Parameters
        ----------
        base : str
            Base name of the fits file. Should be the filename minus '.fits'

        req : dict
            Dictionary of pipeline steps, with boolean entries describing whether
            those steps should be completed or not.

        Returns
        -------
        output_filename : str
            Full path + filename of the output once all the required pipeline steps
            have been run. Should be the true base name + the suffix associated with
            the latest pipeline step to be run.

        true_base : str
            Base name of the fits file. Should be the filename minus '.fits' as
            well as any pipeline-generated suffix
        '''
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

                # In the case where the filenamne has multiple pipeline step names attached,
                # walk back until we find the actual basename of the file
                if self.pipe_step_dict[key] in base:
                    idx = base.index(self.pipe_step_dict[key])
                    if ((idx < baseend) and (idx != -1)):
                        baseend = copy.deepcopy(idx) - 1

        # Remove the entries in skip that are after the last
        # required pipeline step
        stepvals = np.array(list(self.pipe_step_dict.values()))
        if suffix in self.pipe_step_dict.keys():
            lastmatch = np.where(stepvals == suffix)[0][0]
        elif suffix == 'uncal':
            lastmatch = -1
        else:
            raise IndexError("No entry {} in pipeline step dictionary.".format(suffix))
        for sval in stepvals[lastmatch+1:]:
            skip.remove(sval)

        # Find the true basename of the file
        true_base = base[0:baseend]
        true_base = true_base.replace('_uncal', '')
        if self.verbose:
            print("True base name is {}".format(true_base))
            print("Suffix is {}".format(suffix))
            print("Skip vals are {}".format(skip))
        true_base = os.path.split(true_base)[1]

        # Create the output filename by adding the name of the latest
        # pipeline step to be run
        ofile = true_base + '_' + suffix
        # if len(skip) > 0:
        #    for val in skip:
        #        ofile = ofile + '_' + val
        ofile = ofile + '.fits'
        output_filename = os.path.join(self.output_dir, ofile)
        return output_filename, true_base

    def file_search(self, base, generator_object):
        '''Search for all versions of a particular
        file in the given directory

        Parameters
        ----------
        base : str
            Basename of the file to search for. This should be the filename minus
            the .fits as well as any pipeline-created suffix.

        generator_object : obj
            Generator constructed from os.walk on all input search directories

        Returns
        -------
        files : list
            List of found files containing the input base string
        '''
        files = []
        for dirpath, dirnames, fnames in generator_object:
            mch = [f for f in fnames if base in os.path.join(dirpath, f)]
            # or could try:
            # fnmatch.filter(fnames,base)
            for m in mch:
                files.append(os.path.join(dirpath, m))
        return files

    def ORIGINAL_file_search_DELME_PROBABLY(self, base, dir):
        '''Search for all versions of a particular
        file in the given directory

        Parameters
        ----------
        base : str
            Basename of the file to search for. This should be the filename minus
            the .fits as well as any pipeline-created suffix.

        dir : str
            Directory to search. The search will look down the directory tree
            into all subdirectories below this level.

        Returns
        -------
        files : list
            List of found files containing the input base string
        '''
        # This version searches all subdirectories also...
        files = []
        for dirpath, dirnames, fnames in os.walk(dir, topdown=True):
            mch = [f for f in fnames if base in os.path.join(dirpath, f)]
            # or could try:
            # fnmatch.filter(fnames,base)
            for m in mch:
                files.append(os.path.join(dirpath, m))
        return files

    def find_repeats(self, tab):
        '''Find repeated filenames in the input
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
        final_names, findex = np.unique(names, return_index=True)
        final_indices = indices[findex]
        return final_names, final_indices

    def get_file_basename(self, input_file):
        '''Determine a given file's basename.
        Currently this is the full name minus
        the .fits at the end

        Parameters
        ----------
        input_file : str
            Input filename

        Returns
        -------
        base : str
            Base of the filename
        '''
        df = input_file.rfind('.fits')
        if df != -1:
            base = input_file[0:df]
        else:
            raise ValueError(("WARNING: {} is not a valid fits file name."
                              .format(input_file)))
        return base

    def output_exist_check(self, filename):
        '''Check to see if the given file already
        exists. Remove if the user allows it.

        Parameters
        ----------
        filename : str
            Name of file to look for
        '''
        if os.path.isfile(filename):
            if self.overwrite_existing_files:
                try:
                    os.remove(filename)
                except (FileNotFoundError, PermissionError) as error:
                    print(error)
            elif ((self.overwrite_existing_files is False) and
                  (self.output_dir == os.path.dirname(filename))):
                raise ValueError(("WARNING: Output file {} already exists, "
                                  "and overwrite_existing_files set to False."
                                  .format(filename)))

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

        #3print(self.inputs.colnames)
        # print(self.inputs)

        # Create generators of files in self.search_dir, to be searched later for
        # matching filenames
        self.search_dir_decompose()
        search_generator = self.build_search_generator()

        for line in self.inputs:
            file = line['fitsfile']

            #starttime = time.time()
            req_steps = self.step_dict(line['ssbsteps'])
            #req_steps_time = time.time() - starttime
            #print("req_steps_time: {}".format(req_steps_time))

            if self.verbose:
                print("")
                print("File: {}".format(file))
                print("Input required steps: {}".format(line['ssbsteps']))
                print("Required steps: {}".format(req_steps))

            # In order to search for other
            # versions of the file we need
            # to know the basename
            #starttime = time.time()
            basename = self.get_file_basename(file)
            #basename_time = time.time() - starttime
            #print("basename_time: {}".format(basename_time))

            if self.verbose:
                print("Basename: {}".format(basename))

            # Create the output filename based on
            # the required pipeline steps
            #starttime = time.time()
            outname, true_base = self.create_output(basename, req_steps)
            #outname_time = time.time() - starttime
            #print("outname_time: {}".format(outname_time))
            outfiles.append(outname)

            if self.verbose:
                print("Output name: {}".format(outname))

            # Check to see if a file with the same name as the
            # output filename already exists in the output directory.
            # If so and the user allows it, remove the file. If the user
            # has not allowed the removal of files, throw an error.
            # Similarly, if permissions prevent you from successfully
            # removing the file, throw an error
            #starttime = time.time()
            self.output_exist_check(os.path.join(self.output_dir, outname))
            #check_exist_time = time.time() - starttime
            #print("check_exist_time: {}".format(check_exist_time))

            # Search the self.search_dir directories for partially
            # processed files.
            if not self.use_only_given:
                #starttime = time.time()
                #print("true_base:", true_base)
                #print("search_dir", self.search_dir)
                files = self.file_search(true_base, search_generator)
                #if isinstance(self.search_dir, str):
                #    files = self.file_search(true_base, self.search_dir)
                #elif isinstance(self.search_dir, list):
                #    files = []
                #    for searchdir in self.search_dir:
                #        filesdir = self.file_search(true_base, searchdir)
                #        files += filesdir
                #file_search_time = time.time() - starttime
                #print("file_search_time: {}".format(file_search_time))
            else:
                files = [os.path.join(self.search_dir, file)]

            if len(files) == 0:
                print("No matching files found in {} for {}".format(self.search_dir, true_base))
                print("Falling back to the input file: {}".format(file))
                files = [file]

            if self.verbose:
                print("Found files: {}".format(files))

            # Determine the completed calibration steps
            # for each file
            # print("Files to check for completed steps: {}".format(files))

            current_state = {}
            #starttime = time.time()
            for f in files:
                state = self.completed_steps(f)
                current_state[f] = state

                if self.verbose:
                    print("    file {}, current state: {}".format(f, current_state[f]))
            #comp_steps_time = time.time() - starttime
            #print("comp_steps_time: {}".format(comp_steps_time))

            # Select the file to use
            if self.use_only_given:
                input_file = file
            else:
                #starttime = time.time()
                input_file = self.choose_file(files, req_steps, current_state)
                #choosefile_time = time.time() - starttime
                #print("choosefile_time: {}".format(choosefile_time))

            realinput.append(input_file)
            if self.verbose:
                print("File to use: {}".format(input_file))

            # Create list of pipeline steps that must be
            # run on the file
            #starttime = time.time()
            #print("current_State: {}".format(current_state))
            #print("input_file: {}".format(input_file))
            to_run = self.steps_to_run(input_file, req_steps,
                                       current_state[input_file])
            #stepstorun_time = time.time() - starttime
            #print("stepstorun_time: {}".format(stepstorun_time))

            # Add column to input table
            #starttime = time.time()
            tr = [key for key in to_run if to_run[key]]
            trstr = ''
            for s in tr:
                trstr = trstr + ',' + s
            if len(trstr) > 0:
                trstr = trstr[1:]
            else:
                trstr = 'None'
            all_to_run.append(trstr)
            #addcol_time = time.time() - starttime
            #print("addcol_time: {}".format(addcol_time))

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
        realcol = Column(data=realinput, name='real_input_file')
        self.inputs.add_column(realcol)
        outcol = Column(data=outfiles, name='output_name')
        self.inputs.add_column(outcol)
        toruncol = Column(all_to_run, name='steps_to_run')
        self.inputs.add_column(toruncol)

        self.proc_table = copy.deepcopy(self.inputs)

        if self.verbose:
            print("Input table updated.")
            print(self.proc_table)
            ascii.write(self.proc_table, 'test_out.txt', overwrite=True)

        # Turn the table into a series of strun commands
        #starttime = time.time()
        self.strun = self.strun_command(realcol, toruncol, outcol)  # ,reffiles=??)
        #make_command_time = time.time() - starttime
        #print("make_command_time: {}".format(make_command_time))

        if self.verbose:
            print(self.strun)

        c = Column(self.strun, name='strun_command')
        c_tab = Table()
        c_tab.add_column(self.inputs['cmdID'])
        c_tab.add_column(c)

        if self.verbose:
            ascii.write(c_tab, 'test_strun_commands.txt', overwrite=True)

        # Done
        # calling function can now use self.proc_table and
        # self.strun to access results

    def search_dir_decompose(self):
        '''Turn the self.search_dir into a list if necessary'''
        if ',' in self.search_dir:
            self.search_dir = [element.strip() for element in self.search_dir.split(',')]

    def step_dict(self, stepstr):
        '''Translate the input list of required pipeline steps
        into a dictionary with boolean entries

        Parameters
        ----------
        stepstr : str
            Comma-separated list of pipeline step abbreviations. These come from
            the origina config file.

        Returns
        -------
        req : dict
            Dictionary of pipeline steps with boolean values designating whether
            each step was requested in stepstr
        '''
        # Make a copy of pipeline step dictionary and
        # intialize all steps to False
        req = copy.deepcopy(self.pipe_step_dict)
        for key in req:
            req[key] = False
        if stepstr is None:
            return req

        # Strip all whitespace from the list of steps
        # and split into a list
        #pattern = re.compile(r'\s+')
        #stepstr = re.sub(pattern, '', stepstr)
        #stepslist = stepstr.split(',')
        stepslist = [element.strip() for element in stepstr.split(',')]
        for ele in stepslist:
            if ele not in list(req.keys()):
                raise ValueError(("WARNING: unrecognized pipeline step: {}"
                                  .format(ele)))
            try:
                req[ele] = True
            except KeyError as error:
                print(error)
        return req

    def steps_to_run(self, infile, req, current):
        '''Return a list of the pipeline steps that need
        to be run to bring the input file up to the state
        described by req. Return a dictionary of boolean
        values just like the req and current inputs

        Parameters
        ----------
        infile : str
            Input fits filename

        req : dict
            Dictionary of pipeline steps with boolean entries descibing whether
            each step should be complete on the final output file

        current : dict
            Dictionary of pipeline steps with boolean entries describing whether
            each step has already been completed on infile

        Returns
        -------
        torun : dict
            Dictionary of pipeline steps with boolean entires descrbing which
            required steps must still be run on infile
        '''
        torun = copy.deepcopy(current)
        for key in req:
            if req[key] == current[key]:
                torun[key] = False
            elif ((req[key] is True) & (current[key] is False)):
                torun[key] = True
            elif ((req[key] is False) & (current[key] is True)):
                print(("WARNING: Input file {} has had {} step run, "
                       "but the requirements say that it should not "
                       "be. Need a new input file.".format(infile, key)))
        return torun

    def strun_command(self, input, steps_to_run, outfile_name, overrides=[]):
        '''Create the necessary strun command to run the
        appropriate JWST calibration pipeline steps

        Parameters
        ----------
        input : astropy.table.Column
            astropy table column object listing input fits filenames

        steps_to_run : astropy.table.Column
            astropy table column object giving a comma-separated list of pipeline
            steps to run for each file

        outfile_name : astropy.table.Column
            astropy table column object listing pipeline ouput fits file name

        Returns
        -------
        cmds : list
            List of command line calls to the JWST calibration pipeline with appropriate
            flags to run/skip requested steps
        '''

        # Quick fix for ramp fitting, which uses a different
        # output suffix than step name, unlike the other steps
        step_names = copy.deepcopy(self.pipe_step_dict)
        step_names['rampfit'] = 'ramp_fitting'

        cmds = []

        # Determine appropriate reference file overrides
        # something

        initial = 'strun calwebb_detector1.cfg '
        for infile, steps, outfile in zip(input, steps_to_run, outfile_name):
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
                             .format(final_step, outfile))

                # Put the whole command together
                cmd = with_file + skip_text + out_text  # +override_text
            cmds.append(cmd)
        return cmds

#! /usr/bin/env python

'''
Search for a requested file with specificed pipeline steps applied,
within a given directory tree.
If the file exists, return the location. If it does not
exist, search for versions of the file with no or only earlier
pipeline steps applied. If those exist, return the location
and optionally call the appropriate pipeline steps to produce
the requested version of the file.

Assume that in the case where pipeline steps will be run to create
the requested file, that any arguments (such as override_XXXX, and
output_file for any intermediate steps) are placed in the provided 
config files. That way we don't need to keep track of overriding 
reference files here.

Version 0.0 - 20 Jan 2017 - BNH
'''

import os,sys
import argparse
from collections import OrderedDict
import fnmatch
from glob import glob
import copy
from jwst import datamodels
import numpy as np

class Produce:

    def __init__(self):
        self.verbose = False
        self.uncal_appendix = 'uncal'
        self.search_dir = None #top-level dir to search
        self.filename = '' #file to search for
        #make code able to handle directory and filename
        #inputs, or simply the full path filename
        self.checkheader = True
        self.run_pipeline = True
        self.outdir = None
        self.dq_cfg = 'dq_init.cfg'
        self.saturation_cfg = 'saturation.cfg'
        self.superbias_cfg = 'superbias.cfg'
        self.refpix_cfg = 'refpix.cfg'
        self.dark_cfg = 'dark.cfg'
        self.linearity_cfg = 'linearity.cfg'
        self.ipc_cfg = 'ipc.cfg'
        self.jump_cfg = 'jump.cfg'
        self.rampfit_cfg = 'ramp_fit.cfg'

        
    def returnBasename(self,filename):
        #given a filename with various pipeline-added
        #appended text, strip down and return the basename
        pass

    
    def findFile(self,filepatt,topdir):
        #search down through directory tree
        found = []
        for root,dirs,files in os.walk(topdir):
            for file in fnmatch.filter(files,filepatt):
                found.append(os.path.join(root,file))
        return found
    
    
    def findFileLocal(self,filepatt,directory):
        #search only in the given directory
        return glob(os.path.join(directory+filepatt))

    
    def getBaseName(self,filename):
        #strip off appended text and return only the basename of the file
        uc = filename.rfind(self.uncal_appendix)
        slash = filename.find('_',uc)
        return filename[0:uc+len(self.uncal_appendix)],filename[slash+1:]


    def completedSteps(self,file,checkheader=True):
        #identify and return the pipeline steps completed
        #for the input file
        appended = [('dq','dq_init'),('sat','saturation'),('super','superbias'),
                    ('ref','refpix'),('ipc','ipc'),('lin','linearity'),
                    ('dark','dark_sub'),('jump','jump'),('rampfit','rate')]
        appended = OrderedDict(appended)
        done = copy.deepcopy(appended)
        if checkheader == False:
            for key in appended:
                if appended[key] in file:
                    done[key] = True
                else:
                    done[key] = False
        else:
            data = datamodels.open(file)
            status = data.meta.cal_step._instance #returns dict of completed steps
            finsteps = list(status.keys())
            #print('finsteps is {}'.format(finsteps))
            for key in appended:
                if appended[key] in finsteps:
                    done[key] = True
                else:
                    done[key] = False
            #print('done is {}'.format(done))
        return done


    def getPipeSteps(self,files,checkheader):
        #given a list of different versions of a single file
        #(perhaps at various stages of pipeline)
        #find which pipeline steps have been run
        #on each file
        #print('files into getpipesteps are {}'.format(files))
        finished = []
        for file in files:
            #pipedone is an orderedDict [dq_init,saturation,superbias]
            pipedone = self.completedSteps(file,checkheader=checkheader)
            finished.append(pipedone)
            #print('for file {}'.format(file))
            #print('finished is {}'.format(finished))
        return finished

    
    def findClosestFile(self,found,doneSteps,requestedSteps):
        #from a list of input files, find the file that is closest
        #to the requested (searched for) file, in terms of completed
        #pipeline steps
        lastdone = max([i for i,x in enumerate(requestedSteps) if requestedSteps[x]])
        goodFile = []
        goodSteps = []
        goodLastDone = []

        print('List of found files is {}'.format(found))

        
        for foundfile,foundSteps in zip(found,doneSteps):
            compareSteps = [x&y for (x,y) in zip(foundSteps.values(), requestedSteps.values())]
            foundTrue = np.array([i for i,x in enumerate(foundSteps) if foundSteps[x]])
            if len(foundTrue) == 0:
                foundTrue = np.array([-1])
            firstFalse = min([i for i,x in enumerate(compareSteps) if x == False])

            #for testing
            #print('firstFalse is {}'.format(firstFalse))
            #print('lastdone is {}'.format(lastdone))
            #print('foundTrue is {}'.format(foundTrue))
            #print('')
            #print('')
                    
            if ((firstFalse <= lastdone) & (max(foundTrue-firstFalse) <= 0)):
                goodFile.append(foundfile)
                goodSteps.append(foundSteps)
                goodLastDone.append(max(foundTrue))

        #print('goodFile is {}'.format(goodFile))
        if len(goodFile) > 0:
            bestfile = np.where(goodLastDone == np.max(goodLastDone))[0][0]
            bestRunSteps = goodSteps[bestfile]
            print('File to start with is: {}'.format(goodFile[bestfile]))
            return goodFile[bestfile],bestRunSteps
        else:
            return None


    def add_dir_to_cfg(self,configfile,directory):
        #add the given path as the output directory in the config file
        flag = False
        str_to_add = "output_dir = '"+directory+"'"
        with open(configfile,'r+') as f:
            lines = f.readlines()
        for i,line in enumerate(lines):
            if 'output_dir' in line:
                lines[i] = str_to_add
                flag = True
        if flag == False:
            lines.append(str_to_add)

        with open(configfile,'w') as f:
            f.writelines(lines)
            
            

        
    def runPipeline(self,infile,infile_steps,requestedSteps,outdir):
        #beginning with the given file and orderedDict of completed
        #pipeline steps, run the appropriate steps to get the file
        #to a state that matches the requested pipeline steps
        from jwst.dq_init import DQInitStep
        from jwst.saturation import SaturationStep
        from jwst.superbias import SuperBiasStep
        from jwst.refpix import RefPixStep
        from jwst.ipc import IPCStep
        from jwst.linearity import LinearityStep
        from jwst.dark_current import DarkCurrentStep
        from jwst.jump import JumpStep
        from jwst.ramp_fitting import RampFitStep

        stepDict = {'dq':DQInitStep,'sat':SaturationStep,'super':SuperBiasStep,
                    'ref':RefPixStep,'ipc':IPCStep,'lin':LinearityStep,
                    'dark':DarkCurrentStep,'jump':JumpStep,'fit':RampFitStep}
        cfgDict = {'dq':self.dq_cfg,'sat':self.saturation_cfg,'super':self.superbias_cfg,
                    'ref':self.refpix_cfg,'ipc':self.ipc_cfg,'lin':self.linearity_cfg,
                    'dark':self.dark_cfg,'jump':self.jump_cfg,'fit':self.rampfit_cfg}
        appendDict = {'dq':'dq_init','sat':'saturation','super':'superbias',
                    'ref':'refpix','ipc':'ipc','lin':'linearity',
                    'dark':'dark','jump':'jump','rampfit':'rate'}

        
        ramp = datamodels.open(infile)

        infilename = os.path.basename(infile)
        indirname = os.path.dirname(infile)
        outfile = infilename[0:-5]
        
        #set the output directory
        if outdir is None:
            #print('outdir is none. Placing output file into')
            #print('{}, where input file was found.'.format(indirname))
            outdir = indirname

        #find which steps need to be run
        numToRun = 0
        for key in requestedSteps:
            if ((requestedSteps[key] == True) & (infile_steps[key] == False)):
                numToRun += 1
    
        counter = 0
        for key in requestedSteps:
            if ((requestedSteps[key] == True) & (infile_steps[key] == False)):
                counter += 1
                stepToRun = stepDict[key]
                cfgfile = cfgDict[key]

                #if outdir is provided, then open the config file and add it to
                #the output_file there.
                #if self.outdir is not None:
                #    self.add_dir_to_cfg(cfgfile,self.outdir)

                #append the appropriate string to the output filename
                outfile = outfile + '_' + appendDict[key]
                
                ramp = stepToRun.call(ramp,config_file=cfgfile)

                #Save the file after running the last pipeline step
                if counter == numToRun:
                    ramp.save(os.path.join(outdir,outfile+'.fits')) 
                
        
    def run(self):
    
        #break the input filename into a filename and a directory
        searchfile = os.path.basename(self.filename)
        searchpath = os.path.dirname(self.filename)

        #if the filename is not a full path, assume the user
        #wants to search starting from the current directory

        #this is the search for the exact file requested
        if len(searchpath) == 0:
            if self.search_dir is not None:
                searchpath = self.search_dir
            else:
                searchpath = './'
        else:
            if self.search_dir is not None:
                print("WARNING: full path name provided for file to search for,")
                print("({})".format(self.filename))
                print("but search_dir is also provided. ({})".format(self.search_dir)) 
                print("Ignoring file path info and searching in the search_dir.")
                searchpath = self.search_dir

        #Search for the exact file requested
        found = self.findFile(searchfile,searchpath)

        #If the exact file exists, return the location
        #Note that this can be a list if multiple instances
        #of the file exist
        if len(found) > 0:
            return found[0]

        else:
            #If the exact file isn't found, strip any pipeline-added
            #text from the filename and search on the basename with a
            #wildcard, to see what versions, if any, are present.
            filebase,donestring = self.getBaseName(searchfile)
            filebase = filebase + '*.fits'

            #print('searchpath is {}'.format(searchpath))
            found = self.findFile(filebase,searchpath)

            #if some version of the requested file is found, find what
            #pipeline steps need to be done and do them
            if len(found) > 0:

                #get a list of pipeline steps that are requested, based on
                #the original filename searched for
                requestedSteps = self.completedSteps(searchfile,checkheader=False)
                lastdone = max([i for i,x in enumerate(requestedSteps) if requestedSteps[x]])

                #start_file is the version of the file closest to that requested
                doneSteps = self.getPipeSteps(found,checkheader=self.checkheader)

                #compare the requested pipeline steps to the completed pipeline
                #steps and identify the file that is closest to the requested file
                start_file,start_file_steps = self.findClosestFile(found,doneSteps,requestedSteps)

                if start_file is None:
                    print("Search for {}: no files found.".format(self.filename))
                    return start_file
                
                #If requested, create the file we are searching for by running the
                #appropriate pipeline steps.
                if self.run_pipeline:
                    self.runPipeline(start_file,start_file_steps,requestedSteps,self.outdir)
                    #make sure you are returning the correct location of files....
                    if self.outdir is None:
                        self.outdir = './'
                    print("Search for {}: {} found and processed to make the requested file"
                          .format(self.filename,start_file))
                    return self.outdir+self.filename
                else:
                    #if the user chooses not to run the pipeline to make the
                    #requested file, then return the name of the file that is
                    #closest to the requested version
                    print("Search for {}: Best file found is {}.".format(self.filename,start_file))
                    return start_file

            else: #case where no versions of the requested file are found
                print("Search for {}: No files found.".format(self.filename))
                return None

            
    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparge.argument_parser(usage=usage,description='Find requested file and if not present, run pipeline to create it')
        parser.add_argument("filename",help="Name of file to search for")
        parser.add_argument("--search_dir",help="Top-level directory to search within. If not present, assume filename is a full path",default=None)
        parser.add_argument("--uncal_appendix",help="String giving the appendix used for raw filenames.",default="uncal")
        parser.add_argument("--checkheader",help="Look in file headers to see what pipeline steps have been run, rather than relying on appended strings in the filename.",action="store_true")
        parser.add_argument("-v","--verbose",help="Verbose mode.",action="store_true")
        parser.add_argument("--run_pipeline",help="Run pipeline to create the requested version of file",action="store_true")
        parser.add_argument("--outdir",help='Optional output directory for files created by the pipeline run',default=None)
        parser.add_argument("--dq_cfg",help="config file to use for dq_init pipeline step",default=None)
        parser.add_argument("--saturation_cfg",help="config file to use for saturation pipeline step",default=None)
        parser.add_argument("--superbias_cfg",help="config file to use for superbias pipeline step",default=None)
        parser.add_argument("--refpix_cfg",help="config file to use for refpix pipeline step",default=None)
        parser.add_argument("--dark_cfg",help="config file to use for dark subtraction pipeline step",default=None)
        parser.add_argument("--linearity_cfg",help="config file to use for linearity pipeline step",default=None)
        parser.add_argument("--ipc_cfg",help="config file to use for IPC pipeline step",default=None)
        parser.add_argument("--jump_cfg",help="config file to use for jump pipeline step",default=None)
        parser.add_argument("--rampfit_cfg",help="config file to use for rampfit pipeline step",default=None)
        return parser

if __name__ == '__main__':

    usagestring = "USAGE: find_file.py myfile_uncal_dq_init.fits"

    file = Produce()
    parser = file.add_options(usage = usagestring)
    args = parser.parse_args(namespace = file)

    file.run()
    

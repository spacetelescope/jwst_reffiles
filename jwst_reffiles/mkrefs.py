#!/usr/bin/env python
'''
create reference files for JWST instruments
A. Rest
'''
import argparse
import glob
import os,re,sys,types,copy,random
import logging
from logging import StreamHandler
from logging.handlers import RotatingFileHandler

import astropy.io.fits as fits
from astropy.io import ascii
import astropy
import numpy as np
import scipy

from jwst_reffiles.plugin_wrapper import mkrefclass_template
from jwst_reffiles.pipeline.calib_prep import CalibPrep
from jwst_reffiles.pipeline import pipeline_steps
from jwst_reffiles.utils.logging_functions import configure_logging, log_info, log_fail, make_log_file
from jwst_reffiles.utils.tools import astrotableclass, yamlcfgclass
from jwst_reffiles.utils.tools import makepath, executecommand, append2file, rmfile

# get the root dir of the code. This is needed only if the scripts are not installed as a module!
#if 'JWST_MKREFS_SRCDIR' in os.environ:
#    rootdir = os.environ['JWST_MKREFS_SRCDIR']
#    sys.path.extend([rootdir,
#                     '%s/gain' % rootdir,
#                     '%s/badpix_map' % rootdir])
#else:
#rootdir = os.path.dirname(os.path.realpath(__file__))


class cmdsclass(astrotableclass):
    def __init__(self, cmdtype, *args, verbose=0, debug=False, manual_log=None, **kwargs):
        #self.loginfo.append(('info', "Creating cmdsclass instance with cmdtype {}".format(cmdtype)))
        if manual_log is not None:
            self.manual_log = manual_log
            self.manual_log.append(('info', 'Instantiating cmdsclass {} with manual log'.format(cmdtype)))
        else:
            self.manual_log = None
            self.logger = logging.getLogger('jwst_reffiles.mkrefs.cmdsclass')
        astrotableclass.__init__(self, *args, **kwargs)
        self.verbose=verbose
        self.debug=debug
        self.cmdtype=cmdtype
        if not (cmdtype in ['ssb','ref','refmaster']):
            error_string = "cmd type {} not defined yet!".format(cmdtype)
            self.create_log_entry('error', error_string)
            raise RuntimeError(error_string)

    def create_log_entry(self, level, message, traceback=False):
        """Create a log entry for an instance of the class, dependong on
        whether the logger has been defined yet or not"""
        if self.manual_log is not None:
            self.manual_log.append((level, message))
        else:
            if level == 'debug':
                self.logger.debug(message)
            elif level == 'info':
                self.logger.info(message)
            elif level == 'warning':
                self.logger.warning(message)
            elif level == 'error':
                self.logger.error(message, exc_info=traceback)

    def filename_with_suffix(self, filename, addsuffix=None):
        if addsuffix != None:
            if re.search('%s$' % addsuffix,filename)==None:
                if re.search('\.$',filename)==None and re.search('^\.',addsuffix)==None: filename+='.'
                filename += addsuffix
        return(filename)

    def check_if_files_exists(self, file_col='output_name', file_exists_col='file_exists', addsuffix=None):
        """ Check if files in file_col exists, and fill the column file_exists_col with True or False """

        file_exists = np.full((len(self.t)), False)

        for i in range(len(self.t)):
            filename = self.filename_with_suffix(self.t[file_col][i], addsuffix=addsuffix)
            print('############# checking %s' % filename)
            if os.path.isfile(filename):
                file_exists[i] = True

        self.t[file_exists_col] = file_exists
        return(0)

    def check_if_already_in_batch(self, cmd_col='strun_command', already_in_batch_col='already_in_batch'):
        """ check if the reduced input files are already running in batch (Condor), and set 'already_in_batch' to True or False """

        print('\n*************************\n Check for Batch still needs to be implemented here!! \n *************************\n')

        already_in_batch = np.full((len(self.t)),False)
        for i in range(len(self.t)):
            # ADD CHECK HERE!!!
            #self.logger.info('Check for {}'.format(self.t[cmd_col][i]))
            self.create_log_entry('info', 'Check for {}'.format(self.t[cmd_col][i]))
            #already_in_batch[i]=True if in batch
        self.t[already_in_batch_col]=already_in_batch
        return(0)

    def pick_cmds_to_execute(self, execute_cmds=None, force_redo=False, maxNexe=None, execute_cmds_col='execute_cmds'):
        """ pick the commands to be executed!  execute_cmds can
        contain an initial one-dimensional input array of True or
        False. set execute_cmds to False if file_exists or
        already_in_batch=False.  if force_redo, then execute even if
        file_exists, and throw error message if already_in_batch """

        if execute_cmds is None:
            execute_cmds = np.full((len(self.t)), True)
        else:
            # just some error checking...
            if len(self.t) != len(execute_cmds):
                error_string = ("Length of initial execute_commands array is not the same than the length of "
                                "the cmd table ({}!={})").format(len(self.t), len(execute_cmds))
                self.create_log_entry('error', error_string, traceback=True)
                #self.logger.error(error_string)
                raise RuntimeError(error_string)

        Nexe = 0
        for i in range(len(self.t)):
            exeflag = execute_cmds[i]
            # Don't redo it if file already exists and if it is not force_redo
            if (not force_redo) and self.t['file_already_exists'][i]:
                try:
                    self.create_log_entry('info', ('File already exists and force redo is not set: {}')
                                          .format(self.t['output_name'].data[i]))
                    #self.logger.info(('File already exists and force redo is not set: {}')
                    #                 .format(self.t['output_name'].data[i]))
                except KeyError:
                    self.create_log_entry('info', ('File already exists and force redo is not set: {}')
                                          .format(self.t['outbasename'].data[i]))
                    #self.logger.info(('File already exists and force redo is not set: {}')
                    #                 .format(self.t['outbasename'].data[i]))
                exeflag = False

            if self.t['already_in_batch'][i]:
                exeflag = False
                # force_redo? Cannot do that, otherwise all hell would break loose with overwriting files etc
                if force_redo:
                    error_string = ("Cannot use --force_redo when command {} of index {} is still running in "
                                    "batch mode!").format(self.t['strun_command'][i], i)
                    self.create_log_entry('error', error_string, traceback=True)
                    #self.logger.error(error_string, exc_info=True)
                    raise RuntimeError(error_string)

            if maxNexe is not None and Nexe >= maxNexe:
                skip_string = 'Skipping executing command for index {}, since maxNexe={}'.format(i, maxNexe)
                self.create_log_entry('info', skip_string)
                #self.logger.info(skip_string)
                exeflag = False

            if exeflag:
                Nexe += 1

            execute_cmds[i] = exeflag

        self.t[execute_cmds_col] = execute_cmds
        return(0)

    def getlogfilenames(self, filename, logFlag=False, errorlogFlag=False):
        (basename, suffix) = os.path.splitext(filename)
        if suffix != '.fits':
            basename=filename

        if errorlogFlag:
            errorlog = '{}.err.txt'.format(basename)
        else:
            errorlog = None

        if logFlag:
            outlog = '{}.log.txt'.format(basename)
        else:
            outlog = None

        #if suffix != '.fits':
        #    self.create_log_entry('warning', 'It seems like the filename {} is not a fits file.'.format(filename))
            #self.logger.error('WARNING: It seems like the filename {} is not a fits file.'.format(filename))

        return(outlog, errorlog)

    def clean_old_output_files(self, indeces2run, outputfile_col='output_name', addsuffix=None):
        """ cleaning up old output files """

        outputfiles = self.t[outputfile_col][indeces2run]

        for outfile in outputfiles:
            outfile = self.filename_with_suffix(outfile,addsuffix=addsuffix)
            (logfilename, errlogfilename) = self.getlogfilenames(outfile, logFlag=True, errorlogFlag=True)
            rmfile(outfile)
            rmfile(logfilename)
            # Don't clean up errlogfilename, we want to keep track of old errors...
        return(0)

    def dummy_batch(self, batchfilename):
        """ dummy call to batch mode to return reasonable values """
        print("### run commands in batch mode: NOT YET IMPLEMENTED!!!")
        # this are just dummy return values until batch mode is implemented
        (errorflag, batchID) = (1, random.randint(1, 100000))

        return(errorflag, batchID)

    def run_cmds_batchmode(self, indeces2run, batchfilename, batchIDfilename=None,
                           logFlag=False, errorlogFlag=False,
                           cmds_col='strun_command', cmds_executed_col='strun_executed',
                           outputfile_col='output_name', addsuffix='fits'):

        # submit the batch here, return errorflag in case there is an issue submitting the batch
        errorflag = 0
        strun_list = self.t[cmds_col][indeces2run]

        # remove the old batchfile if exists!
        if os.path.isfile(batchfilename):
            os.remove(batchfilename)
            if os.path.isfile(batchfilename):
                raise RuntimeError('ERROR: could not save batch file {}'.format(batchfilename))

        # Save the cmds into the batch file
        self.create_log_entry('info', 'Saving cmds to {} for batch'.format(batchfilename))
        #self.logger.info('Saving cmds to {} for batch'.format(batchfilename))
        strun_list = [l+'\n' for l in strun_list]
        if open(batchfilename, 'w').writelines(strun_list):
            self.create_log_entry('error', 'ERROR: could not save batch file {}'.format(batchfilename))
            raise RuntimeError('ERROR: could not save batch file {}'.format(batchfilename))

        # this are just a dummy call to get an errorflag and random batchID. this will be replaces with the real batch call
        # ***************
        # need to be implemented: log and errorlog based on logFlag, errorlogFlag!
        #
        #for i in indeces2run:
        #    outfile = self.filename_with_suffix(self.t[outputfile_col][i],addsuffix=addsuffix)
        #    (logfilename,errlogfilename) = self.getlogfilenames(outfile, logFlag=logFlag, errorlogFlag=errorlogFlag)
        (errorflag, batchID) = self.dummy_batch(batchfilename)

        # append the batch ID into the batch ID file. This allows later on to test if batch jobs are still running!
        self.create_log_entry('info', 'Appending batch ID to {}'.format(batchIDfilename))
        #self.logger.info('Appending batch ID to {}'.format(batchIDfilename))
        if open(batchIDfilename, 'a').writelines(['{}\n'.format(batchID)]):
            self.create_log_entry('error', 'ERROR: could not save batch ID file {}'.format(batchIDfilename))
            raise RuntimeError('ERROR: could not save batch ID file {}'.format(batchIDfilename))

        # mark the entries that are submitted to the batch
        if errorflag:
            self.t[cmds_executed_col][indeces2run]='ERROR{}_batch'.format(errorflag)
        else:
            self.t[cmds_executed_col][indeces2run]='batch'

        return(0)

    def run_cmds_serial(self, indeces2run, logFlag=False, errorlogFlag=False,
                        cmds_col='strun_command', cmds_executed_col='strun_executed',
                        outputfile_col='output_name', addsuffix='fits'):

        if self.debug:
            self.create_log_entry('debug', '********* DEBUG!!!! ******\nJust touching the output files!!!')
            #self.logger.debug('********* DEBUG!!!! ******\nJust touching the output files!!!')
            for i in indeces2run:
                # outfile = self.filename_with_suffix(self.t[outputfile_col][i],addsuffix=addsuffix)
                outfile = self.t[outputfile_col][i]
                os.system("touch %s" % outfile)
                self.t[cmds_executed_col][i]='serial'
            return(0)

        for i in indeces2run:
            cmd = self.t[cmds_col][i]

            # decide if logfiles, depending on config filea
            #outfile = self.filename_with_suffix(self.t[outputfile_col][i],addsuffix=addsuffix)
            outfile = self.t[outputfile_col][i]
            (logfilename, errlogfilename) = self.getlogfilenames(outfile, logFlag=logFlag, errorlogFlag=errorlogFlag)

            self.create_log_entry('info', "Logfile to create for {} will be {}".format(outfile, logfilename))
            self.create_log_entry('info', "Error log to create for {} will be {}".format(outfile, errlogfilename))
            self.create_log_entry('info', '### Executing command for index {}: {}'.format(i, cmd))
            #self.logger.info("Logfile to create for {} will be {}".format(outfile, logfilename))
            #self.logger.info("Error log to create for {} will be {}".format(outfile, errlogfilename))
            #self.logger.info('### Executing command for index {}: {}'.format(i, cmd))

            (errorflag) = executecommand(cmd, '', cmdlog=logfilename, errorlog=errlogfilename)

            # extra error checking
            if not os.path.isfile(self.filename_with_suffix(outfile,addsuffix='fits')):
                self.create_log_entry('warning', 'File {} did not get created with strun command!'.format(self.filename_with_suffix(outfile,addsuffix='fits')))
                #self.logger.warning('ERROR: file {} did not get created with strun command!'.format(outfile))
                errorflag |= 2

            # fill table with results.
            if errorflag:
                self.t[cmds_executed_col][i]='ERROR{}_serial'.format(errorflag)
                self.create_log_entry('warning', ("In run_cmds_serial, Error Flag: {}. There was an error "
                                                  "running this command.".format(errorflag)))
                #self.logger.warning("Error Flag: {}. There was an error running this command.".format(errorflag))
            else:
                self.create_log_entry('info', "Error Flag: 0. Command completed successfully.")
                #self.logger.info("Error Flag: 0. Command completed successfully.")
                self.t[cmds_executed_col][i]='serial'

        return(0)


class mkrefsclass(astrotableclass):
    def __init__(self, logfile=False,counter=0):
        astrotableclass.__init__(self)
        self.loginfo = []

        # config file
        self.cfg = yamlcfgclass()

        self.verbose = 0
        self.debug = False
        self.onlyshow = False

        self.basedir = None
        self.basename = None
        self.outsubdir = None
        self.ssbdir = None
        self.runID = None
        self.runlabel = None

        # this allows access to some of the base functions (e.g., loading cfg file, default_optional_arguments)
        self.mkref_template = mkrefclass_template()

        self.imtable = astrotableclass()
        self.images4ssb = astrotableclass()
        #self.darks = None
        #self.flats = None

        self.ssbtable_excludecols4saving = ['repeat_of_index_number', 'index_contained_within']

        #
        self.DDtable = astrotableclass()
        self.FFtable = astrotableclass()
        self.DDFFtable = astrotableclass()

        self.ssbcmdtable = cmdsclass('ssb', manual_log=self.loginfo)
        self.refcmdtable = cmdsclass('ref', manual_log=self.loginfo)
        self.refmastercmdtable = cmdsclass('refmaster', manual_log=self.loginfo)

        # self.allowed_reflabels = ['example_bpm', 'bpm', 'rdnoise_nircam', 'gain_armin']
        mkref_file_list = glob.glob('{}/*/mkref_*.py'.format(os.path.dirname(__file__)))
        mkref_file_basenames = [os.path.basename(entry) for entry in mkref_file_list]
        self.allowed_reflabels = [entry.replace('mkref_', '').replace('.py', '') for entry in mkref_file_basenames]
        self.reflabel_directories = [os.path.dirname(entry) for entry in mkref_file_list]
        self.mkref_file_list = copy.deepcopy(mkref_file_list)

        print('Allowed reflabels are {}'.format(self.allowed_reflabels))
        #print('reflabel_directories:',self.reflabel_directories)
        #print('mkref_file_list:',self.mkref_file_list)

        self.loginfo.append(('info', 'Allowed reflabels are {}'.format(self.allowed_reflabels)))
        self.loginfo.append(('info', 'Reflabel directories are {}'.format(self.reflabel_directories)))

    def define_options(self, parser=None, usage=None, conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
        parser.add_argument("reflabels_and_imagelist", nargs='+', help=("list of ref types to be done"
                                                                        "and image (or file patterns) lists"))
        parser.add_argument('-b', '--batchmode', help="run the commands in batch mode (default=%(default)s)",
                            action="store_true", default=False)
        parser.add_argument('--onlyshow', help="only show what would be done, but don't do it."
                                               "(default=%(default)s)", action="store_true", default=False)

        parser.add_argument('--outrootdir', default=None,
                            help=('output root directory. outrootdir[/outsubdir][/runID]/reflabel/'
                                  'reflabel[_outsubdir][_runID][.addsuffix].cmdID.reftype.fits '
                                  '(default=%(default)s)'))
        parser.add_argument('--outsubdir', default=None,
                            help=('subdir added to the output root directory (and filename) '
                                  '(default=%(default)s)'))
        parser.add_argument('--addsuffix', default=None,
                            help='suffix added to the output basename (default=%(default)s)')

        parser.add_argument('--runID', default=None,
                            help=('runID subdir added to output root directory (and filename) '
                                  '(default=%(default)s)'))
        parser.add_argument('--runIDNdigits', default=None,
                            help='how many leading zeros in runID (default=%(default)s)')
        parser.add_argument('-n', '--newrunID', help="get the next available runID",
                            action="store_true", default=False)
        parser.add_argument('--skip_runID_ssbdir', help="Don't use runID in the default ssb dir",
                            action="store_true", default=None)

        parser.add_argument('--ssbdir', default=None,
                            help='Directory into which ssb files are stored (default=%(default)s)')

        parser.add_argument('--force_redo_strun', help="Redo the ssb reduction with strun even if file already exists",
                            action="store_true", default=False)
        parser.add_argument('--maxNstrun', type=int, default=None,
                            help='limit the number of strun commands to be run')

        parser.add_argument('--force_redo_ref', help="Redo the reference file even if file already exists",
                            action="store_true", default=False)
        parser.add_argument('--maxNref', type=int, default=None,
                            help='limit the number of reference file commands to be run')
        parser.add_argument('--logFileLevel', help="Minimum message level to add to log file", default="info")
        parser.add_argument('--logScreenLevel', help="Minimum message level to display to screen", default="warning")

        ###
        ### get the options from the individual mkref_X.py commands!
        ###
        # this gets the default optional paramters (verbose, debug, cfg file etc) from the mkref_template class
        self.mkref_template.default_optional_arguments(parser)

        # Loop through all allowed reflabels, and add the extra options
        self.mkrefpackages = {}
        for reflabel, refdir in zip(self.allowed_reflabels, self.reflabel_directories):
            dir_end = refdir.split('/')[-1]
            mkref_label_name = 'mkref_{}'.format(reflabel)
            self.mkrefpackages[mkref_label_name] = __import__("jwst_reffiles.{}.{}".
                                                              format(dir_end, mkref_label_name),
                                                              globals(), locals(), ['mkrefclass'], 0)
            mkref = self.mkrefpackages[mkref_label_name].mkrefclass()
            if reflabel != mkref.reflabel:
                self.loginfo.append(('warning', ('reflabel={} in script {} is inconsistent with mkref '
                                                 'script name {}!'.format(mkref.reflabel, mkref_label_name,
                                                                          mkref_label_name))))
                raise RuntimeError('reflabel={} in script {} is inconsistent with mkref script name {}!'
                                   .format(mkref.reflabel, mkref_label_name, mkref_label_name))
            # get the reflabel specific options
            mkref.extra_optional_arguments(parser)
        return(parser)

    def loadcfgfiles(self, *pargs, **kwargs):
        return(self.mkref_template.loadcfgfiles(*pargs, **kwargs))

    def getbasedir(self, outrootdir=None, outsubdir=None, runlabel=None):
        """ basedir based on options and cfg file: outrootdir[/outsubdir][/runID]"""
        if outrootdir is None:
            if self.cfg.params['output']['outrootdir'] == '':
                raise RuntimeError("No output rootdir!")
            else:
                outrootdir = self.cfg.params['output']['outrootdir']

        basedir = os.path.abspath(os.path.expandvars(outrootdir))
        basename = ''

        # outsubdir: test if outsubdir is passed or if it is specified in the config file ...
        if outsubdir is None and self.cfg.params['output']['outsubdir'] != '':
            outsubdir = self.cfg.params['output']['outsubdir']
        # ... and if yes, add it
        if outsubdir is not None:
            basedir += '/%s' % outsubdir
            if basename != '':
                basename += '_'
            basename += outsubdir

        if runlabel is not None:
            basedir += '/%s' % runlabel
            if basename != '':
                basename += '_'
            basename += runlabel

        if basename == '':
            basename = 'B'

        return(basedir, '%s/%s' % (basedir, basename), outsubdir)

    def getbasename(self):
        s = self.basedir+'/'
        if self.outsubdir is not None:
            s += '%s/' % self.outsubdir
        if self.runID is not None:
            s += '%s/' % self.outsubdir

    def getssbdir(self, ssbdir, basedir, skip_runID_ssbdir=False):
        """ get the ssbdir. It can be hardcoded to a directory with
        --ssbdir or in the cfg file. If it is not specified, then
        basedir/ssb is used"""

        if skip_runID_ssbdir==None: skip_runID_ssbdir = self.cfg.params['output']['skip_runID_ssbdir']

        print('CCCC',skip_runID_ssbdir)

        if ssbdir is not None and ssbdir != '':
            return(ssbdir)
        elif self.cfg.params['output']['ssbdir'] is not None and self.cfg.params['output']['ssbdir'] != '':
            return(self.cfg.params['output']['ssbdir'])
        else:
            if skip_runID_ssbdir:
                print('CCC',self.basedir)
                dir = re.sub('\/run\d+$','',self.basedir)
                print('CCC',dir)
                return('%s/ssb' % dir)
            else:
                return('%s/ssb' % self.basedir)
        return(None)

    def get_highest_runID(self, basedir):
        """ find all subdirs witn run\d+ within the passed basedir, and extract the highest runID"""
        files = glob.glob('%s/run*' % basedir)
        if len(files) == 0:
            self.loginfo.append(('info', 'In get_highest_runID: No runIDs yet!'))
            return(None)
        runIDs = []
        runIDpattern = re.compile('run(\d+)$')
        for fname in files:
            # print(os.path.basename(fname))
            m = runIDpattern.search(fname)
            if m is not None:
                runID = int(m.groups()[0])
                runIDs.append(runID)
        if len(runIDs) > 0:
            runID = max(runIDs)
        else:
            runID = None
        return(runID)

    def set_runID(self, outrootdir=None, outsubdir=None, runID=None, runIDNdigits=None, newrunID=False):
        """ set the runID depending on the options and config file """

        # if no runID is specified with the options, then use the one from the cfg file
        if runID is None:
            runID = mkrefs.cfg.params['output']['runID']

        if runID == '' or runID is None:
            # If no runID specified, make sure it is set to None!
            self.runID = None
        elif isinstance(runID, int) or runID.isdigit():
            # if runID is an integer, set it to it!
            self.runID = int(runID)
        else:
            # if runID is AUTO, get the current highest runID, and increment if --new_runID (-n)
            if runID.upper() == 'AUTO':
                (basedir, basename, subdir) = self.getbasedir(outrootdir=outrootdir, outsubdir=outsubdir,
                                                              runlabel=None)
                self.runID = self.get_highest_runID(basedir)
                if self.runID is None:
                    self.runID = 0
                else:
                    if newrunID:
                        self.runID += 1
            else:
                raise RuntimeError("Incorrect runID=%s requested!" % runID)

        if self.runID is not None:
            if runIDNdigits is None:
                runIDNdigits = self.cfg.params['output']['runIDNdigits']
            if runIDNdigits is None:
                runIDNdigits = 1

            b = 'run%0'+'%d' % runIDNdigits+'d'
            self.runlabel = b % (self.runID)
        else:
            self.runlabel = ''

        return(self.runID, self.runlabel)

    def set_dirs(self, outrootdir=None, outsubdir=None, runID=None, runIDNdigits=1, newrunID=False,
                 ssbdir=None, skip_runID_ssbdir=False):
        # get the runID
        mkrefs.set_runID(outrootdir=outrootdir, outsubdir=outsubdir,
                         runID=runID, runIDNdigits=runIDNdigits, newrunID=newrunID)

        (self.basedir, self.basename, self.outsubdir) = self.getbasedir(outrootdir=outrootdir,
                                                                        outsubdir=outsubdir,
                                                                        runlabel=self.runlabel)

        self.ssbdir = self.getssbdir(ssbdir, self.basedir, skip_runID_ssbdir=skip_runID_ssbdir)

        # Create the directories if they don't already exist
        if not self.onlyshow:
            state_basedir = makepath(self.basedir)
            state_ssbdir = makepath(self.ssbdir)

        self.loginfo.append(('info', '#### Directories:'))
        self.loginfo.append(('info', 'basedir: {}'.format(self.basedir)))
        self.loginfo.append(('info', 'basename: {}'.format(self.basename)))
        self.loginfo.append(('info', 'ssbdir: {}'.format(self.ssbdir)))
        self.loginfo.append(('info', 'outsubdir: {}'.format(self.outsubdir)))
        self.loginfo.append(('info', 'runID: {}'.format(self.runID)))
        self.loginfo.append(('info', 'runlabel: {}'.format(self.runlabel)))

        return(0)

    def getrefoutbasename(self, reflabel, reftype, cmdID):
        outbasename = '{}/{}/{}.{}.cmd{}_{}'.format(self.basedir,reflabel,reflabel,os.path.basename(self.basename),cmdID,reftype)
        return(outbasename)

    def trim_imagelist(self, imagelist, basenamepattern=None):
        '''
        only keep fits file that match basenamepattern
        '''
        self.logger.info("#### Trimming images...")
        # print imagelist
        if imagelist is None or len(imagelist) == 0:
            self.logger.info('Nothing to do, no images!!!')
            return(imagelist)

        if basenamepattern is not None:
            newimagelist = []
            basenamepattern_compiled = re.compile(basenamepattern)
            # make sure the input files are all ok!
            for i in range(len(imagelist)):
                m = basenamepattern_compiled.search(imagelist[i])
                if m is None:
                    self.logger.info('SKIPPING {} in trim_imagelist'.format(imagelist[i]))
                else:
                    #if len(m.groups())==0:
                    #    raise RuntimeError,"%s is matching %s, but does not return basename" % (basenamepattern_compiled.pattern,imagelist[i])
                    #basename = m.groups()[0]
                    #print 'VVV',imagelist[i]
                    newimagelist.append(imagelist[i])
            if len(imagelist) != len(newimagelist):
                if self.verbose > 1:
                    self.logger.info('skipping {} out of {} images, {} left'.format(len(imagelist)-len(newimagelist),
                                                                         len(imagelist), len(newimagelist)))
            imagelist = newimagelist

        return(imagelist)

    def parse_reflabels_images(self, reflabels_and_imagelist, basenamepattern=None):
        reflabellist = []
        imagelist = []
        for s in reflabels_and_imagelist:
            if s in self.cfg.params['reflabels']:
                reflabellist.append(s)
            else:
                if not os.path.isfile(s):
                    self.logger.error("ERROR: file {} does not exist, thus not a viable input file".format(s))
                    raise RuntimeError("ERROR: file {} does not exist, thus not a viable input file".format(s))
                imagelist.append(os.path.abspath(s))

        imagelist = self.trim_imagelist(imagelist, basenamepattern)
        return(reflabellist, imagelist)

    def getimtypes(self):
        if not ('imtype' in self.imtable.t.colnames):
            self.imtable.t['imtype'] = None

        darkpattern = re.compile(self.cfg.params['inputfiles']['dark_pattern'])
        flatpattern = re.compile(self.cfg.params['inputfiles']['flat_pattern'])

        for i in range(len(self.imtable.t)):
            shortfilename = os.path.basename(self.imtable.t['fitsfile'][i])
            if darkpattern.search(shortfilename):
                self.imtable.t['imtype'][i] = 'dark'
            elif flatpattern.search(shortfilename):
                self.imtable.t['imtype'][i] = 'flat'
            else:
                self.logger.warning("ERROR: image type of image {} is unknown!".format(shortfilename))
                raise RuntimeError('ERROR: image type of image %s is unknown!' % shortfilename)

    def getimageinfo(self, imagelist, dateobsfitskey=None, timeobsfitskey=None, mjdobsfitskey=None):
        self.logger.info("Getting image info")

        # self.imtable['fitsfile'].format('%s')
        self.imtable.t['fitsfile'] = imagelist
        self.imtable.t['imID'] = list(range(len(imagelist)))
        self.imtable.t['imtype'] = None
        self.imtable.t['skip'] = False

        self.logger.info("Table self.imtable created")

        requiredfitskeys = self.cfg.params['inputfiles']['requiredfitskeys']
        if requiredfitskeys is None:
            requiredfitskeys = []
        if type(requiredfitskeys) == bytes:
            requiredfitskeys = [requiredfitskeys]
        if dateobsfitskey is not None:
            requiredfitskeys.append(dateobsfitskey)
        if timeobsfitskey is not None:
            requiredfitskeys.append(timeobsfitskey)
        if mjdobsfitskey is not None:
            requiredfitskeys.append(mjdobsfitskey)
            mjdcol = mjdobsfitskey
        else:
            mjdcol = 'MJD'

        self.imtable.fitsheader2table('fitsfile',
                                      requiredfitskeys=requiredfitskeys,
                                      optionalfitskey=self.cfg.params['inputfiles']['optionalfitskeys'],
                                      raiseError=True, skipcolname='skip')
        self.imtable.dateobs2mjd(dateobsfitskey, mjdcol, mjdobscol=mjdobsfitskey, timeobscol=timeobsfitskey)

        self.getimtypes()
        # sort by MJD
        self.imtable.t.sort('MJD')

        #self.darks = self.imtable.t[np.where(self.imtable.t['imtype']=='dark')]
        #self.flats = self.imtable.t[np.where(self.imtable.t['imtype']=='flat')]
        return(0)

    #@log_fail
    #@log_info
    def organize_inputfiles(self, reflabels_and_imagelist):
        # parse teh command line arguments for reflabels and images
        (reflabellist, imagelist) = self.parse_reflabels_images(reflabels_and_imagelist,
                                                                basenamepattern=self.cfg.params['inputfiles']['basenamepattern'])
        self.reflabellist = reflabellist

        # make the image table and populate it with info. Also get teh darks and flats table
        self.getimageinfo(imagelist,
                          dateobsfitskey=self.cfg.params['inputfiles']['dateobs_fitskey'],
                          timeobsfitskey=self.cfg.params['inputfiles']['timeobs_fitskey'],
                          mjdobsfitskey=self.cfg.params['inputfiles']['mjdobs_fitskey'])

        self.detectors = set(self.imtable.t['DETECTOR'])

        #if self.verbose:
        #    print('#################\n### %d images found!' % len(self.imtable.t))
        #    print('### %d darks, %d flats' % (len(np.where(self.imtable.t['imtype'] == 'dark')[0]),
        #                                      len(np.where(self.imtable.t['imtype'] == 'flat')[0])))
        #    print('### %d detectors:' % (len(self.detectors)), ", ".join(self.detectors))
        #    if self.verbose > 1:
        #        print(self.imtable.t)

        self.logger.info('#################\n### {} images found!'.format(len(self.imtable.t)))
        self.logger.info('### {} darks, {} flats'.format(len(np.where(self.imtable.t['imtype'] == 'dark')[0]),
                                                     len(np.where(self.imtable.t['imtype'] == 'flat')[0])))
        self.logger.info('### {} detectors:'.format((len(self.detectors)), ", ".join(self.detectors)))
        columns = ''
        for col in self.imtable.t.colnames:
            columns = columns + '  {}'.format(col)
        self.logger.info(columns)
        for row in self.imtable.t:
            row_str = ''
            for col in self.imtable.t.colnames:
                row_str = row_str + '  {}'.format(row[col])
            self.logger.info(row_str)

        return(0)

    def check_inputfiles(self):
        self.logger.info('basedir: {}'.format(self.basedir))
        self.logger.info('basename: {}'.format(self.basename))

        # test if output dir already exists.
        #if not os.path.isdir(self.basedir):
        #    # if it doesn't exist: write the image list into it if !onlyshow!
        #    if not self.onlyshow:
        #        if self.verbose > 1:
        #            print('Creating directory %s' % self.basedir)
        #        os.makedirs(self.basedir)
        #        if not os.path.isdir(self.basedir):
        #            raise RuntimeError('Cannot create directory %s' % self.basedir)
        #    else:
        #        print('*** only showing: directory %s would be created here!' % self.basedir)

        imtbl_filename = '%s.im.txt' % self.basename
        self.logger.info('Image table filename: {}'.format(imtbl_filename))

        # check: does the saved inputim file has the same input images
        # than the current list? If yes, then all is good, if not
        # throw an error!
        if os.path.isfile(imtbl_filename):
            imtable_before = astrotableclass()
            imtable_before.load(imtbl_filename)

            if (len(self.imtable.t) != len(imtable_before.t)):
                warning_string = ('Length {} of new list of input images is different than the one from '
                                  'the saved table ({}), cannot proceed! either get a new runID with -n, '
                                  'or remove the files in {}'.format(len(self.imtable.t),
                                                                     len(imtable_before.t),
                                                                     self.basedir))
                self.logger.error(warning_string)
                raise RuntimeError(warning_string)
            for i in range(len(self.imtable.t['fitsfile'])):
                if self.imtable.t['fitsfile'][i] != imtable_before.t['fitsfile'][i]:
                    warning_string = ('Input file {} in new list is different than previous input '
                                      'file {}'.format(self.imtable.t['fitsfile'][i],
                                                       imtable_before.t['fitsfile'][i]))
                    self.logger.error(warning_string)
                    raise RuntimeError(warning_string)

            self.logger.info('New input list matches previous input list! Continuing...')

        self.imtable.write(imtbl_filename, verbose=True, clobber=True)

        # just some paranoia...
        if not os.path.isfile(imtbl_filename):
            error_string = 'Could not write {}! permission issues?'.format(imtbl_filename)
            self.logger.error(error_string)
            raise RuntimeError(error_string)

        return(0)

    def get_optional_arguments(self, args, sysargv):
        fitspattern = re.compile('\.fits$')

        opt_arg_list = []
        for i in range(1, len(sysargv)):
            if sysargv[i] in args.reflabels_and_imagelist:
                if sysargv[i] in self.reflabellist:
                    self.logger.debug('this is a reflabel', sysargv[i])
                else:
                    # test if it is an image
                    if not fitspattern.search(sysargv[i]):
                        self.logger.debug('not a fits file! filepattern?', sysargv[i])
            else:
                opt_arg_list.append(sysargv[i])
        self.logger.debug('optional arguments:', opt_arg_list)
        return(opt_arg_list)

    def getDlist(self, detector):
        '''
        returns list of Dark indeces, where the indices refer to the self.darks table
        '''

        #self.imtable.t['skip'][7]=True
        # indices for dark frames
        dindex, = np.where(self.imtable.t['imtype'] == 'dark')
        # indices for detector and not skipped
        dindex4detector = dindex[np.where(np.logical_and(self.imtable.t['DETECTOR'][dindex] == detector,
                                                         np.logical_not(self.imtable.t['skip'][dindex])))]
        self.logger.info('Possible {} Darks for detector {}'.format(len(dindex4detector), detector))
        self.logger.info(self.imtable.t[dindex4detector])

        #D = astrotableclass(names=('D1index','D1imID'),dtype=('i4', 'i4'))
        D = astrotableclass()
        D.t['D1index'] = dindex4detector
        D.t['D1imID'] = self.imtable.t['imID'][dindex4detector]
        D.t['D1index', 'D1imID'].format = '%d'
        #print D.t
        #sys.exit(0)
        return(D, ('D1',))

    def getDDlist(self, detector, max_Delta_MJD=None):
        '''
        returns list of Dark-Dark pair indeces,  where the indices refer to the self.imtable table
        '''
        self.logger.debug('# Getting DD list')
        #self.imtable.t['skip'][7]=True
        #print  self.imtable.t[6:11]

        # indices for dark frames
        dindex, = np.where(self.imtable.t['imtype'] == 'dark')
        # indices for detector and not skipped
        dindex4detector = dindex[np.where(np.logical_and(self.imtable.t['DETECTOR'][dindex] == detector,
                                                         np.logical_not(self.imtable.t['skip'][dindex])))]
        self.logger.info('Possible %d Darks for detector %s' % (len(dindex4detector), detector))
        self.logger.info(self.imtable.t[dindex4detector])

        DD = astrotableclass(names=('D1index', 'D2index', 'D1imID', 'D2imID'), dtype=('i4', 'i4', 'i4', 'i4'))
        i = 0
        while i < len(dindex4detector)-1:
            if max_Delta_MJD is not None:
                self.logger.info('Checking if imID={} and {} can be DD'.format(self.imtable.t['imID'][dindex4detector[i]],
                                                                           self.imtable.t['imID'][dindex4detector[i+1]]))
                dMJD = self.imtable.t['MJD'][dindex4detector[i+1]]-self.imtable.t['MJD'][dindex4detector[i]]
                self.logger.debug('dMJD:', dMJD)
                if dMJD > max_Delta_MJD:
                    self.logger.info(('Skipping imID={} (MJD={}) since imID={} is not within timelimit (Delta '
                                  'MJD = {}>{})!'.format(self.imtable.t['imID'][dindex4detector[i]],
                                                         self.imtable.t['MJD'][i],
                                                         self.imtable.t['imID'][dindex4detector[i+1]],
                                                         dMJD, max_Delta_MJD)))
                    i += 1
                    continue
            self.logger.info('Adding DD pair with imID={} and {}'.format(self.imtable.t['imID'][dindex4detector[i]],
                                                                     self.imtable.t['imID'][dindex4detector[i+1]]))
            DD.t.add_row({'D1index': dindex4detector[i], 'D2index': dindex4detector[i+1],
                          'D1imID': self.imtable.t['imID'][dindex4detector[i]],
                          'D2imID': self.imtable.t['imID'][dindex4detector[i+1]]})
            i += 2

        self.logger.debug(DD.t)
        return(DD, ('D1', 'D2'))

    def getDpluslist(self, detector, max_Delta_MJD=None):
        '''
        returns all Dark indeces, where the indices refer to the self.imtable table
        '''
        self.logger.info('# Getting D+ list')

        # indices for dark frames
        dindex, = np.where(self.imtable.t['imtype'] == 'dark')
        # indices for detector and not skipped
        dindex4detector = dindex[np.where(np.logical_and(self.imtable.t['DETECTOR'][dindex] == detector,
                                                         np.logical_not(self.imtable.t['skip'][dindex])))]
        self.logger.info('{} Darks for detector {}'.format(len(dindex4detector), detector))
        for goodrow in self.imtable.t[dindex4detector]:
            row_str = ''
            for col in self.imtable.t.colnames:
                row_str = row_str + '  {}'.format(goodrow[col])
            self.logger.info(row_str)

        if len(dindex4detector)<1:

            Dplus = astrotableclass()
            return(Dplus,())

        cols = ['D%dindex' % i for i in range(1,len(dindex4detector)+1)]
        cols.extend(['D%dimID' % i for i in range(1,len(dindex4detector)+1)])
        dtypes = ['i4' for i in range(len(cols))]
        imagelabels = ['D%d' % i for i in range(1,len(dindex4detector)+1)]

        Dplus = astrotableclass(names=cols, dtype=dtypes)
        r =  Dplus.t.add_row({})

        for i in range(1,len(dindex4detector)+1):
            Dplus.t['D%dindex' % i][r]=dindex4detector[i-1]
            Dplus.t['D%dimID' % i][r]=self.imtable.t['imID'][dindex4detector[i-1]]

        return(Dplus,imagelabels)

    def getFpluslist(self, detector, max_Delta_MJD=None):
        '''
        returns all Flat indeces, where the indices refer to the self.imtable table
        '''
        self.logger.info('# Getting FF list')

        # indices for dark frames
        findex, = np.where(self.imtable.t['imtype'] == 'flat')
        # indices for detector and not skipped
        findex4detector = findex[np.where(np.logical_and(self.imtable.t['DETECTOR'][findex] == detector,
                                                         np.logical_not(self.imtable.t['skip'][findex])))]
        self.logger.info('{} Flats for detector {}'.format(len(findex4detector), detector))
        for goodrow in self.imtable.t[findex4detector]:
            row_str = ''
            for col in self.imtable.t.colnames:
                row_str = row_str + '  {}'.format(goodrow[col])
            self.logger.info(row_str)

        if len(findex4detector)<1:

            Fplus = astrotableclass()
            return(Fplus,())

        cols = ['F%dindex' % i for i in range(1,len(findex4detector)+1)]
        cols.extend(['F%dimID' % i for i in range(1,len(findex4detector)+1)])
        dtypes = ['i4' for i in range(len(cols))]
        imagelabels = ['F%d' % i for i in range(1,len(findex4detector)+1)]

        Fplus = astrotableclass(names=cols, dtype=dtypes)
        r =  Fplus.t.add_row({})

        for i in range(1,len(findex4detector)+1):
            Fplus.t['F%dindex' % i][r]=findex4detector[i-1]
            Fplus.t['F%dimID' % i][r]=self.imtable.t['imID'][findex4detector[i-1]]

        return(Fplus,imagelabels)

    def getFFlist(self, detector, max_Delta_MJD=None):
        '''
        returns list of Flat-Flat pair indeces, where the indices refer to the self.imtable table
        '''
        self.logger.info('# Getting FF list')

        # indices for dark frames
        findex, = np.where(self.imtable.t['imtype'] == 'flat')
        # indices for detector and not skipped
        findex4detector = findex[np.where(np.logical_and(self.imtable.t['DETECTOR'][findex] == detector,
                                                         np.logical_not(self.imtable.t['skip'][findex])))]
        self.logger.info('Possible {} Flats for detector {}'.format(len(findex4detector), detector))
        self.logger.info(self.imtable.t[findex4detector])

        FF = astrotableclass(names=('F1index', 'F2index', 'F1imID', 'F2imID'), dtype=('i4', 'i4', 'i4', 'i4'))
        i = 0
        while i < len(findex4detector)-1:
            if max_Delta_MJD is not None:
                self.logger.info('Checking if imID={} and {} can be FF'.format(self.imtable.t['imID'][findex4detector[i]],
                                                                           self.imtable.t['imID'][findex4detector[i+1]]))
                dMJD = self.imtable.t['MJD'][findex4detector[i+1]]-self.imtable.t['MJD'][findex4detector[i]]
                self.logger.debug('dMJD:', dMJD)
                if dMJD > max_Delta_MJD:
                    self.logger.info(('Skipping imID={} (MJD={}) since imID={} is not within '
                                  'timelimit (Delta MJD = {}>{})!').format(self.imtable.t['imID'][findex4detector[i]],
                                                                           self.imtable.t['MJD'][i],
                                                                           self.imtable.t['imID'][findex4detector[i+1]],
                                                                           dMJD, max_Delta_MJD))
                    i += 1
                    continue
            self.logger.info('Adding FF pair with imID={} and {}'.format(self.imtable.t['imID'][findex4detector[i]],
                                                                     self.imtable.t['imID'][findex4detector[i+1]]))
            FF.t.add_row({'F1index': findex4detector[i], 'F2index': findex4detector[i+1],
                          'F1imID': self.imtable.t['imID'][findex4detector[i]],
                          'F2imID': self.imtable.t['imID'][findex4detector[i+1]]})
            i += 2

        self.logger.debug(FF.t)
        return(FF, ('F1', 'F2'))

    def getDDFFlist(self, detector, DD_max_Delta_MJD=None, FF_max_Delta_MJD=None, DDFF_max_Delta_MJD=None):
        '''
        returns list of Flat-Flat pair indeces, where the indices refer to the self.imtable table
        '''
        logger.info('\n### Getting DDFF list')
        DD, DDimtypes = self.getDDlist(detector, max_Delta_MJD=DD_max_Delta_MJD)
        FF, FFimtypes = self.getFFlist(detector, max_Delta_MJD=FF_max_Delta_MJD)
        logger.info('### DD and FF lists created, now matching them!!!')

        DDFF = astrotableclass(names=('F1index', 'F2index', 'F1imID', 'F2imID', 'D1index', 'D2index', 'D1imID',
                                      'D2imID'), dtype=('i4', 'i4', 'i4', 'i4', 'i4', 'i4', 'i4', 'i4'))
        ddcount = np.zeros(len(DD.t))

        for f in range(len(FF.t)):
            logger.debug('# Finding DD pair for FF pair with imID=%d and %d' % (FF.t['F1imID'][f], FF.t['F2imID'][f]))
            if DDFF_max_Delta_MJD is not None:
                FF_MJDmin = np.amin(np.array((self.imtable.t['MJD'][FF.t['F1index'][f]],
                                              self.imtable.t['MJD'][FF.t['F2index'][f]])))
                FF_MJDmax = np.amax(np.array((self.imtable.t['MJD'][FF.t['F1index'][f]],
                                              self.imtable.t['MJD'][FF.t['F2index'][f]])))
                logger.debug('FF MJDs: {} {}'.format(FF_MJDmin, FF_MJDmax))

            d_best = None
            ddcountmin = None
            for d in range(len(DD.t)):
                if ddcountmin is None or ddcount[d] < ddcountmin:
                    if DDFF_max_Delta_MJD is not None:
                        logger.info('# Testing DD pair with imID=%d and %d' % (DD.t['D1imID'][d], DD.t['D2imID'][d]))
                        DD_MJDmin = np.amin(np.array((self.imtable.t['MJD'][DD.t['D1index'][d]],
                                                      self.imtable.t['MJD'][DD.t['D2index'][d]])))
                        DD_MJDmax = np.amax(np.array((self.imtable.t['MJD'][DD.t['D1index'][d]],
                                                      self.imtable.t['MJD'][DD.t['D2index'][d]])))
                        logger.info('DD MJD range', DD_MJDmin, DD_MJDmax)
                        dMJD = np.fabs([DD_MJDmax-FF_MJDmin, DD_MJDmin-FF_MJDmax, DD_MJDmax-FF_MJDmax,
                                        DD_MJDmin-FF_MJDmin])
                        max_dMJD = np.amax(dMJD)
                        if max_dMJD > DDFF_max_Delta_MJD:
                            logger.info('DD pair with imID=%d and %d cannot be used, dMJD=%f>%f' % (DD.t['D1imID'][d],
                                                                                                    DD.t['D2imID'][d],
                                                                                                    max_dMJD, DDFF_max_Delta_MJD))
                            continue
                    ddcountmin = ddcount[d]
                    d_best = d

            if d_best is None:
                logger.info('SKIPPING following FF pair since there is no matching DD pair!')
                logger.info(FF.t[f])
            else:
                logger.info('SUCCESS! DD pair with imID %d and %d found for FF pair with imID %d and %d' % (DD.t['D1imID'][d_best],
                                                                                                             DD.t['D2imID'][d_best],
                                                                                                             FF.t['F1imID'][f],
                                                                                                             FF.t['F2imID'][f]))
                ddcount[d_best] += 1
                DDFF.t.add_row({'D1index': DD.t['D1index'][d_best],
                                'D2index': DD.t['D2index'][d_best],
                                'D1imID': DD.t['D1imID'][d_best],
                                'D2imID': DD.t['D2imID'][d_best],
                                'F1index': FF.t['F1index'][f],
                                'F2index': FF.t['F2index'][f],
                                'F1imID': FF.t['F1imID'][f],
                                'F2imID': FF.t['F2imID'][f]})

        return(DDFF, ('D1', 'D2', 'F1', 'F2'))

    def get_inputimage_sets(self, reflabel, detector, DD_max_Delta_MJD=None, FF_max_Delta_MJD=None,
                            DDFF_max_Delta_MJD=None):
        imtypes = self.cfg.params[reflabel]['imtypes']
        self.logger.debug("imtypes is {}".format(imtypes))
        imagesets = []
        if imtypes == 'D':
            imagesets, imagelabels = self.getDlist(detector)
        elif imtypes == 'DD':
            imagesets, imagelabels = self.getDDlist(detector, max_Delta_MJD=DD_max_Delta_MJD)
        elif imtypes == 'D+':
            imagesets, imagelabels = self.getDpluslist(detector)
        elif imtypes == 'FF':
            imagesets, imagelabels = self.getFFlist(detector, max_Delta_MJD=FF_max_Delta_MJD)
        elif imtypes == 'F+':
            imagesets, imagelabels = self.getFpluslist(detector)
        elif imtypes == 'DDFF' or imtypes == 'FFDD':
            imagesets, imagelabels = self.getDDFFlist(detector, DD_max_Delta_MJD=DD_max_Delta_MJD,
                                                      FF_max_Delta_MJD=FF_max_Delta_MJD,
                                                      DDFF_max_Delta_MJD=DDFF_max_Delta_MJD)
        else:
            raise RuntimeError("ERROR: imtypes=%s not yet implemented!" % imtypes)
        return(imagesets, imagelabels)

    def check_if_primary_strun(self, primary_strun_col='primary_strun'):
        """ check which entries are primary strun commands. """
        primary_strun = np.full((len(self.ssbcmdtable.t)),True)
        for i in range(len(self.ssbcmdtable.t)):
            if self.ssbcmdtable.t['index_contained_within'][i][0]!=-1:
                primary_strun[i]=False
            if self.ssbcmdtable.t['repeat_of_index_number'][i][0]!=-1:
                primary_strun[i]=False
        self.ssbcmdtable.t[primary_strun_col]=primary_strun
        return(0)

    #@log_fail
    #@log_info
    def mk_ssb_cmds(self, force_redo_strun=False, maxNstrun=None):
        '''
        Construct the reflabel commands, get the input files
        '''
        self.logger.info('\n##################################\n### Constructing commands\n##################################')

        self.ssbcmdtable.verbose = self.verbose
        self.ssbcmdtable.debug = self.debug

        # Expand ssbstep list if a shorthand is used
        for reflabel in self.reflabellist:
            step_name = self.cfg.params[reflabel]['ssbsteps'].strip()
            if step_name[-1] == '-':
                self.logger.info('Expanding step name {} for {}'.format(step_name, self.cfg.params['instrument']))
                self.cfg.params[reflabel]['ssbsteps'] = pipeline_steps.step_minus(step_name, self.cfg.params['instrument'])
            elif step_name[-1] == '+':
                self.logger.info('Expanding step name {} for {}'.format(step_name, self.cfg.params['instrument']))
                self.cfg.params[reflabel]['ssbsteps'] = pipeline_steps.step_plus(step_name, self.cfg.params['instrument'])

        # this is just to get the correct dtype for the reflabel  columns
        #dummyreflabel = astrotableclass()
        #dummyreflabel.t['reflabel'] = self.reflabellist
        #dummyreflabel.t['ssbsteps'] = [self.cfg.params[reflabel]['ssbsteps'] for reflabel in self.reflabellist]

        self.refcmdtable = cmdsclass('ref',
                                     verbose=self.verbose, debug=self.debug,
                                     names=('reflabel', 'detector', 'cmdID', 'Nim', 'outbasename'),
                                     dtype=(np.dtype(object),np.dtype(object), 'i8', 'i8',np.dtype(object)))
        #self.refcmdtable = astrotableclass(names=('reflabel','detector','cmdID','Nim'))
        self.refcmdtable.t['cmdID', 'Nim'].format = '%5d'
        self.refcmdtable.t['reflabel', 'detector', 'outbasename'].format = '%s'

        #print(self.imtable.t['DETECTOR'].dtype)
        #sys.exit(0)

        self.inputimagestable = astrotableclass(names=('cmdID', 'reflabel', 'imlabel', 'imtype', 'detector',
                                                       'ssbsteps', 'imindex', 'imID', 'MJD', 'fitsfile'),
#                                                dtype=('i8', dummyreflabel.t['reflabel'].dtype, 'S40',
                                                dtype=('i8', np.dtype(object), np.dtype(object),
                                                       self.imtable.t['imtype'].dtype,
                                                       self.imtable.t['DETECTOR'].dtype,
                                                       np.dtype(object),
                                                       'i8', 'i8', 'f8', self.imtable.t['fitsfile'].dtype))
        self.inputimagestable.t['cmdID', 'imindex', 'imID'].format = '%5d'
        self.inputimagestable.t['MJD'].format = '%.8f'

        cmdID = 0
        for reflabel in self.reflabellist:

            self.logger.info('SSB steps for reflabel {}: {}'.format(reflabel, self.cfg.params[reflabel]['ssbsteps']))

            counter = 0
            for detector in self.detectors:
                self.logger.info('\n#############################\n### Constructing {} commands for detector {}'
                             .format(reflabel, detector))
                inputimagesets, inputimagelabels = self.get_inputimage_sets(reflabel, detector,
                                                                            DD_max_Delta_MJD=self.cfg.params['DD']['max_Delta_MJD'],
                                                                            FF_max_Delta_MJD=self.cfg.params['FF']['max_Delta_MJD'],
                                                                            DDFF_max_Delta_MJD=self.cfg.params['DDFF']['max_Delta_MJD'])

                self.logger.info('### {} image sets for {} for detector {}'.format(len(inputimagesets.t),
                                                                               reflabel, detector))
                self.logger.debug(inputimagesets.t)
                for i in range(len(inputimagesets.t)):
                    self.logger.debug('image set index {}'.format(i))
                    for inputimagelabel in inputimagelabels:

                        # This is the index to the image in imtable
                        imindex = inputimagesets.t['%sindex' % inputimagelabel][i]

                        # sanity test: imindex and imID need to agree!!
                        if inputimagesets.t['%simID' % inputimagelabel][i] != self.imtable.t['imID'][imindex]:
                            runtime_string = ('BUG!!! This should not happen! the imIDs {} and {} need to '
                                              'match!').format(inputimagesets.t['{}imID'.format(inputimagelabel)][i],
                                                               self.imtable.t['imID'][imindex])
                            self.logger.error(runtime_string)
                            raise RuntimeError(runtime_string)

                        self.logger.debug(self.imtable.t[inputimagesets.t['{}index'.format(inputimagelabel)][i]])

                        dict2add = {'cmdID': cmdID, 'imindex': imindex, 'reflabel': reflabel,
                                    'detector': detector, 'imlabel': inputimagelabel,
                                    'ssbsteps': self.cfg.params[reflabel]['ssbsteps']}
                        for key in ['fitsfile', 'imID', 'imtype', 'MJD']:
                            dict2add[key] = self.imtable.t[key][imindex]
                        self.inputimagestable.t.add_row(dict2add)

                    self.refcmdtable.t.add_row({'reflabel': reflabel,
                                                'detector': detector,
                                                'cmdID': cmdID,
                                                'Nim': len(inputimagelabels),
                                                'outbasename': self.getrefoutbasename(reflabel,
                                                                                      self.cfg.params[reflabel]['reftype'],
                                                                                      cmdID)})
                    cmdID += 1

        self.logger.debug('\n### COMMANDS:')
        self.logger.debug('\n',self.refcmdtable.t)
        self.logger.debug('\n### INPUT FILES:')
        self.logger.debug('\n',self.inputimagestable.t)

        print('**** ADD: additional check if input images are the same!! ****')
        self.inputimagestable.write('%s.inputim.txt' % self.basename, verbose=True, clobber=True)

        mmm = CalibPrep(self.cfg.params['instrument'])
        mmm.inputs = self.inputimagestable.t

        print('Bryan: we need to look for ssb files in the ssb output dir, and then also in the optional pipeline_prod_search_dir. Let me know how I should pass this inof!')
        mmm.search_dir = self.ssbdir
        # add additional search dirs!
        if self.cfg.params['output']['pipeline_prod_search_dir'] is not None:
            mmm.search_dir += ',%s' % (self.cfg.params['output']['pipeline_prod_search_dir'])

        mmm.output_dir = self.ssbdir
        mmm.prepare()

        #print('BACK IN MKREFS:')
        # print(mmm.proc_table['index', 'cmdID', 'reflabel', 'steps_to_run', 'repeat_of_index_number', 'index_contained_within'])
        #print(mmm.proc_table['strun_command'][-1])

        self.ssbcmdtable.verbose = self.verbose
        self.ssbcmdtable.t = mmm.proc_table

        # check which entries are primary strun commands. Only primary
        # strun commands need to be executed, secondary strun commands
        # are already covered by primary strun commands
        self.check_if_primary_strun()

        # check if the reduced input files exist or not, and fill
        # 'file_already_exists' column with True or False
        self.ssbcmdtable.check_if_files_exists(file_col='output_name', file_exists_col='file_already_exists')

        # check if the reduced input files are already running in
        # batch (Condor), and set 'already_in_batch' to True or False
        self.ssbcmdtable.check_if_already_in_batch()

        # pick the strun commands that need to be executed
        self.ssbcmdtable.pick_cmds_to_execute(execute_cmds=copy.deepcopy(self.ssbcmdtable.t['primary_strun']),
                                              force_redo=force_redo_strun, maxNexe=maxNstrun,
                                              execute_cmds_col='execute_strun')

        # set strun_executed to None
        self.ssbcmdtable.t['strun_executed'] = np.full((len(self.ssbcmdtable.t)), None)

        self.logger.info('ssb table colnames: {}'.format(self.ssbcmdtable.t.colnames))
        self.ssbcmdtable.write('%s.ssbcmds.txt' % self.basename, verbose=True, clobber=True,
                               exclude_names=self.ssbtable_excludecols4saving)

        #self.refcmdtable.write('%s.refcmds.txt' % self.basename,verbose=True,clobber=True)

        #self.logger.info(self.ssbcmdtable.t['index', 'real_input_file', 'ssbsteps', 'repeat_of_index_number',
        #                                    'index_contained_within', 'primary_strun', 'file_already_exists',
        #                                    'already_in_batch', 'execute_strun', 'strun_executed'])

        #print('strun commands:')
        #print(mmm.strun)
        print(self.ssbcmdtable.t['index', 'real_input_file', 'ssbsteps', 'repeat_of_index_number',
                                 'index_contained_within', 'primary_strun', 'file_already_exists',
                                 'already_in_batch', 'execute_strun', 'strun_executed'])

        return(0)

    #@log_fail
    #@log_info
    def run_ssb_cmds(self, batchmode=False, ssblogFlag=None, ssberrorlogFlag=None):

        if self.onlyshow:
            self.logger.info('\n*** ONLYSHOW: skipping running the strun commands!\n')
            return(0)

        if ssblogFlag==None: ssblogFlag = self.cfg.params['output']['ssblogFlag']
        if ssberrorlogFlag==None: ssberrorlogFlag = self.cfg.params['output']['ssberrorlogFlag']

        # indeces of the commands that need to be run
        indeces2run, = np.where(self.ssbcmdtable.t['execute_strun'])

        # get rid of old output files
        self.ssbcmdtable.clean_old_output_files(indeces2run, outputfile_col='output_name')

        if batchmode:
            self.ssbcmdtable.run_cmds_batchmode(indeces2run, '%s.ssb_batch.txt' % self.basename,
                                                batchIDfilename='%s.ssb_batchID.txt' % self.basename,
                                                logFlag=ssblogFlag, errorlogFlag=ssberrorlogFlag,
                                                cmds_col='strun_command', cmds_executed_col='strun_executed')
        else:
            self.ssbcmdtable.run_cmds_serial(indeces2run,
                                             logFlag=ssblogFlag, errorlogFlag=ssberrorlogFlag,
                                             cmds_col='strun_command', cmds_executed_col='strun_executed')

        self.ssbcmdtable.check_if_files_exists(file_col='output_name', file_exists_col='file_exists')

        self.logger.debug('SSB Command Table:')
        self.logger.debug(self.ssbcmdtable.t['index', 'real_input_file', 'repeat_of_index_number',
                                             'index_contained_within', 'primary_strun', 'file_already_exists',
                                             'already_in_batch', 'execute_strun', 'strun_executed', 'file_exists'])

        # Save the ssb cmd table with the new columns
        self.ssbcmdtable.write('%s.ssbcmds.txt' % self.basename, verbose=True, clobber=True, exclude_names=self.ssbtable_excludecols4saving)

        return(0)

    #@log_fail
    #@log_info
    def mk_ref_cmds(self, force_redo_refcmds=False, maxNrefcmds=None):

        self.refcmdtable.t['refcmd'] = None
        self.logger.info('refcmd table colnames: {}'.format(self.refcmdtable.t.colnames))
        self.logger.debug(self.refcmdtable.t)

        # loop through refcomds table, and create the individual commands
        for i in range(len(self.refcmdtable.t)):
            reflabel = self.refcmdtable.t['reflabel'][i]
            self.logger.info('### cmd ID {}: building cmd for {}'.format(self.refcmdtable.t['cmdID'][i], reflabel))

            module_name = 'mkref_{}.py'.format(reflabel)
            module_and_path = [entry for entry in self.mkref_file_list if module_name in entry][0]
            refcmd = '{} {}.fits'.format(module_and_path, self.refcmdtable.t['outbasename'][i])

            # (1) get the input images from self.ssbcmdtable
            # (2) parse through the options

            # these are the input images for this reference file command
            indeces2run, = np.where(self.ssbcmdtable.t['cmdID'] == self.refcmdtable.t['cmdID'][i])
            # print the images and some error checking
            t_inputimages = len(set(self.inputimagestable.t['fitsfile']))
            if (t_inputimages != self.refcmdtable.t['Nim'][i]):
                self.logger.debug('Images found for this ref command:')
                self.logger.debug(self.ssbcmdtable.t['index', 'cmdID', 'reflabel', 'imlabel', 'imtype',
                                                 'imindex', 'imID', 'fitsfile'][indeces2run])

                # some error checking: this should never be true
                if (len(indeces2run)!=self.refcmdtable.t['Nim'][i]):
                    raise RuntimeError("Expected {} images, but {} found!".format(self.refcmdtable.t['Nim'][i], len(indeces2run)))

            inputlist_filename = '{}.imlist.txt'.format(self.refcmdtable.t['outbasename'][i])
            self.ssbcmdtable.write(inputlist_filename, indeces=indeces2run, verbose=self.verbose, clobber=True, exclude_names=self.ssbtable_excludecols4saving)
            refcmd += ' {}'.format(inputlist_filename)

            # now get the options specific to the mkref_X.py command
            parser4mkref = argparse.ArgumentParser(conflict_handler='resolve')
            mkref_reflabel = self.mkrefpackages["mkref_{}".format(reflabel)].mkrefclass()
            mkref_reflabel.default_optional_arguments(parser4mkref)
            mkref_reflabel.extra_optional_arguments(parser4mkref)
            # parse_known_args goes through all the arguments and
            # figures out which one are allowed for that particular
            # parser
            args_reflabel = parser4mkref.parse_known_args()

            # now get the allowed arguments and their values
            allowed_args_dict = vars(args_reflabel[0])
            # we want to have the sorted list of keys to this
            # dictionary, so that the commands always look the same!
            allowed_args = sorted(allowed_args_dict.keys())

            # loop through allowed options
            optionstring = ''
            for arg in allowed_args:
                # skip None options. These should be defaults!
                if allowed_args_dict[arg] is None:
                    continue

                # skip adding the cfgfile option if it is the same as the default!
                if arg == 'cfgfile':
                    if allowed_args_dict[arg] == parser4mkref.get_default('cfgfile'):
                        continue

                #the argument with -- is the key for the argparse dictionary!
                #argparse info: https://svn.python.org/projects/python/trunk/Lib/argparse.py
                arg_dashes = '--{}'.format(arg)
                if isinstance(parser4mkref._option_string_actions[arg_dashes], argparse._CountAction):
                    # 1st special case: counter
                    #print(arg,'CountAction',allowed_args_dict[arg],parser4mkref._option_string_actions[arg_dashes].nargs)
                    for n in range(allowed_args_dict[arg]):
                        optionstring+=' {}'.format(parser4mkref._option_string_actions[arg_dashes].option_strings[-1])
                elif isinstance(parser4mkref._option_string_actions[arg_dashes],argparse._AppendAction) and (parser4mkref._option_string_actions[arg_dashes].nargs==None):
                    # 2nd special case: append with nargs=None: argument is a simple list, not a list of lists if nargs!=None!
                    for n in range(len(allowed_args_dict[arg])):
                       optionstring+=' {} {}'.format(parser4mkref._option_string_actions[arg_dashes].option_strings[-1],allowed_args_dict[arg][n])
                else:
                    # there seem to be three cases for the normal argument values: single values, lists, and lists of lists
                    if isinstance(allowed_args_dict[arg], list):
                        if isinstance(allowed_args_dict[arg][0], list):
                            for arg2 in allowed_args_dict[arg]:
                                optionstring += ' {} {}'.format(parser4mkref._option_string_actions[arg_dashes].option_strings[-1],' '.join(arg2))
                        else:
                            optionstring += ' {} {}'.format(parser4mkref._option_string_actions[arg_dashes].option_strings[-1],' '.join(allowed_args_dict[arg]))
                    else:
                        optionstring += ' {} {}'.format(parser4mkref._option_string_actions[arg_dashes].option_strings[-1],allowed_args_dict[arg])

            refcmd += optionstring
            self.refcmdtable.t['refcmd'][i] = refcmd

            self.logger.info('refcmd is : {}'.format(refcmd))

        #if self.verbose>1:
        #    print(self.refcmdtable.t['reflabel', 'cmdID', 'refcmd'])

        # check if the reduced input files exist or not, and fill
        # 'file_already_exists' column with True or False
        self.refcmdtable.check_if_files_exists(file_col='outbasename', file_exists_col='file_already_exists', addsuffix='fits')

        # check if the reduced input files are already running in
        # batch (Condor), and set 'already_in_batch' to True or False
        self.refcmdtable.check_if_already_in_batch(cmd_col='refcmd')

        # pick the strun commands that need to be executed
        self.refcmdtable.pick_cmds_to_execute(execute_cmds=None, force_redo=force_redo_refcmds, maxNexe=maxNrefcmds, execute_cmds_col='execute_refcmd')

        # set strun_executed to None
        self.refcmdtable.t['refcmd_executed'] = np.full((len(self.refcmdtable.t)),None)

        self.refcmdtable.write('%s.refcmds.txt' % self.basename,verbose=True,clobber=True)

        self.logger.debug('ref table colnames:', self.refcmdtable.t.colnames)
        #self.logger.info(self.refcmdtable.t['reflabel', 'detector', 'cmdID', 'outbasename',
        #                                'file_already_exists', 'already_in_batch',
        #                                'execute_refcmd', 'refcmd_executed'])
        print(self.refcmdtable.t['reflabel', 'detector', 'cmdID', 'outbasename',
                                'file_already_exists', 'already_in_batch',
                                 'execute_refcmd', 'refcmd_executed'])

    #@log_fail
    #@log_info
    def run_ref_cmds(self, batchmode=False, reflogFlag=None, referrorlogFlag=None):

        if self.onlyshow:
            self.logger.error('\n*** ONLYSHOW: skipping running the ref commands!\n')
            return(0)

        if reflogFlag==None: reflogFlag = self.cfg.params['output']['reflogFlag']
        if referrorlogFlag==None: referrorlogFlag = self.cfg.params['output']['referrorlogFlag']

        # indeces of the commands that need to be run
        indeces2run, = np.where(self.refcmdtable.t['execute_refcmd'])

        # get rid of old output files
        self.refcmdtable.clean_old_output_files(indeces2run, outputfile_col='outbasename', addsuffix='fits')

        if batchmode:
            # execute in batch mode
            self.refcmdtable.run_cmds_batchmode(indeces2run, '%s.ref_batch.txt' % self.basename,
                                                batchIDfilename='%s.ref_batchID.txt' % self.basename,
                                                logFlag=reflogFlag, errorlogFlag=referrorlogFlag,
                                                cmds_col='refcmd', cmds_executed_col='refcmd_executed',
                                                outputfile_col='outbasename', addsuffix='fits')
        else:
            # execute serially
            self.refcmdtable.debug=False
            self.refcmdtable.run_cmds_serial(indeces2run,
                                             logFlag=reflogFlag, errorlogFlag=referrorlogFlag,
                                             cmds_col='refcmd', cmds_executed_col='refcmd_executed',
                                             outputfile_col='outbasename', addsuffix='fits')

        self.refcmdtable.check_if_files_exists(file_col='outbasename', file_exists_col='file_exists', addsuffix='fits')
        self.logger.debug(self.refcmdtable.t['reflabel', 'detector', 'cmdID', 'Nim', 'outbasename',
                                         'file_already_exists', 'already_in_batch', 'execute_refcmd',
                                         'refcmd_executed', 'file_exists'])

        # Save the ref cmd table with the new columns
        self.refcmdtable.write('%s.refcmds.txt' % self.basename, verbose=True, clobber=True)

        sys.exit(0)

        return(0)

    @log_fail
    @log_info
    def combinerefs(self):
        print("### combinerefs: NOT YET IMPLEMENTED!!!")
        sys.exit(0)

    @log_fail
    @log_info
    def overview(self):
        print("### overview: NOT YET IMPLEMENTED!!!")
        sys.exit(0)

    @log_fail
    @log_info
    def submitbatch(self):
        print("### submitbatch: NOT YET IMPLEMENTED!!!")
        sys.exit(0)


if __name__ == '__main__':

    mkrefs = mkrefsclass(counter=1)

    # Set up a list to hold early logging info. Each element in the
    # list will be a tuple composed of (level, message), where level
    # is one of 'warning', 'error', 'info', 'debug'. The output
    # directory is not defined until later, so we cannot create
    # the logfile yet.
    # mkrefs.loginfo = []

    parser = mkrefs.define_options()
    args = parser.parse_args()

    mkrefs.loginfo.append(('info', 'Input files: {}'.format(args.reflabels_and_imagelist)))

    # set verbose, debug, and onlyshow level
    mkrefs.verbose = args.verbose
    mkrefs.debug = args.debug
    mkrefs.onlyshow = args.onlyshow

    # Load config files
    mkrefs.cfg.loadcfgfiles(args.cfgfile,
                            extracfgfiles=args.extracfgfile,
                            params=args.params,
                            params4all=args.pall,
                            params4sections=args.pp,
                            verbose=args.verbose)


    # set the basedir, basename, ssbdir, outsubdir, runID, runlabel
    mkrefs.set_dirs(outrootdir=args.outrootdir, outsubdir=args.outsubdir,
                    runID=args.runID, runIDNdigits=args.runIDNdigits, newrunID=args.newrunID,
                    ssbdir=args.ssbdir, skip_runID_ssbdir=args.skip_runID_ssbdir)

    # Set up logging
    main_log_file = configure_logging("mkrefs", path=mkrefs.basedir, log_file_level=args.logFileLevel.lower(),
                                      log_screen_level=args.logScreenLevel.lower())
    print("LOG FILE: {}".format(main_log_file))
    mkrefs.logger = logging.getLogger('jwst_reffiles.mkrefs.mkrefsclass')
    mkrefs.logger.info('Log file configured.')

    # Add the earlier logged information to the logger
    for line in mkrefs.loginfo:
        level = line[0]
        if level == 'critical':
            mkrefs.logger.critical(line[1])
        if level == 'warning':
            mkrefs.logger.warning(line[1])
        elif level == 'error':
            mkrefs.logger.error(line[1])
        elif level == 'info':
            mkrefs.logger.info(line[1])
        elif level == 'debug':
            mkrefs.logger.debug(line[1])
        else:
            raise ValueError('Unrecognized logging level: {}'.format(level))


    # get the inputfile list and reflabel list. For input files, get into!
    mkrefs.organize_inputfiles(args.reflabels_and_imagelist)

    # check if current input files are consistent with previous ionput files of same output directory!
    mkrefs.check_inputfiles()

    # initialize the refcmds table. This table does not contain yet the refcmds...
    #mkrefs.mk_refcmds_table()

    # create the ssb commands. this also initializes the refcmds table
    mkrefs.mk_ssb_cmds(force_redo_strun=args.force_redo_strun, maxNstrun=args.maxNstrun)

    # run the ssb commands
    mkrefs.run_ssb_cmds(batchmode=args.batchmode)

    # create the reference file commands
    mkrefs.mk_ref_cmds(force_redo_refcmds=args.force_redo_ref, maxNrefcmds=args.maxNref)

    # run the reference file  commands
    mkrefs.run_ref_cmds(batchmode=args.batchmode)
    sys.exit(0)

    mkrefs.combinerefs()

    mkrefs.overview()

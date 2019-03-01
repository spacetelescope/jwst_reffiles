#!/usr/bin/env python
'''
wrapper around recarray with convenience functions to ease handling of tables
In addition: config file using yaml
A. Rest
'''

import os,io
import re
import sys
import types
import subprocess
import pytz

import astropy
from astropy.io import ascii
from astropy.time import Time
from datetime import datetime
import astropy.io.fits as fits
import numpy as np
from numpy import recarray
import numpy.lib.recfunctions as nprecf
import yaml


def makepath(path, raiseError=True):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError:
                raise RuntimeError('ERROR: Cannot create directory %s' % path)
            else:
                return(1)
    return(0)


def makepath4file(filename, raiseError=True):
    path = os.path.dirname(filename)
    if not os.path.isdir(path):
        return(makepath(path, raiseError=raiseError))
    else:
        return(0)

def save2file(filename,lines,verbose=0):
    if os.path.lexists(filename):
        os.remove(filename)
    if os.path.isfile(filename):
        raise RuntimeError('ERROR: Cannot remove %s' % filename)
    errorcode = append2file(filename,lines,verbose=verbose)
    return(errorcode)

def append2file(filename,lines,verbose=0):
    if filename==None:
        if verbose>2: print('No filename, returning')
        return(0)

    if isinstance(lines,str):
        lines = [lines,]
    if lines==None:
        return(0)
    if len(lines)==0:
        return(0)

    if os.path.isfile(filename):
        if verbose:print('Appending to file %s' % filename)
        buff = open(filename, 'a')
    else:
        if verbose:print('Writing to file %s' % filename)
        buff = open(filename, 'w')

    r=re.compile('\n$')
    #if first line does not have a \n, assume none of them have and add it!
    if not r.search(lines[0]):
       for i in range(len(lines)):
           lines[i]+='\n'

    buff.writelines(lines)
    buff.close()

    return(0)

def rmfile(filename,raiseError=1):
    " if file exists, remove it "
    if os.path.lexists(filename):
        os.remove(filename)
        if os.path.isfile(filename):
            if raiseError == 1:
                raise RuntimeError('Cannot remove {}'.format(filename))
            else:
                print('ERROR: Cannot remove {}'.format(filename))
                return(1)
    return(0)

def executecommand(cmd, successword, shell=True, errorlog=None, cmdlog=None, clobbercmdlog=False,
                   cmdlog_access_mode='w', errorlog_access_mode='a', cmdlog_save_as_chunk_flag=True,
                   log_buffering=-1, return_output=False,
                   verbose=2,
                   timezone='US/Eastern'):
    """ execute the given command cmd as a shell command.

    cmdlog and errorlog can be a filename or a file handle.

    clobbercmdlog=True and cmdlog_access_mode only apply if cmdlog is
    a filename. log_access_mode can be all allowed access_mode values
    for subprocess.Popen.  if clobbercmdlog=True, the cmdlog file is
    deleted before logging begins, even if cmdlog_access_mode='a'.

    if cmdlog!=None: As default, cmdlog_save_as_chunk_flag=True, which
    means that the output of the cmd is save in cmdlog *after* the
    command is finished. If there is an error, it could be that then
    the output is not saved in the cmdlog, thus for debugging, the
    user can set cmdlog_save_as_chunk_flag=False, and then each
    individual output line is immediately saved into cmdlog

    The output of the command is shown on screen in real-time if
    verbose>1.

    The errorlog captures the stderr output. By default, old errors
    are not deleted (i.e. errorlog_access_mode='a'). If there is an
    error, this routine returns 1

    if return_output=True, then the stdout and stderr are returned as
    lists: (returncode,stdout_lines,stderr_lines). Otherwise, only
    retruncode is returned

    """

    errorcode = 0

    time_cmdstart = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if verbose:
        print('{} executing: {}'.format(time_cmdstart, cmd))

    # check if logs are strings or file handles
    if cmdlog is not None:
        if isinstance(cmdlog, str):
            if clobbercmdlog:
                rmfile(cmdlog)
            cmdlog_filehandle = open(cmdlog, cmdlog_access_mode, log_buffering)
        elif isinstance(cmdlog, io.IOBase):
            if clobbercmdlog:
                raise RuntimeError('Cannot clobber a filehandle (cmdlog)')
            cmdlog_filehandle = cmdlog
        else:
            raise RuntimeError('Can\'t understand cmdlog')
        cmdlog_filehandle.write('{} executing: {}\n'.format(time_cmdstart, cmd))
    else:
        cmdlog_filehandle = None

    if errorlog is not None:
        if isinstance(errorlog, str):
            errorlog_filehandle = open(errorlog, errorlog_access_mode, log_buffering)
        elif isinstance(errorlog, io.IOBase):
            errorlog_filehandle = errorlog
        else:
            raise RuntimeError('Can\'t understand errorlog')
    else:
        errorlog_filehandle = None

    # start the subprocess...
    p = subprocess.Popen(cmd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)

    # capture the output while the subprocess is running
    if verbose > 1:
        print('Output:')
    stdout_lines = []

    # Use a regular expression to search for a word that indicate that the command executed successfully:
    # set successflag=0, and compile the search expression.
    if successword != '':
        successflag = 0
        m = re.compile(successword)
    else:
        successflag = 1
        m = None

    while p.poll() is None:
            l1 = p.stdout.readline()
            if l1 == b'':
                continue
            l = l1.decode('UTF-8').strip()
            if verbose > 1:
                print(l)
            if cmdlog_filehandle is not None or errorlog_filehandle is not None:
                if l != '':
                    if cmdlog_save_as_chunk_flag or errorlog_filehandle is not None or return_output:
                        stdout_lines.append(l+'\n')
                    if not cmdlog_save_as_chunk_flag:
                        #cmdlog_filehandle.write(l1)
                        cmdlog_filehandle.write(l+'\n')

            # check if the success word is in the output.
            if m is not None and successflag == 0:
                if m.search(l):
                    successflag = 1

    # Done!
    time_cmdend = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    endstring = 'execution finished ({} to {}): {}'.format(time_cmdstart, time_cmdend, cmd)
    if verbose:
        print(endstring)

    # if wanted, save the output into the cmdlog
    if cmdlog_filehandle is not None:
        if cmdlog_save_as_chunk_flag:
            cmdlog_filehandle.writelines(stdout_lines)
        cmdlog_filehandle.write(endstring+'\n')

    # check for errors!
    stderr_lines = p.stderr.readlines()

    # Coulnd't find success word? add to the error messages...
    if successflag != 1:
        errline = 'ERROR: could not find success expression {} in the output!'.format(successword)
        print(errline)
        stderr_lines.append('{}\n'.format(errline).encode())

    # if errors, print them again and save in errorlog
    if len(stderr_lines) > 0:

        errline = '\n* THERE WERE ERRORS IN THE EXECUTION OF CMD ({} to {}): {}'.format(time_cmdstart, time_cmdend, cmd)
        print(errline)
        if errorlog_filehandle is not None:
            #errorlog_filehandle.write('{}\n'.format(errline).encode())
            errorlog_filehandle.write(errline+'\n')
        print('* Error messages:')
        for l1 in stderr_lines:
            l = l1.decode('UTF-8').strip()
            print('*', l)
            if errorlog_filehandle is not None:
                errorlog_filehandle.write(l+'\n')
        errorcode = 1

    # close the files
    if cmdlog_filehandle is not None:
        # if the file handle was opened within routine, close it. if the filehandle was passed, don't close it
        if isinstance(cmdlog, str):
            cmdlog_filehandle.close()
    if errorlog_filehandle is not None:
        # if the file handle was opened within routine, close it. if the filehandle was passed, don't close it
        if isinstance(errorlog, str):
            errorlog_filehandle.close()

    if return_output:
        return(errorcode, stdout_lines, stderr_lines)
    else:
        return(errorcode)


class yamlcfgclass:
    def __init__(self):
        self.params = None
        self.errorflag = False

        self.sections = None
        self.envvarpattern = re.compile('\$(\w+)')

    def error(self, errormsg, raiseErrorFlag=True):
        self.errorflag = True
        if raiseErrorFlag:
            raise RuntimeError(errormsg)
        else:
            print('WARNING:', errormsg)
        return(0)

    def getsections(self):
        self.sections = []
        for section in self.params:
            if type(self.params[section]) is dict:
                self.sections.append(section)

    def subenvvarplaceholder(self, paramdict):
        """ Loop through all string parameters and substitute environment variables """
        for param in paramdict:
#            if type(paramdict[param]) is types.StringType:
            if isinstance(paramdict[param],str):
                envvarnames = self.envvarpattern.findall(paramdict[param])
                if envvarnames:
                    for name in envvarnames:
                        if not (name in os.environ):
                            raise RuntimeError("environment variable %s used in config file, but not set!" % name)
                        envval = os.environ[name]
                        subpattern = '\$%s' % (name)
                        paramdict[param] = re.sub(subpattern, envval, paramdict[param])
            elif isinstance(paramdict[param],dict):
                # recursive: sub environment variables down the dictiionary
                self.subenvvarplaceholder(paramdict[param])
        return(0)

    def addnewparams(self, newparams, requireParamExists=True):
        """ check if all the parameters in the newparams already
        exist. Throw a Runtime Error if requireParamExists=True and
        parameter does not exist. Then add them """

        errorflag = False
        for param in newparams:
            if not (param in self.params):
                self.error('parameter %s does not exist yet!' % (param), raiseErrorFlag=requireParamExists)
                self.params[param] = newparams[param]
                continue

            if type(newparams[param]) is dict:
                # sanity test
                if not (type(self.params[param])) is dict:
                    errorflag = self.error(('parameter %s is a section in the new parameters, '
                                            'but a value in the previous parameters! SKIPPING') %
                                           (param), raiseErrorFlag=requireParamExists)
                    print('SKIPPING %s!!' % param)
                    continue

                # we need to go through each parameter instead of just adding the new
                # dictionary in order to avoid deleting entries in that section that
                # are not getting redefined!!
                for param1 in newparams[param]:
                    if not (param1 in self.params[param]):
                        errorflag = self.error('parameter %s in section %s does not exist yet!' %
                                               (param1, param), raiseErrorFlag=requireParamExists)
                    self.params[param][param1] = newparams[param][param1]
            else:
                # sanity test
                if type(self.params[param]) is dict:
                    errorflag = self.error(('parameter %s is a value in the new parameters, but a '
                                            'section in the previous parameters! ') % (param),
                                           raiseErrorFlag=requireParamExists)
                    print('SKIPPING %s!!' % param)
                    continue

                self.params[param] = newparams[param]

        return(errorflag)

    def setval_all(self, param, val, allflag=True, requireParamExists=True):
        foundflag = False

        # contains the string to be loaded. We use the yaml.load(string) function so that
        # integers and floats are correctly converted!
        paramstring = ''

        # Check if it is a global parameter, but not a section
        if param in self.params:
            if type(self.params[param]) is dict:
                self.error("parameter %s is a section, thus cannot set it to a value" %
                           (param), raiseErrorFlag=requireParamExists)
            else:
                paramstring += '%s: %s\n' % (param, val)
                foundflag = True

        # go through the sections, and set if exists
        if allflag or (not foundflag):
            for section in self.sections:
                if param in self.params[section]:
                    paramstring += '%s:\n    %s: %s\n' % (section, param, val)
                    foundflag = True
                    if not allflag:
                        break

        if (not foundflag):
            self.error("Could not find parameter=%s" % (param), raiseErrorFlag=requireParamExists)
            return(1)
        else:
            # get the new parameters
            newparams = yaml.load(paramstring)
            # add them
            errorflag = self.addnewparams(newparams, requireParamExists=requireParamExists)
            return(errorflag)

    def setvals_all(self, param_val_list, allflag=True, requireParamExists=True):
        if param_val_list is None:
            return(0)
        errorflag = False
        for (param, val) in param_val_list:
            errorflag |= self.setval_all(param, val, allflag=allflag, requireParamExists=requireParamExists)
        return(errorflag)

    def setvals_section(self, section_param_val_list, requireParamExists=True):

        if section_param_val_list is None or len(section_param_val_list) == 0:
            return(0)

        errorflag = False

        # add all section,param,val tuples to the parameter string
        for (section, param, val) in section_param_val_list:
            paramstring = '%s:\n    %s: %s\n' % (section, param, val)

            # get the new parameters
            newparams = yaml.load(paramstring)

            # check and add them
            errorflag |= self.addnewparams(newparams, requireParamExists=requireParamExists)

        # make sure the parameters already exist
        # errorflag = self.checknewparams(newparams,requireParamExists=requireParamExists)

        if errorflag:
            print('!!! There are errors when setting the config parameters !!!!!')

        return(errorflag)

    def setvals(self, param_val_list, requireParamExists=True):

        if param_val_list is None or len(param_val_list) == 0:
            return(0)

        errorflag = False

        # add all section,param,val tuples to the parameter string
        for (param, val) in param_val_list:
            paramstring = '%s: %s\n' % (param, val)

            # get the new parameters
            newparams = yaml.load(paramstring)

            # check and add them
            errorflag |= self.addnewparams(newparams, requireParamExists=requireParamExists)

        # make sure the parameters already exist
        # errorflag = self.checknewparams(newparams,requireParamExists=requireParamExists)

        if errorflag:
            print('!!! There are errors when setting the config parameters !!!!!')

        return(errorflag)

    def addcfgfile(self, filename, requireParamExists=True, verbose=0):
        if self.params is None:
            raise RuntimeError("Cannot add extra cfg file %s before the main cfg file is loaded" % (filename))

        if verbose > 0:
            print('loading extra cfg file:', filename)

        # get the new parameters
        newparams = yaml.load(open(filename, 'r'))

        # check and add them
        errorflag = self.addnewparams(newparams, requireParamExists=requireParamExists)

        if errorflag:
            print('!!! There are errors in the config file %s !!!!!' % filename)
        return(errorflag)

    def loadcfgfiles(self, maincfgfile, extracfgfiles=None, params=None, params4all=None,
                     params4sections=None, requireParamExists=True, verbose=0):
        if maincfgfile is None:
            raise RuntimeError("Main config file is not specified!")

        # Read the main config file
        if verbose > 0:
            print('loading main cfg file:', maincfgfile)
        self.params = yaml.load(open(maincfgfile, 'r'))

        # identify all the sections in the config file
        self.getsections()

        # read in the additional config files
        if extracfgfiles is not None:
            if isinstance(extracfgfiles,str):
                self.addcfgfile(extracfgfiles, requireParamExists=requireParamExists, verbose=verbose)
            elif isinstance(extracfgfiles,list):
                for filename in extracfgfiles:
                    self.addcfgfile(filename, requireParamExists=requireParamExists, verbose=verbose)
            else:
                print('ERROR: this is the extra cfg filelist:', extracfgfiles)
                raise RuntimeError("ERROR: Don't know what to do with this filelist!!")

        # change the configs based on -p and --pp options
        self.setvals(params, requireParamExists=requireParamExists)
        self.setvals_all(params4all, requireParamExists=requireParamExists)
        self.setvals_section(params4sections, requireParamExists=requireParamExists)

        # replace environment variables in the configs
        self.subenvvarplaceholder(self.params)

        if self.errorflag:
            print('!!! There were errors when getting the config parameters!! !!!!!')

        if verbose > 4:
            print(yaml.dump(self.params))

        return(self.errorflag)


# http://docs.astropy.org/en/stable/table/index.html#astropy-table
# table groups!! http://docs.astropy.org/en/stable/table/operations.html#table-operations
class astrotableclass:
    def __init__(self, **kwargs):
        self.t = astropy.table.Table(**kwargs)
        self.verbose = 0

    def load(self, filename, namesMapping={}, formatMapping={}, **kwargs):
        #self.t = ascii.read(filename,format='commented_header',delimiter='\s',fill_values=[('-',0),('--',0)])
        try:
            self.t = ascii.read(filename, delimiter='\s')
        except:
            print('ERROR: could not read %s!' % filename)
            return(1)

        if len(self.t.colnames) < 1:
            print('ERROR: no data in %s?' % filename)
            return(1)

        if self.t.colnames[0] == 'col1':
            if self.verbose:
                print('WARNING: it looks like there is no header info!')

        if len(namesMapping) > 0:
            for name in self.t.colnames:
                if name in namesMapping:
                    self.t.rename_column(name, namesMapping[name])
        if len(formatMapping) > 0:
            for name in formatMapping:
                if name in self.t.colnames:
                    self.t[name].forma = formatMapping[name]
                else:
                    print('WARNING! col %s does not exist, so cannot format it!' % name)

    def write(self, filename, indeces=None, clobber=True, verbose=False, format='commented_header',
              makepathFlag=True, **kwargs):
        if verbose:
            print('Saving %s' % filename)

        # make the path to the file if necessary
        if makepathFlag:
            if makepath4file(filename, raiseError=False):
                print('ERROR: could not make directory for %s' % filename)
                return(1)

        # if clobber, then remove the old file first
        if clobber and os.path.isfile(filename):
            os.remove(filename)
            if os.path.isfile(filename):
                print('ERROR: could not save %s' % filename)
                return(2)

        if indeces is None:
            ascii.write(self.t, filename, format=format, **kwargs)
        else:
            ascii.write(self.t[indeces], filename, format=format, **kwargs)

        return(0)


    def fitsheader2table(self, fitsfilecolname, rowindices=None, requiredfitskeys=None,
                         optionalfitskey=None, raiseError=True, skipcolname=None, headercol=None):
        if rowindices is None:
            rowindices = range(len(self.t))

        if len(rowindices) == 0:
            print('no files!')
            return(0)

        # initialize columns if necessary
        if requiredfitskeys is not None:
            for fitskey in requiredfitskeys:
                if not (fitskey in self.t.colnames):
                    self.t[fitskey] = None
        if optionalfitskey is not None:
            for fitskey in optionalfitskey:
                if not (fitskey in self.t.colnames):
                    self.t[fitskey] = None

        if headercol is not None and (not (headercol in self.t.colnames)):
            self.t[headercol] = None

        # loop through the images
        for rowindex in rowindices:
            header = fits.getheader(self.t[fitsfilecolname][rowindex])
            if headercol is not None:
                self.t[headercol] = header

            if requiredfitskeys is not None:
                for fitskey in requiredfitskeys:
                    if fitskey in header:
                        self.t[fitskey][rowindex] = header[fitskey]
                    else:
                        self.t[fitskey][rowindex] = None
                        if raiseError:
                            raise RuntimeError("fits key %s does not exist in file %s" %
                                               (fitskey, self.t[fitsfilecolname][rowindex]))
                        else:
                            if skipcolname is not None:
                                self.t[skipcolname][rowindex] = True
            if optionalfitskey is not None:
                for fitskey in optionalfitskey:
                    if fitskey in header:
                        self.t[fitskey][rowindex] = header[fitskey]
                    else:
                        self.t[fitskey][rowindex] = None

    def dateobs2mjd(self, dateobscol, mjdcol, mjdobscol=None, timeobscol=None):
        if not (mjdcol in self.t.colnames):
            self.t[mjdcol] = None

        #print(self.t[dateobscol])
        #print(self.t[timeobscol])

        if timeobscol is not None:
            dateobslist = list(self.t[dateobscol]+'T'+self.t[timeobscol])
        else:
            dateobslist = list(self.t[dateobscol])

        dateobjects = Time(dateobslist,  format='isot', scale='utc')
        mjds = dateobjects.mjd

        self.t[mjdcol] = mjds

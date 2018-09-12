#!/usr/bin/env python
'''
wrapper around recarray with convenience functions to ease handling of tables
In addition: config file using yaml
A. Rest
'''
import sys,os,re,types
import numpy as np
from numpy import recarray
import numpy.lib.recfunctions as nprecf

from astropy.io import ascii
from astropy.time import Time
import astropy.io.fits as fits
import astropy
import yaml

def makepath(path,raiseError=True):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError:
                raise RuntimeError, 'ERROR: Cannot create directory %s' % path
            else:
                return(1)
    return(0)

def makepath4file(filename,raiseError=True):
    path = os.path.dirname(filename)
    if not os.path.isdir(path):
        return(makepath(path,raiseError=raiseError))
    else:
        return(0)

class yamlcfgclass:
    def __init__(self):
        self.params = None
        self.errorflag=False

        self.sections = None
        self.envvarpattern=re.compile('\$(\w+)')
        
    def error(self,errormsg,raiseErrorFlag=True):
        self.errorflag=True
        if raiseErrorFlag:
            raise RuntimeError,errormsg
        else:
            print 'WARNING:',errormsg
        return(0)
        
    def getsections(self):
        self.sections = []
        for section in self.params:
            if type(self.params[section]) is types.DictType:
                self.sections.append(section)

    def subenvvarplaceholder(self,dict):
        """ Loop through all string parameters and substitute environment variables """

        for param in dict:
            if type(dict[param]) is types.StringType:
                envvarnames=self.envvarpattern.findall(dict[param])
                if envvarnames:
                    for name in envvarnames:
                        if not (name in os.environ): 
                            raise RuntimeError,"environment variable %s used in config file, but not set!" % name
                        envval=os.environ[name]
                        subpattern='\$%s' % (name)
                        dict[param] = re.sub(subpattern,envval,dict[param])
            elif type(dict[param]) is types.DictType:
                # recursive: sub environment variables down the dictiionary 
                self.subenvvarplaceholder(dict[param])
        return(0)
                
    def addnewparams(self,newparams,requireParamExists=True):
        """ check if all the parameters in the newparams already
        exist. Throw a Runtime Error if requireParamExists=True and
        parameter does not exist. Then add them """

        errorflag=False
        for param in newparams:
            if not (param in self.params):
                self.error('parameter %s does not exist yet!' % (param),raiseErrorFlag=requireParamExists)
                self.params[param]=newparams[param]
                continue

            if type(newparams[param]) is types.DictType:
                # sanity test
                if not (type(self.params[param])) is types.DictType:
                    errorflag=self.error('parameter %s is a section in the new parameters, but a value in the previous parameters! SKIPPING' % (param),raiseErrorFlag=requireParamExists)
                    print 'SKIPPING %s!!' % param
                    continue

                # we need to go through each parameter instead of just adding the new dictionary in order to avoid deleting entries in that section that are not getting redefined!!
                for param1 in newparams[param]:
                    if not (param1 in self.params[param]):
                        errorflag=self.error('parameter %s in section %s does not exist yet!' % (param1,param),raiseErrorFlag=requireParamExists)
                    self.params[param][param1]=newparams[param][param1]
            else:
                # sanity test
                if type(self.params[param]) is types.DictType:
                    errorflag=self.error('parameter %s is a value in the new parameters, but a section in the previous parameters! ' % (param),raiseErrorFlag=requireParamExists)
                    print 'SKIPPING %s!!' % param
                    continue

                self.params[param]=newparams[param]
                
        return(errorflag)
    
    def setval_all(self,param,val,allflag=True,requireParamExists=True):
        foundflag=False

        # contains the string to be loaded. We use the yaml.load(string) function so that integers and floats are correctly converted!
        paramstring = ''
        
        # Check if it is a global parameter, but not a section
        if param in self.params:
            if type(self.params[param]) is types.DictType:
                self.error("parameter %s is a section, thus cannot set it to a value" % (param),raiseErrorFlag=requireParamExists)
            else:
                paramstring+='%s: %s\n' % (param,val)
                foundflag=True
                
        # go through the sections, and set if exists
        if allflag or (not foundflag):
            for section in self.sections:
                if param in self.params[section]:
                    paramstring+='%s:\n    %s: %s\n' % (section,param,val)
                    foundflag=True
                    if not allflag:
                        break
        
        if (not foundflag):
            self.error("Could not find parameter=%s" % (param),raiseErrorFlag=requireParamExists)
            return(1)
        else:
            # get the new parameters
            newparams = yaml.load(paramstring)
            # add them
            errorflag = self.addnewparams(newparams,requireParamExists=requireParamExists)
            return(errorflag)

    def setvals_all(self,param_val_list,allflag=True,requireParamExists=True):
        if param_val_list == None: return(0)
        errorflag = False
        for (param,val) in param_val_list:
            errorflag |= self.setval_all(param,val,allflag=allflag,requireParamExists=requireParamExists)
        return(errorflag)

    def setvals_section(self,section_param_val_list,requireParamExists=True):

        if section_param_val_list==None or len(section_param_val_list)==0:
            return(0)

        errorflag=False
        
        # add all section,param,val tuples to the parameter string
        for (section,param,val) in section_param_val_list:
            paramstring = '%s:\n    %s: %s\n' % (section,param,val)

            # get the new parameters
            newparams = yaml.load(paramstring)

            # check and add them
            errorflag |= self.addnewparams(newparams,requireParamExists=requireParamExists)

        # make sure the parameters already exist
        # errorflag = self.checknewparams(newparams,requireParamExists=requireParamExists)

        if errorflag:
            print '!!! There are errors when setting the config parameters !!!!!'
            
        return(errorflag)
    
    def setvals(self,param_val_list,requireParamExists=True):

        if param_val_list==None or len(param_val_list)==0:
            return(0)

        errorflag=False
        
        # add all section,param,val tuples to the parameter string
        for (param,val) in param_val_list:
            paramstring = '%s: %s\n' % (param,val)

            # get the new parameters
            newparams = yaml.load(paramstring)

            # check and add them
            errorflag |= self.addnewparams(newparams,requireParamExists=requireParamExists)

        # make sure the parameters already exist
        # errorflag = self.checknewparams(newparams,requireParamExists=requireParamExists)

        if errorflag:
            print '!!! There are errors when setting the config parameters !!!!!'
            
        return(errorflag)
    
    def addcfgfile(self,filename,requireParamExists=True,verbose=0):
        if self.params == None:
            raise RuntimeError,"Cannot add extra cfg file %s before the main cfg file is loaded" % (filename)
        
        if verbose>0:
            print 'loading extra cfg file:',filename

        # get the new parameters
        newparams = yaml.load(open(filename,'r'))

        # check and add them
        errorflag = self.addnewparams(newparams,requireParamExists=requireParamExists)

        if errorflag:
            print '!!! There are errors in the config file %s !!!!!' % filename
        return(errorflag)

    def loadcfgfiles(self,maincfgfile,extracfgfiles=None,params=None,params4all=None,params4sections=None,requireParamExists=True,verbose=0):
        if maincfgfile==None:
            raise RuntimeError,"Main config file is not specified!"        

        # Read the main config file
        if verbose>0:
            print 'loading main cfg file:',maincfgfile
        self.params = yaml.load(open(maincfgfile,'r'))
        
        # identify all the sections in the config file
        self.getsections()

        # read in the additional config files
        if extracfgfiles!=None:
            if type(extracfgfiles) is types.StringType:
                self.addcfgfile(extracfgfiles,requireParamExists=requireParamExists,verbose=verbose)
            elif type(extracfgfiles) is types.ListType:
                for filename in extracfgfiles:
                    self.addcfgfile(filename,requireParamExists=requireParamExists,verbose=verbose)
            else:
                if raiseErrorFlag:
                    print 'ERROR: this is the extra cfg filelist:',extracfgfiles
                    raise RuntimeError,"ERROR: Don't know what to do with this filelist!!"
                return(1)

        # change the configs based on -p and --pp options
        self.setvals(params,requireParamExists=requireParamExists)
        self.setvals_all(params4all,requireParamExists=requireParamExists)
        self.setvals_section(params4sections,requireParamExists=requireParamExists)        

        # replace environment variables in the configs
        self.subenvvarplaceholder(self.params)
        
        if self.errorflag:
            print '!!! There were errors when getting the config parameters!! !!!!!'

        if verbose>4:
            print yaml.dump(self.params)
            
        return(self.errorflag)

    
#http://docs.astropy.org/en/stable/table/index.html#astropy-table
#table groups!! http://docs.astropy.org/en/stable/table/operations.html#table-operations
class astrotableclass:
    def __init__(self, **kwargs):
        self.t = astropy.table.Table(**kwargs)
        self.verbose = 0
        
    def load(self,filename,namesMapping={},formatMapping={},**kwargs):
        #self.t = ascii.read(filename,format='commented_header',delimiter='\s',fill_values=[('-',0),('--',0)])
        try:
            self.t = ascii.read(filename,delimiter='\s')
        except:
            print 'ERROR: could not read %s!' % filename
            return(1)
        
        if len(self.t.colnames)<1:
            print 'ERROR: no data in %s?' % filename
            return(1)

        if self.t.colnames[0]=='col1':
            if self.verbose: print 'WARNING: it looks like there is no header info!'
            
        if len(namesMapping)>0:
            for name in self.t.colnames:
                if name in namesMapping:
                    self.t.rename_column(name,namesMapping[name])
        if len(formatMapping)>0:
            for name in formatMapping:
                if name in self.t.colnames:
                    self.t[name].format=formatMapping[name]
                else:
                    print 'WARNING! col %s does not exist, so cannot format it!' % name

    def write(self,filename,clobber=True,verbose=False,format='commented_header',makepathFlag=True,**kwargs):
        if verbose:
            print 'Saving %s' % filename

        # make the path to the file if necessary
        if makepathFlag:
            if makepath4file(filename,raiseError=False):
                print 'ERROR: could not make directory for %s' % filename
                return(1)

        # if clobber, then remove the old file first
        if clobber and os.path.isfile(filename):
            os.remove(filename)
            if os.path.isfile(filename):
                print 'ERROR: could not save %s' % filename
                return(2)
            
        ascii.write(self.t,filename,format=format,**kwargs)

    def fitsheader2table(self,fitsfilecolname,rowindices=None,requiredfitskeys=None,optionalfitskey=None,raiseError=True,skipcolname=None,headercol=None):
        if rowindices==None:
            rowindices = xrange(len(self.t))

        if len(rowindices)==0:
            print 'no files!'
            return(0)

        # initialize columns if necessary
        if requiredfitskeys!=None:
            for fitskey in requiredfitskeys:
                if not (fitskey in self.t.colnames):
                    self.t[fitskey]=None
        if optionalfitskey!=None:
            for fitskey in optionalfitskey:
                if not (fitskey in self.t.colnames):
                    self.t[fitskey]=None

        if headercol!=None and (not (headercol in self.t.colnames)):
            self.t[headercol]=None

        # loop through the images
        for rowindex in rowindices:
            header = fits.getheader(self.t[fitsfilecolname][rowindex])
            if headercol!=None:
                self.t[headercol]=header
                
            if requiredfitskeys!=None:
                for fitskey in requiredfitskeys:
                    if fitskey in header:
                        self.t[fitskey][rowindex]=header[fitskey]
                    else:
                        self.t[fitskey][rowindex]=None
                        if raiseError:
                            raise RuntimeError,"fits key %s does not exist in file %s" % (fitskey,self.t[fitsfilecolname][rowindex])
                        else:
                            if skipcolname!=None:
                                 self.t[skipcolname][rowindex]=True
            if optionalfitskey!=None:
                for fitskey in optionalfitskey:
                    if fitskey in header:
                        self.t[fitskey][rowindex]=header[fitskey]
                    else:
                        self.t[fitskey][rowindex]=None

    def dateobs2mjd(self,dateobscol,mjdcol,mjdobscol=None,timeobscol=None):
        if not (mjdcol in self.t.colnames):
            self.t[mjdcol]=None

        print self.t[dateobscol]
        print self.t[timeobscol]
            
        if timeobscol!=None:
            dateobslist = list(self.t[dateobscol]+'T'+self.t[timeobscol])
        else:
            dateobslist = list(self.t[dateobscol])

        dateobjects = Time(dateobslist,  format='isot', scale='utc')
        mjds = dateobjects.mjd

        self.t[mjdcol]=mjds



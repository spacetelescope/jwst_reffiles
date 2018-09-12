#!/usr/bin/env python
'''                                                                                                                                                                                            
create reference files for JWST instruments
A. Rest
'''

import sys, os,re,types,glob
import scipy,argparse
import numpy as np
import astropy.io.fits as fits


#from jwst_lib.models import RampModel

# get the root dir of the code. This is needed only if the scripts are not installed as a module!
if 'JWST_MKREFS_SRCDIR' in os.environ:
    rootdir = os.environ['JWST_MKREFS_SRCDIR']
    sys.path.extend([rootdir,
                     '%s/gain' % rootdir,
                     '%s/badpix_map' % rootdir])
#else:
#    rootdir = os.path.dirname(os.path.realpath(__file__))

from tools import astrotableclass,yamlcfgclass


class mkrefclass:
    def __init__(self):
        self.cfg = None
        self.outbasename = None

        self.inputfiles = None
    
    def defaultrefoptions(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage)
        parser.add_argument('--verbose', '-v', action='count')
        parser.add_argument('-d','--debug',help="debug", action='count')

        parser.add_argument('-o','--outbasename', default=None,
                            help='full output basename, e.g. for photometry table. Supercedes all other output options, e.g. --outrootdir and --outsubdir (default=%(default)s)')
        parser.add_argument('--outrootdir', default=None,
                            help='output root directory. If not specified, then the directory of input fits file is used.  (default=%(default)s)')
        parser.add_argument('--outsubdir', default=None,
                            help='subdir added to the output root directory (default=%(default)s)')
        parser.add_argument('--addsuffix', default=None,
                            help='suffix added to the output basename (default=%(default)s)')
        parser.add_argument('--imdates2subdir',help="add firstdate_lastdate subdir of the input images used (default=%(default)s)",action="store_true",default=False)

        # options for config file
        if 'JWST_MKREFS_CONFIGFILE' in os.environ and os.environ['JWST_MKREFS_CONFIGFILE']!='':
            cfgfile = os.environ['JWST_MKREFS_CONFIGFILE']
        else:
            cfgfile = None
        parser.add_argument('-c','--cfgfile', default=cfgfile,help='main config file. (default=%(default)s)')
        parser.add_argument('-e','--extracfgfile', action='append', default=None, 
                            help='additional config file. These cfg files do not need to have all parameters. They overwrite the parameters in the main cfg file.')
        parser.add_argument('-p', '--params', action='append', default=None, nargs=2,
                            help='"param val": change parameter in config file (not in section, only main part) (default=%(default)s)')
        parser.add_argument('--pall', action='append', default=None, nargs=2,
                            help='"param val". change parameter in all sections of config file  (section independent) (default=%(default)s)')
        parser.add_argument('--pp', action='append', default=None, nargs=3,
                            help='"section param val". change parameters in given section of config file (default=%(default)s)')
        
        return parser

    def gainoptions(self, parser,only_optional_arguments=False):
        if not only_optional_arguments:
            parser.add_argument('dark1', 
                                help='specify the first dark exposure (default=%(default)s)')
            parser.add_argument('dark2', 
                                help='specify the second dark exposure (default=%(default)s)')
            parser.add_argument('flat1', 
                                help='specify the first flat exposure (default=%(default)s)')
            parser.add_argument('flat2', 
                                help='specify the second flat exposure (default=%(default)s)')
        #parser.add_argument('--dummyoption1', 
        #                    help='testoption1 gain')
        parser = self.defaultrefoptions(parser=parser)
        return parser

    def readnoiseoptions(self, parser,only_optional_arguments=False):
        if not only_optional_arguments:
            parser.add_argument('dark', 
                                help='specify the dark exposure (default=%(default)s)')
        parser.add_argument('--dummyoption1', 
                            help='testoption1 readnoise')
        parser.add_argument('--dummyoption2', 
                            help='testoption2 readnoise')
        parser = self.defaultrefoptions(parser=parser)
        return parser

    def bpmoptions(self, parser,only_optional_arguments=False):
        if not only_optional_arguments:
            parser.add_argument('dark', 
                                help='specify the dark exposure (default=%(default)s)')
        parser = self.defaultrefoptions(parser=parser)
        return parser

        
    def refoptions(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage)
        self.defaultrefoptions(parser=parser, usage=None)
        subparsers = parser.add_subparsers(title='options for individual reference files',dest='reftype')
                                           #description='valid subcommands',
                                           #help='additional help')
        subparserlist = {}

        subparserlist['gain']=subparsers.add_parser('gain')
        self.gainoptions(subparserlist['gain'])
        
        subparserlist['readnoise']=subparsers.add_parser('readnoise')
        self.readnoiseoptions(subparserlist['readnoise'])

        subparserlist['bpm']=subparsers.add_parser('bpm')
        self.readnoiseoptions(subparserlist['bpm'])

        return(parser,subparserlist)
        
    def refoptions4mkrefs(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler='resolve')
        self.defaultrefoptions(parser=parser, usage=None)

        self.gainoptions(parser,only_optional_arguments=True)
        self.readnoiseoptions(parser,only_optional_arguments=True)

        return(parser)

    def setoutbasename(self,regtype,imagelist,outbasename=None,outrootdir=None,outsubdir=None,addsuffix=None,imdates2subdir=False):
        if outbasename!=None:
            self.outbasename = outbasename
            return(outbasename)

        if len(imagelist)>1:
            print '### more than 1 input file in setoutbasename not yet implemented!!!! exiting .... ###'
            sys.exit(0)
        
        # if outdir is passed or specified, use it, if not use the directiroy of the input filename
        (imdir,imbasename)=os.path.split(os.path.abspath(imagelist[0]))
        if outrootdir == None:
            if self.cfg.params['output']['outrootdir']=='':
                outrootdir = imdir
            else:
                outrootdir = self.cfg.params['output']['outrootdir']

        outbasename = '%s' % outrootdir

        # outsubdir: test if outsubdir is passed or if it is specified in the config file ...
        if outsubdir == None and self.cfg.params['output']['outsubdir']!='':
            outsubdir = self.cfg.params['output']['outsubdir']
        # ... and if yes, add it
        if outsubdir!=None:
            outbasename +='/%s' % outsubdir

        outbasename +='/%s_%s' % (regtype,re.sub('\.fits$','',imbasename))

        # addsuffix: test if addsuffix is passed or if it is specified in the config file ...
        if addsuffix == None and self.cfg.params['output']['addsuffix']!='':
            addsuffix = self.cfg.params['output']['addsuffix']
        # ... and if yes, add it
        if addsuffix!=None:
            if not re.search('\.',addsuffix):outbasename +='.'
            outbasename +='%s' % addsuffix
        
        if self.verbose>1:
            print 'Output basenbame:',outbasename

        self.outbasename = outbasename
            
        return(outbasename)

    def loadcfgfiles(self,maincfgfile,extracfgfiles=None,params=None,params4all=None,params4sections=None,requireParamExists=True):
        if self.cfg == None:
            self.cfg = yamlcfgclass()
        if self.cfg.loadcfgfiles(maincfgfile,extracfgfiles=extracfgfiles,
                                 params=params,params4all=params4all,params4sections=params4sections,
                                 requireParamExists=requireParamExists,verbose=self.verbose):
            raise RuntimeError,"Something went wrong when loading config files!"
        return(0)
    
    def checkssbinputfiles(self,reftype,args):
        print reftype,args
        print "### checkssbinputfiles: not yet implemented!!"

    def call_algorithm(self,reftype,args):
        print "### call_algorithm: not yet implemented!!"

    def validation(self,reftype,args):
        print "### validation: not yet implemented!!"

    def mkwebpage(self,reftype):
        print "### mkwebpage: not yet implemented!!"

    def mkref(self,arglist,onlyinit=False):
        (parser,subparserlist) = self.refoptions()
        args = parser.parse_args(arglist)

        # set verbose level
        self.verbose = args.verbose
        self.debug = args.debug

        # Load config files
        self.loadcfgfiles(args.cfgfile,
                          extracfgfiles=args.extracfgfile,
                          params=args.params,
                          params4all=args.pall,
                          params4sections=args.pp)


        
        self.setoutbasename(args.reftype,[args.dark],outrootdir=args.outrootdir,outsubdir=args.outsubdir,addsuffix=args.addsuffix,imdates2subdir=args.imdates2subdir)

        if onlyinit:
            return(0)

        self.checkssbinputfiles(args.reftype,args)

        self.call_algorithm(args.reftype,args)

        self.validation(args.reftype,args)

        self.mkwebpage(args.reftype)
        
if __name__=='__main__':

    mkref=mkrefclass()
    mkref.mkref(sys.argv[1:])


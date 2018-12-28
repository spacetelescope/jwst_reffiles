#!/usr/bin/env python
'''
create reference files for JWST instruments
A. Rest
'''

import argparse
import glob
import os
import re
import sys
import types

import astropy.io.fits as fits
import numpy as np
import scipy

from jwst_reffiles.utils.tools import astrotableclass,yamlcfgclass


class mkrefclass:
    def __init__(self):
        self.cfg = None
        self.outbasename = None

        self.inputfiles = None

    def default_optional_arguments(self, parser=None, usage=None, conflict_handler=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

        parser.add_argument('--verbose', '-v', action='count')
        parser.add_argument('-d', '--debug', help="debug", action='count')


        # options for config file
        if 'JWST_MKREFS_CONFIGFILE' in os.environ and os.environ['JWST_MKREFS_CONFIGFILE'] != '':
            cfgfile = os.environ['JWST_MKREFS_CONFIGFILE']
        else:
            cfgfile = None
        parser.add_argument('-c', '--cfgfile', default=cfgfile, help='main config file. (default=%(default)s)')
        parser.add_argument('-e', '--extracfgfile', action='append', default=None,
                            help=('additional config file. These cfg files do not need to have all '
                                  'parameters. They overwrite the parameters in the main cfg file.'))
        parser.add_argument('-p', '--params', action='append', default=None, nargs=2,
                            help=('"param val": change parameter in config file (not in section, only '
                                  'main part) (default=%(default)s)'))
        parser.add_argument('--pall', action='append', default=None, nargs=2,
                            help=('"param val". change parameter in all sections of config file '
                                  '(section independent) (default=%(default)s)'))
        parser.add_argument('--pp', action='append', default=None, nargs=3,
                            help=('"section param val". change parameters in given section of '
                                  'config file (default=%(default)s)'))

        return parser

    def positional_arguments(self, parser):
        """ 
        The filename of the output reference file is the first argument.
        The filename of the input image list is the second argument
        """
        parser.add_argument("outputreffilename", nargs=1, help=("The filename of the output reference file"))
        parser.add_argument("imageslist_filename", nargs=1, help=("The filename of the input image list"))
        return(0)

    def extra_optional_arguments(self, parser, only_optional_arguments=False):
        """ add here the extra options for a given reference file maker """
        print("PLACEHOLDER for extraoptions")
        return(0)

    def allargs(self,  parser=None, subparsers=None, subparserlist=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler='resolve')

        self.default_optional_arguments(parser=parser)

        self.extra_optional_arguments(parser)
        self.inputfileoptions(parser)

        return(parser)

    def callagorithm(self,args):
        """ add here the call to your algorithm, passing the appropriate args parameters """
        
        # import myscript from myfile
        # myscript.myroutine(args.D1,args.D2, args.F1, args.F2, myparameter=args.whateveroption)

        
    def loadcfgfiles(self, maincfgfile, extracfgfiles=None, params=None, params4all=None,
                     params4sections=None, requireParamExists=True):
        if self.cfg is None:
            self.cfg = yamlcfgclass()
        if self.cfg.loadcfgfiles(maincfgfile, extracfgfiles=extracfgfiles,
                                 params=params, params4all=params4all, params4sections=params4sections,
                                 requireParamExists=requireParamExists, verbose=self.verbose):
            raise RuntimeError("Something went wrong when loading config files!")
        return(0)

    def mkref(self, arglist, onlyinit=False):
        (parser, subparser, subparserlist) = self.refoptions()
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



if __name__ == '__main__':
    mkref = mkrefclassX()
    (parser) = mkref.refoptions()
    args = parser.parse_args()
    mkref.callagorithm(args)

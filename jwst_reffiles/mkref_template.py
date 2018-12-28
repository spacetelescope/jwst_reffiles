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


class mkrefclass_template:
    def __init__(self):
        self.cfg = yamlcfgclass()

################################################################################
# The following functions (until the next section) should only be modified
# in exceptional circumstances in the child classes of mkref_X.py
################################################################################
    
        
    def positional_arguments(self, parser):
        """ 
        The filename of the output reference file is the first argument.
        The filename of the input image list is the second argument
        """
        parser.add_argument("outputreffilename", nargs=1, help=("The filename of the output reference file"))
        parser.add_argument("imageslist_filename", nargs=1, help=("The filename of the input image list"))
        return(0)

    def default_optional_arguments(self, parser):

        parser.add_argument('--verbose', '-v', action='count')
        parser.add_argument('-d', '--debug', help="debug", action='count')


        # options for config file
        if 'JWST_MKREFS_CONFIGFILE' in os.environ and os.environ['JWST_MKREFS_CONFIGFILE'] != '':
            cfgfile = os.environ['JWST_MKREFS_CONFIGFILE']
        else:
            cfgfile = None
        parser.add_argument('-c', '--cfgfile', default=cfgfile, help='main config file. (default=%(default)s)')
        parser.add_argument('-e', '--extracfgfile', default=None,
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

        return(0)

    def allargs(self,  parser=None, usage=None, conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

        self.positional_arguments(parser)

        self.default_optional_arguments(parser)

        self.extra_optional_arguments(parser)

        return(parser)

    def DELMEloadcfgfiles(self, maincfgfile, extracfgfiles=None, params=None, params4all=None,
                     params4sections=None, requireParamExists=True, verbose=0):
        if self.cfg is None:
            self.cfg = yamlcfgclass()
        if self.cfg.loadcfgfiles(maincfgfile, extracfgfiles=extracfgfiles,
                                 params=params, params4all=params4all, params4sections=params4sections,
                                 requireParamExists=requireParamExists, verbose=verbose):
            raise RuntimeError("Something went wrong when loading config files!")
        return(0)

################################################################################
# The functions below can be modified in the child classes of mkref_X.py
################################################################################

    def extra_optional_arguments(self, parser):
        """ add here the extra options for a given reference file maker """
        print("PLACEHOLDER for extraoptions")
        return(0)

    def callagorithm(self,args):
        """ add here the call to your algorithm, passing the appropriate args parameters """
        
        # import myscript from myfile
        # myscript.myroutine(args.D1,args.D2, args.F1, args.F2, myparameter=args.whateveroption)




if __name__ == '__main__':
    mkref = mkrefclass_template()
    parser = mkref.allargs()
    args = parser.parse_args()
    mkref.callagorithm(args)

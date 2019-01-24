#!/usr/bin/env python
'''
This class contains the basic infrastructure needed when plugging in a
reference file-generating module. This class and the functions within
it should not have to be changed by the user when wrapping a user-
written module. All functions which need to be modified will be in the
mkref_your_new_module_name.py file, which should be based on
plugin_template.py.
'''

import argparse
import glob
import os
import re
import sys
import types

from astropy.io import ascii
import astropy.io.fits as fits
import numpy as np

from jwst_reffiles.utils.tools import astrotableclass, yamlcfgclass


class mkrefclass_template:
    def __init__(self):
        self.cfg = yamlcfgclass()

        self.reflabel = None
        self.reftype = None

        self.args = None
        self.parser = None

        self.inputimages = None

################################################################################
# The following functions (until the next section) should only be modified
# in exceptional circumstances in the child classes of mkref_X.py
################################################################################

    def positional_arguments(self, parser):
        """Add input and output filenames to the parser.
        The filename of the output reference file is the first argument.
        The filename of the input image list is the second argument

        Parameters
        ----------
        parser : argparse.parser
        """
        parser.add_argument("outputreffilename", help=("The filename of the output reference file"))
        parser.add_argument("imageslist_filename", help=("The filename of the input image list"))
        return(0)

    def default_optional_arguments(self, parser):
        """Optional arguments common to all plugged in modules

        Parameters
        ----------
        parser : argparse.parser
        """
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

        return(0)

    def define_args(self, parser=None, usage=None, conflict_handler='resolve'):
        """Create argparser and populate with arguments

        Parameters
        ----------
        parser : argparse.parser

        usage : str
            Usage string for command line

        conflict_handler : str
            Method of resolving conflicting argument names. Input to argparse

        Returns
        -------
        parser : argparse.parser
        """
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

        self.positional_arguments(parser)

        self.default_optional_arguments(parser)

        self.extra_optional_arguments(parser)

        return(parser)

    def determine_parameters(self):
        """Loop over the parameters that are present in both the config file
        and the extra_optional_arguments function above. If the parameter
        is set in extra_optional_arguments, then use that value in the call
        to the reference file generator. If not, use the value in the config
        ile.
        """
        key_list = set(list(self.cfg.params[self.reflabel].keys()) + list(self.args.__dict__.keys()))
        self.parameters = {}
        for key in key_list:
            if key in self.cfg.params[self.reflabel]:
                self.parameters[key] = self.cfg.params[self.reflabel][key]
            if key in self.args.__dict__:
                if self.args.__dict__[key] is not None:
                    self.parameters[key] = self.args.__dict__[key]

    def initialize(self, argument_list=None):
        """Read in and define all arguments from the command line as well
        as the config file
        """
        # first, parse the arguments
        self.parser = self.define_args()
        self.args = self.parser.parse_args()

        # Load config files
        self.cfg.loadcfgfiles(self.args.cfgfile,
                              extracfgfiles=self.args.extracfgfile,
                              params=self.args.params,
                              params4all=self.args.pall,
                              params4sections=self.args.pp,
                              verbose=self.args.verbose)

        # load the input images
        self.inputimages = ascii.read(self.args.imageslist_filename, delimiter='\s')

    def wrap_it_up(self, reffilename=None):
        #todo, add testing
        # (1) load the reference file into the approprite data model
        # (2) run it through ssb and see if it breaks

        #if self.cfg.params[self.reflabel]['test4datamodel']:
        #    print('To be implemented: test for data models!!')

        #if self.cfg.params[self.reflabel]['test4ssb']:
        #    print('To be implemented: run it through ssb!!')

        return(0)

    def make_reference_file(self):

        self.initialize()

        self.determine_parameters()

        self.callalgorithm()

        self.wrap_it_up()

        print('{} successfully finished!'.format(self.reflabel))

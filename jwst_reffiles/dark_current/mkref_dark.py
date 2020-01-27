#! /usr/bin/env python

"""
Plug-in script for the readnoise reffile creation module:
readnoise.py, which uses the readnoise generation algorithms
decided upon by the JWST reference file generation working group.
This class is based on that in the
template file: jwst_reffiles/templates/plugin_template.py
"""

import argparse
import copy
import os
import re
import sys
import types

from jwst_reffiles.plugin_wrapper import mkrefclass_template

# import the dark current reference file script
from jwst_reffiles.dark_current import dark_reffile

class mkrefclass(mkrefclass_template):
    def __init__(self, *args, **kwargs):
        mkrefclass_template.__init__(self, *args, **kwargs)

        # Set the reflabel as the name of the imported module
        self.reflabel = 'dark'

        # Set the reftype
        self.reftype = 'dark'

    def extra_optional_arguments(self, parser):
        """Any arguments added here will give the option of overriding
        the default argument values present in the config file. To override,
        call these arguments from the command line in the call to mkrefs.py
        """

        parser.add_argument('--max_equivalent_groups', help=('The maximum number of groups to have open '
                                                             'at one time. This is to avoid running into '
                                                             'memory issues during execution.'))
        parser.add_argument('--sigma_threshold', help=('The number of sigma to use when calculating the iterative '
                                                       'sigma-clipped mean and stdev dark values.'))
        parser.add_argument('--output_file', help=('Name of fits file to use for the final reference file'))
        parser.add_argument('--contribution_file', help=('Name of PDF file to contain maps of the number of '
                                                         'values contributing to the final dark current value '
                                                         'for each pixel/group.'))
        parser.add_argument('--pedigree', help='Pedigree for the data used to create the reffile (e.g. GROUND)')
        parser.add_argument('--use_after', help='Use after date/time for the data used to create the reffile (e.g. 2037-09-30)')
        parser.add_argument('--author', help='Author of the reference file')
        parser.add_argument('--descrip', help='Descrip value for the reffile header')
        parser.add_argument('--history', help='Text to include as the HISTORY entry of the reference file.')
        return(0)

    def callalgorithm(self):
        """Call the dark current algorithm. The only requirement is that the output
        reference file is saved as self.args.outputreffilename

        mkrefs.py will supply the input files in self.inputimages['output_name'].
        This will be a list containing the filenames to use as input. The
        file types (e.g. dark, flat) associated with each filename are
        contained in self.inputimages['imtype']. From this, you can specify
        the appropriate file names in the call to your module.
        """

        # Call the wrapped module and provide the proper arguments from the
        # self.parameters dictionary.
        dark_reffile.dark(file_list=self.inputimages,
                          max_equivalent_groups=self.parameters['max_equivalent_groups'],
                          sigma_threshold=self.parameters['sigma_threshold'],
                          pedigree=self.parameters['pedigree'],
                          descrip=self.parameters['descrip'],
                          use_after=self.parameters['use_after'],
                          author=self.parameters['author'],
                          history=self.parameters['history'],
                          output_file=self.args.outputreffilename,
                          contribution_file=self.parameters['contribution_file'])
        return(0)


if __name__ == '__main__':
    """This should not need to be changed. This will read in the config
    files, import the script, generate the self.parameters dictionary, and
    run the argument parser above.
    """
    mkref = mkrefclass()
    mkref.make_reference_file()

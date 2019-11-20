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

# import the readnoise script
from jwst_reffiles.readnoise import readnoise

class mkrefclass(mkrefclass_template):
    def __init__(self, *args, **kwargs):
        mkrefclass_template.__init__(self, *args, **kwargs)

        # Set the reflabel as the name of the imported module
        self.reflabel = 'readnoise'

        # Set the reftype
        self.reftype = 'rdn'

    def extra_optional_arguments(self, parser):
        """Any arguments added here will give the option of overriding
        the default argument values present in the config file. To override,
        call these arguments from the command line in the call to mkrefs.py
        """

        parser.add_argument('--method', help=('The method to use when calculating the readnoise. '
                                              'Options are: stack and ramp'))
        parser.add_argument('--group_diff_type', help=('The method for calculating group differences. '
                                                       'Options are: independent and consecutive'))
        parser.add_argument('--clipping_sigma', help=('Number of sigma to use when sigma-clipping.'))
        parser.add_argument('--max_clipping_iters', help=('Maximum number of iterations to use when '
                                                          'sigma-clipping.'))
        parser.add_argument('--nproc', help=('The number of processes to use during multiprocessing. '))
        parser.add_argument('--slice_width', help=('The width (in pixels) of the image slice to use '
                                                   'during multiprocessing. The readnoise of each slice '
                                                   'is calculated separately during multiprocessing and '
                                                   'combined together at the end of processing. Only '
                                                   'relevant if method==stack.'))
        parser.add_argument('--author', help=('CRDS-required name of the reference file author, to be '
                                              'placed in the referece file header.'))
        parser.add_argument('--description', help=('CRDS-required description of the reference file, to '
                                                   'be placed in the reference file header.'))
        parser.add_argument('--pedigree', help=('CRDS-required pedigree of the data used to create the '
                                                'reference file.'))
        parser.add_argument('--useafter', help=('CRDS-required date of earliest data with which this '
                                                'reffile should be used. (e.g. 2019-04-01T00:00:00).'))
        parser.add_argument('--history', help=('CRDS-required history section to place in the reference '
                                               'file header.'))
        parser.add_argument('--subarray', help=('CRDS-required subarray for which to use this reference '
                                                'file for.'))
        parser.add_argument('--readpatt', help=('CRDS-required read pattern for which to use this '
                                                'reference file for.'))
        parser.add_argument('--save_tmp', help=('Option to save the final readnoise map before turning it '
                                                'into CRDS format. This is useful if the CRDS transformation '
                                                'fails; in this scenario, you wont lose all of the final '
                                                'readnoise results so all of the previous processing wasnt '
                                                'for nothing.'))

        return(0)

    def callalgorithm(self):
        """Call the readnoise algorithm. The only requirement is that the output
        reference file is saved as self.args.outputreffilename

        mkrefs.py will supply the input files in self.inputimages['output_name'].
        This will be a list containing the filenames to use as input. The
        file types (e.g. dark, flat) associated with each filename are
        contained in self.inputimages['imtype']. From this, you can specify
        the appropriate file names in the call to your module.
        """

        # Call the wrapped module and provide the proper arguments from the
        # self.parameters dictionary.
        readnoise.make_readnoise(self.inputimages,
                                 method=self.parameters['method'],
                                 group_diff_type=self.parameters['group_diff_type'], 
                                 clipping_sigma=self.parameters['clipping_sigma'], 
                                 max_clipping_iters=self.parameters['max_clipping_iters'], 
                                 nproc=self.parameters['nproc'], 
                                 slice_width=self.parameters['slice_width'], 
                                 outfile=self.args.outputreffilename, 
                                 author=self.parameters['author'], 
                                 description=self.parameters['description'], 
                                 pedigree=self.parameters['pedigree'], 
                                 useafter=self.parameters['useafter'], 
                                 history=self.parameters['history'], 
                                 subarray=self.parameters['subarray'], 
                                 readpatt=self.parameters['readpatt'], 
                                 save_tmp=self.parameters['save_tmp'])

        return(0)


if __name__ == '__main__':
    """This should not need to be changed. This will read in the config
    files, import the script, generate the self.parameters dictionary, and
    run the argument parser above.
    """
    mkref = mkrefclass()
    mkref.make_reference_file()

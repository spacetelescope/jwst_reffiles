#!/usr/bin/env python
'''
Plug-in script for a user-generated reference file creation module.

All arguments/keyword inputs to the reference-file creation module
must be included in either the config file or the extra_optional_arguments
function below, or both.

Inputs for the parameters in the call to your
        module can come from two different locations:

        1) the config file
        2) the argument parser

You may add any parameters for your module to the mkrefs config file.
The values for these parameters will then be read in from the config
file when mkrefs.py is called. If you wish, you may also add the
options to the argument parser in this module. This will allow you
to override the values set in the config file when calling mkrefs.py
via the command line.
'''

import argparse
import os
import re
import sys
import types

from jwst_reffiles.plugin_wrapper import mkrefclass_template

# import your script!
from jwst_reffiles.example.example_readnoise_module import MyReadnoise


class mkrefclass(mkrefclass_template):
    def __init__(self, *args, **kwargs):
        mkrefclass_template.__init__(self, *args, **kwargs)

        # You need to set this correctly
        self.reflabel = 'example_readnoise_module'
        self.reftype = 'rdnoise'

    def extra_optional_arguments(self, parser):
        """Any arguments added here will give the option of overriding
        the default argument values present in the config file.
        """
        #parser.add_argument('--verbose', help='Boolean controlling information printed to the screen')
        parser.add_argument('--boxsize',
                            help='Size in pixels of boxes in grid across detector')
        parser.add_argument('--sigma_threshold',
                            help='Sigma-clipping threshold to use when calculating stdev')
        return(0)

    def callalgorithm(self):
        """Call your algorithm. The only requirement is that the output
        reference file is saved as self.args.outputreffilename

        mkrefs.py will supply the input files in self.inputimages['name'].
        This will be a list containing the filenames to use as input. The
        file types (e.g. dark, flat) associated with each filename are
        contained in self.inputimages['imtype']. From this, you can specify
        the appropriate file names in the call to your module.
        """
        print('*** This is the output filename of the reference file: {}'
              .format(self.args.outputreffilename))
        print('*** This is the file with the list of the input images: {}'
              .format(self.args.imageslist_filename))
        print('*** The table in this file has been loaded into self.inputimages')

        print('*** These are the filenames of the input images, and some more info:')
        for i in range(len(self.inputimages)):
            print('image {} has the type {} and taken at MJD={}'.format(self.inputimages['output_name'][i],
                                                                        self.inputimages['imtype'][i],
                                                                        self.inputimages['MJD'][i]))

        # Loop over the parameters that are present in both the config file
        # and the extra_optional_arguments function above. If the parameter
        # is set in extra_optional_arguments, then use that value in the call
        # to the reference file generator. If not, use the value in the config
        # file.
        params = self.determine_parameters()

        # Call the wrapped module and provide the proper arguments
        ron = MyReadnoise(self.inputimages['output_name'].data[0],
                          boxsize=params['boxsize'], sigma_threshold=params['sigma_threshold'],
                          output_name=self.args.outputreffilename)

        return(0)


if __name__ == '__main__':
    """This should not need to be changed. This will read in the config
    files, import the script, and run the argument parser above.
    """
    mkref = mkrefclass()
    mkref.make_reference_file()

#!/usr/bin/env python
'''
Example plug-in script for the user-generated readnoise creation module:
example_readnoise_module.py. This class is based on that in the
template file: jwst_reffiles/templates/plugin_template.py
'''

import argparse
import os
import re
import sys
import types

from jwst_reffiles.plugin_wrapper import mkrefclass_template

# import the readnoise script
from jwst_reffiles.example.example_readnoise_module import MyReadnoise


class mkrefclass(mkrefclass_template):
    def __init__(self, *args, **kwargs):
        mkrefclass_template.__init__(self, *args, **kwargs)

        # Set the reflabel as the name of the imported module
        self.reflabel = 'example_readnoise_module'

        # Set the reftype as rdnoise since example_readnoise_modeul is a
        # readnoise reference file generator.
        self.reftype = 'rdnoise'

    def extra_optional_arguments(self, parser):
        """Any arguments added here will give the option of overriding
        the default argument values present in the config file. To override,
        call these arguments from teh command line in the call to mkrefs.py
        """
        parser.add_argument('--boxsize',
                            help='Size in pixels of boxes in grid across detector')
        parser.add_argument('--sigma_threshold',
                            help='Sigma-clipping threshold to use when calculating stdev')
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

        # Call the wrapped module and provide the proper arguments from the
        # self.parameters dictionary.
        ron = MyReadnoise(self.inputimages['output_name'].data[0],
                          boxsize=self.parameters['boxsize'],
                          sigma_threshold=self.parameters['sigma_threshold'],
                          output_name=self.args.outputreffilename)
        return(0)


if __name__ == '__main__':
    """This should not need to be changed. This will read in the config
    files, import the script, generate the self.parameters dictionary, and
    run the argument parser above.
    """
    mkref = mkrefclass()
    mkref.make_reference_file()

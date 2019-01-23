#!/usr/bin/env python
'''
Example plug-in script for readnoise
A. Rest
'''

import argparse
import os
import re
import sys
import types

from plugin_wrapper import mkrefclass_template

# import your script!
from jwst_reffiles.example.example_readnoise_module import MyReadnoise


class mkrefclass(mkrefclass_template):
    def __init__(self, *args, **kwargs):
        mkrefclass_template.__init__(self, *args, **kwargs)

        # You need to set this correctly
        self.reflabel = 'example_readnoise_module_to_plug_in'
        self.reftype = 'rdnoise'

    def extra_optional_arguments(self, parser):
        # you can add whatever optional arguments here!
        parser.add_argument('--verbose', help='Boolean controlling information printed to the screen')
        parser.add_argument('--boxsize',
                            help='Size in pixels of boxes in grid across detector')
        parser.add_argument('--sigma_threshold',
                            help='Sigma-clipping threshold to use when calculating stdev')
        parser.add_argument('--output_directory', help="Directory into which output file is saved")
        return(0)

    def callalgorithm(self):
        # call your algorithm. The only requirement is that the output
        # reference file is saved as self.args.outputreffilename

        print('*** This is the output filename of the reference file: {}'.format(self.args.outputreffilename))
        print('*** This is the file with the list of the input images: {}'.format(self.args.imageslist_filename))
        print('*** The table in this file got loaded into self.inputimages')

        print('*** These are the filenames of the input images, and some more info:')
        for i in range(len(self.inputimages)):
            print('image {} has the type {} and taken at MJD={}'.format(self.inputimages['output_name'][i],
                                                                        self.inputimages['imtype'][i],
                                                                        self.inputimages['MJD'][i]))

        print('*** The config file is also loaded. This is how you can access it. There can be a section for each reflabel.')
        print('*** section {}, parameter {}: {}'.format(self.reflabel,'random_option',self.cfg.params[self.reflabel]['random_option']))

        print('*** You can use all this info to call whatever algorithm you like to apply!')
        ron = MyReadnoise(self.inputimages['output_name'], verbose=self.verbose, boxsize=self.boxsize,
                          sigma_threshold=self.sigma_threshold, output_directory=self.output_directory)

        so wait...the input options here should be from parser here, or from config file???

            self.inputimages['output_name'][0],self.inputimages['output_name'][1],
                                self.inputimages['output_name'][0],self.inputimages['output_name'][1],
                                whateveroption1=True,whateveroption2=self.cfg.params[self.reflabel]['random_option'])

        print('*** Now save it to the output reference filename {}'.format(self.args.outputreffilename))
        mymask.save(self.args.outputreffilename)

        return(0)


if __name__ == '__main__':
    mkref = mkrefclass()
    # mkref.make_reference_file()
    mkref.extra_optional_arguments()
    mkref.callalgorithm()

#!/usr/bin/env python
'''
Plug-in script for the bad pixel mask creation module:
bad_pixel_mask.py, which uses the bad pxiel mask geneartion algorithms
decided upon by the JWST reference file generation working group.
This class is based on that in the
template file: jwst_reffiles/templates/plugin_template.py
'''

import argparse
import copy
import os
import re
import sys
import types

from jwst_reffiles.plugin_wrapper import mkrefclass_template

# import the bad pixel mask script
from jwst_reffiles.bad_pixel_mask import bad_pixel_mask as bpm
from jwst_reffiles.utils.definitions import PIPE_STEPS


class mkrefclass(mkrefclass_template):
    def __init__(self, *args, **kwargs):
        mkrefclass_template.__init__(self, *args, **kwargs)

        # Set the reflabel as the name of the imported module
        self.reflabel = 'bad_pixel_mask'

        # Set the reftype
        self.reftype = 'bpm'

    def extra_optional_arguments(self, parser):
        """Any arguments added here will give the option of overriding
        the default argument values present in the config file. To override,
        call these arguments from the command line in the call to mkrefs.py
        """
        parser.add_argument('--dead_search', help='Whether or not to search for DEAD pixels')

        parser.add_argument('--low_qe_and_open_search', help=('Whether or not to search for LOW_QE, OPEN, '
                                                            'and ADJ_OPEN pixels'))
        parser.add_argument('--dead_search_type', help=('Type of search to use when looking for dead pixels. '
                                                      'Options are: sigma_rate, absolute_rate, and '
                                                      'zero_signal'))
        parser.add_argument('--sigma_threshold', help=('Number of standard deviations to use when sigma-'
                                                     'clipping to calculate the mean slope image or the mean '
                                                     'across the detector'))
        parser.add_argument('--normalization_method', help=('Specify how the mean image is normalized prior '
                                                            'to searching for bad pixels.'))
        parser.add_argument('--smoothing_box_width', help=('Width in pixels of the box kernel to use to '
                                                         'compute the smoothed mean image'))
        parser.add_argument('--smoothing_type', help='Type of smoothing to do ``Box2D `` or ``median`` filtering')
        parser.add_argument('--dead_sigma_threshold', help=('Number of standard deviations below the mean at '
                                                          'which a pixel is considered dead.'))
        parser.add_argument('--max_dead_norm_signal', help=('Maximum normalized signal rate of a pixel that is '
                                                          'considered dead'))
        parser.add_argument('--run_dead_flux_check', help=('Whether or not to check for dead pixels using an absolute flux value'))
        parser.add_argument('--dead_flux_check_files', nargs='+', help=('List of ramp (uncalibrated) files to use to check the '
                                                                        'flux of average of last 4 groups. If None then the '
                                                                        'ramp files are not read in and no flux_check is done.'))
        parser.add_argument('--flux_check', type=int, help=('Tolerance on average signal in last 4 groups. If dead_flux_check is '
                                                  'a list of uncalibrated files, then the average of the last four groups '
                                                  'for all the integrations is determined. If this average > flux_check '
                                                  'then this pixel is not a dead pixel.'))
        parser.add_argument('--max_low_qe_norm_signal', help=('The maximum normalized signal a pixel can have '
                                                            'and be considered low QE.'))
        parser.add_argument('--max_open_adj_norm_signal', help=('The maximum normalized signal a pixel '
                                                              'adjacent to a low QE pixel can have in order '
                                                              'for the low QE pixel to be reclassified as '
                                                              'OPEN'))
        parser.add_argument('--manual_flag_file', help=(('Name of file containing list of pixels to be added manually')))
        parser.add_argument('--do_not_use', help=('List of bad pixel types where the DO_NOT_USE flag should '
                                                  'also be applied (e.g. ["DEAD", "LOW_QE"])'))
        parser.add_argument('--author', help=('CRDS-required name of the reference file author, to be placed '
                                            'in the referece file header'))
        parser.add_argument('--description', help=('CRDS-required description of the reference file, to be '
                                                 'placed in the reference file header'))
        parser.add_argument('--pedigree', help=('CRDS-required pedigree of the data used to create the '
                                              'reference file'))
        parser.add_argument('--useafter', help=('CRDS-required date of earliest data with which this reference '
                                              'file should be used. (e.g. "2019-04-01 00:00:00"'))
        parser.add_argument('--history', help='Text to be placed in the HISTORY keyword of the output reference file')
        parser.add_argument('--quality_check', help=("If True, the pipeline is run using the output reference "
                                                   "file to be sure the pipeline doens't crash"))
        return(0)

    def callalgorithm(self):
        """Call the bpm algorithm. The only requirement is that the output
        reference file is saved as self.args.outputreffilename

        mkrefs.py will supply the input files in self.inputimages['output_name'].
        This will be a list containing the filenames to use as input. The
        file types (e.g. dark, flat) associated with each filename are
        contained in self.inputimages['imtype']. From this, you can specify
        the appropriate file names in the call to your module.
        """
        # Call the wrapped module and provide the proper arguments from the
        # self.parameters dictionary.

        # A set of uncal files is needed in the case that ``run_dead_flux_check``
        # is True. If dead_flux_check_files is None or 'none' then a list of uncal
        # files is generated from the input list of rate files.
        if self.parameters['run_dead_flux_check']:
            if self.parameters['dead_flux_check_files'] is None:
                self.parameters['dead_flux_check_files'] = 'none'

            input_file_directory, input_file_name = os.path.split(self.inputimagetable['fitsfile'][0])

            if isinstance(self.parameters['dead_flux_check_files'], str):
                possible_suffixes = copy.deepcopy(PIPE_STEPS)
                possible_suffixes.extend(['0_ramp_fit', '1_ramp_fit'])
                uncal_files = []
                for filename in self.inputimages:
                    directory, file = os.path.split(filename)
                    for suffix in possible_suffixes:
                        if suffix in file:
                            file = file.replace('_{}'.format(suffix), '')
                    uncal_files.append(os.path.join(input_file_directory, file.replace('.fits', '_uncal.fits')))
            elif isinstance(self.parameters['dead_flux_check_files'], list):
                uncal_files = self.parameters['dead_flux_check_files']
            else:
                raise ValueError('ERROR: dead_flux_check must be either a list, None, or "none"')
        else:
            uncal_files = self.parameters['dead_flux_check_files']

        bpm.find_bad_pix(self.inputimages,
                         dead_search=self.parameters['dead_search'],
                         low_qe_and_open_search=self.parameters['low_qe_and_open_search'],
                         dead_search_type=self.parameters['dead_search_type'],
                         sigma_threshold=self.parameters['sigma_threshold'],
                         normalization_method=self.parameters['normalization_method'],
                         smoothing_box_width=self.parameters['smoothing_box_width'],
                         smoothing_type=self.parameters['smoothing_type'],
                         dead_sigma_threshold=self.parameters['dead_sigma_threshold'],
                         max_dead_norm_signal=self.parameters['max_dead_norm_signal'],
                         run_dead_flux_check=self.parameters['run_dead_flux_check'],
                         dead_flux_check_files=uncal_files,
                         flux_check=self.parameters['flux_check'],
                         max_low_qe_norm_signal=self.parameters['max_low_qe_norm_signal'],
                         max_open_adj_norm_signal=self.parameters['max_open_adj_norm_signal'],
                         manual_flag_file=self.parameters['manual_flag_file'],
                         do_not_use=self.parameters['do_not_use'],
                         output_file=self.args.outputreffilename,
                         author=self.parameters['author'],
                         description=self.parameters['description'],
                         pedigree=self.parameters['pedigree'],
                         useafter=self.parameters['useafter'],
                         history=self.parameters['history'],
                         quality_check=self.parameters['quality_check'])
        return(0)


if __name__ == '__main__':
    """This should not need to be changed. This will read in the config
    files, import the script, generate the self.parameters dictionary, and
    run the argument parser above.
    """
    mkref = mkrefclass()
    mkref.make_reference_file()

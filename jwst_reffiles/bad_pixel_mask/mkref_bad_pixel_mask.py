#!/usr/bin/env python
"""
Plug-in script for the bad pixel mask creation module:
bad_pixel_mask.py, which uses the bad pxiel mask geneartion algorithms
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

from jwst_reffiles.bad_pixel_mask import bad_pixel_mask as bpm
from jwst_reffiles.utils.constants import RATE_FILE_SUFFIXES
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
        parser.add_argument('--flat_mean_sigma_threshold', help=('Number of standard deviations to use when sigma-'
                                                                'clipping to calculate the mean slope image or the mean '
                                                                'across the detector'))
        parser.add_argument('--flat_mean_normalization_method', help=('Specify how the mean image is normalized prior '
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
        parser.add_argument('--flat_do_not_use', help=('List of bad pixel types where the DO_NOT_USE flag should '
                                                  'also be applied (e.g. "["DEAD", "LOW_QE"]")'))

        parser.add_argument('--dark_stdev_clipping_sigma', help=('Number of sigma to use when sigma-clipping the 2D array of '
                                                                 'standard deviation values.'))

        parser.add_argument('--dark_max_clipping_iters', type=int, help=('Maximum number of iterations to use when sigma '
                                                                         'clipping to find the mean and standard deviation '
                                                                         'values used when locating noisy pixels.'))

        parser.add_argument('--dark_noisy_threshold', help=('Number of sigma above the mean noise (associated with the slope) '
                                                            'to use as a threshold for identifying noisy pixels.'))

        parser.add_argument('--max_saturated_fraction', help=('Fraction of integrations within which a pixel must be fully '
                                                              'saturated before flagging it as HOT.'))

        parser.add_argument('--max_jump_limit', type=int, help=('Maximum number of jumps a pixel can have in an integration '
                                                                'before it is flagged as a ``high jump`` pixel (which may be '
                                                                'flagged as noisy later).'))

        parser.add_argument('--jump_ratio_threshold', help=('Cutoff for the ratio of jumps early in the ramp to jumps later in '
                                                            'the ramp. Pixels with a ratio greater than this value (and which '
                                                            'also have a high total number of jumps) will be flagged as potential '
                                                            '(I)RC pixels.'))

        parser.add_argument('--early_cutoff_fraction', help=('Fraction of the integration to use when comparing the jump rate '
                                                             'early in the integration to that across the entire integration. '
                                                             'Must be <= 0.5'))

        parser.add_argument('--pedestal_sigma_threshold', help=('Used when searching for RC pixels via the pedestal image. Pixels '
                                                                'with pedestal values more than ``pedestal_sigma_threshold`` above '
                                                                'the mean are flagged as potential RC pixels.'))

        parser.add_argument('--rc_fraction_threshold', help=('Fraction of input files within which the pixel must be identified as '
                                                             'an RC pixel before it will be flagged as a permanent RC pixel.'))

        parser.add_argument('--low_pedestal_fraction', help=('Fraction of input files within which a pixel must be identified as '
                                                             'a low pedestal pixel before it will be flagged as a permanent low '
                                                             'pedestal pixel.'))

        parser.add_argument('--high_cr_fraction', help=('Fraction of input files within which a pixel must be flagged as having a '
                                                        'high number of jumps before it will be flagged as permanently noisy.'))

        parser.add_argument('--flag_values', help=('Dictionary mapping the types of bad pixels searched for to the flag mnemonics '
                                                   'to use when creating the bad pixel file. Keys are the types of bad pixels searched '
                                                   'for, and values are lists that include mnemonics recognized by the jwst calibration '
                                                   'pipeline e.g. {"hot": ["HOT"], "rc": ["RC"], "low_pedestal": ["OTHER_BAD_PIXEL"], "high_cr": ["TELEGRAPH"]}'))

        parser.add_argument('--dark_do_not_use', help=('List of bad pixel types to be flagged as DO_NOT_USE e.g. ["hot", "rc", "low_pedestal", "high_cr"]'))

        parser.add_argument('--plot', help=('If True, produce and save intermediate results from noisy pixel search'))

        parser.add_argument('--output_file', help=('Name of the CRDS-formatted bad pixel reference file to save the final bad pixel map into'))

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
        # Organize the input files into a group of darks and a group of
        # flats
        flatfiles = []
        darkfiles = []
        for row in self.inputimagestable:
            if row['imlabel'] == 'D':
                darkfiles.append(row['fitsfile'])
            elif row['imlabel'] == 'F':
                flatfiles.append(row['fitsfile'])

        # Since this module is called after calib_prep, all of the requested
        # outputs from the various pipeline steps should be present in the
        # output directory. Create lists of these files. The files listed in
        # self.inputimagestable['fitsfile'] should all be the same in terms
        # of their calibration state. So we should only have to check one in
        # order to know what state all of them are.
        dark_slope_files = []
        dark_uncal_files = []
        dark_jump_files = []
        dark_fitopt_files = []

        directory, filename = os.path.split(darkfiles[0])

        # Get the suffix of the input file so we know the calibration state
        suffix = None
        for ramp_suffix in RATE_FILE_SUFFIXES:
            if ramp_suffix in filename:
                dark_slope_files = copy.deepcopy(darkfiles)
                suffix = ramp_suffix
        if suffix is None:
            suffix = filename.split('_')[-1]
            suffix = suffix.replace('.fits', '')
            if suffix == 'uncal':
                dark_uncal_files = copy.deepcopy(darkfiles)
            elif suffix == 'jump':
                dark_jump_files = copy.deepcopy(darkfiles)
            else:
                raise ValueError('Unexpected suffixes for input dark files.')

        # Create lists of the needed calibration state files
        if len(dark_slope_files) > 0:
            dark_uncal_files = [elem.replace(suffix, '_uncal') for elem in dark_slope_files]
            dark_jump_files = [elem.replace(suffix, '_jump') for elem in dark_slope_files]
            dark_fitopt_files = [elem.replace(suffix, '_fitopt') for elem in dark_slope_files]
        elif len(dark_uncal_files) > 0:
            dark_slope_files = [elem.replace(suffix, '_1_ramp_fit') for elem in dark_uncal_files]
            dark_jump_files = [elem.replace(suffix, '_jump') for elem in dark_uncal_files]
            dark_fitopt_files = [elem.replace(suffix, '_fitopt') for elem in dark_uncal_files]
        elif len(dark_jump_files) > 0:
            dark_uncal_files = [elem.replace(suffix, '_uncal') for elem in dark_jump_files]
            dark_slope_files = [elem.replace(suffix, '_1_ramp_fit') for elem in dark_jump_files]
            dark_fitopt_files = [elem.replace(suffix, '_fitopt') for elem in dark_jump_files]

        # Repeat for flat field files
        flat_slope_files = []
        flat_uncal_files = []

        directory, filename = os.path.split(flatfiles[0])

        # Get the suffix of the input file so we know the calibration state
        suffix = None
        for ramp_suffix in RATE_FILE_SUFFIXES:
            if ramp_suffix in filename:
                flat_slope_files = copy.deepcopy(flatfiles)
                suffix = ramp_suffix
        if suffix is None:
            suffix = filename.split('_')[-1]
            suffix = suffix.replace('.fits', '')
            if suffix == 'uncal':
                flat_uncal_files = copy.deepcopy(flatfiles)
            else:
                raise ValueError('Unexpected suffixes for input flat field files.')

        # Create lists of the needed calibration state files
        if len(flat_slope_files) > 0:
            flat_uncal_files = [elem.replace(suffix, '_uncal') for elem in flat_slope_files]
        elif len(flat_uncal_files) > 0:
            flat_slope_files = [elem.replace(suffix, '_1_ramp_fit') for elem in flat_uncal_files]


        # The bad pixel mask module needs to use the file with the individual
        # slopes (_1_ramp_fit.fits), rather than the mean slope (_0_ramp_fit.fits),
        # for exposures with more than one integration per exposure. But for
        # exposures with only one integration, only the _0_ramp_fit file will be
        # produced. So go through the lists of slope files and check to see
        # which versions are present, and adjust the lists accordingly.




        # Call the wrapped module and provide the proper arguments from the
        # self.parameters dictionary.
        bpm.bad_pixels(flat_slope_files=flat_slope_files,
                       dead_search=self.parameters['dead_search'],
                       low_qe_and_open_search=self.parameters['low_qe_and_open_search'],
                       dead_search_type=self.parameters['dead_search_type'],
                       flat_mean_sigma_threshold=self.parameters['flat_mean_sigma_threshold'],
                       flat_mean_normalization_method=self.parameters['flat_mean_normalization_method'],
                       smoothing_box_width=self.parameters['smoothing_box_width'],
                       smoothing_type=self.parameters['smoothing_type'],
                       dead_sigma_threshold=self.parameters['dead_sigma_threshold'],
                       max_dead_norm_signal=self.parameters['max_dead_norm_signal'],
                       run_dead_flux_check=self.parameters['run_dead_flux_check'],
                       dead_flux_check_files=flat_uncal_files,
                       flux_check=self.parameters['flux_check'],
                       max_low_qe_norm_signal=self.parameters['max_low_qe_norm_signal'],
                       max_open_adj_norm_signal=self.parameters['max_open_adj_norm_signal'],
                       manual_flag_file=self.parameters['manual_flag_file'],
                       flat_do_not_use=self.parameters['flat_do_not_use'],
                       dark_slope_files=dark_slope_files,
                       dark_uncal_files=dark_uncal_files,
                       dark_jump_files=dark_jump_files,
                       dark_fitopt_files=dark_fitopt_files,
                       dark_stdev_clipping_sigma=self.parameters['dark_stdev_clipping_sigma'],
                       dark_max_clipping_iters=self.parameters['dark_max_clipping_iters'],
                       dark_noisy_threshold=self.parameters['dark_noisy_threshold'],
                       max_saturated_fraction=self.parameters['max_saturated_fraction'],
                       max_jump_limit=self.parameters['max_jump_limit'],
                       jump_ratio_threshold=self.parameters['jump_ratio_threshold'],
                       early_cutoff_fraction=self.parameters['early_cutoff_fraction'],
                       pedestal_sigma_threshold=self.parameters['pedestal_sigma_threshold'],
                       rc_fraction_threshold=self.parameters['rc_fraction_threshold'],
                       low_pedestal_fraction=self.parameters['low_pedestal_fraction'],
                       high_cr_fraction=self.parameters['high_cr_fraction'],
                       flag_values=self.parameters['flag_values'],
                       dark_do_not_use=self.parameters['dark_do_not_use'],
                       plot=self.parameters['plot'],
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

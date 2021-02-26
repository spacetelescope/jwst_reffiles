#! /usr/bin/env python

"""Wrapper script that calls badpix_from_flats.py and badpix_from_darks.py.
The results are combined to create a single bad pixel reference file.
"""
import copy
import datetime
import os

from astropy.io import fits
from jwst.datamodels import MaskModel, util
import numpy as np

from jwst_reffiles.bad_pixel_mask import badpix_from_flats
from jwst_reffiles.dark_current import badpix_from_darks

# Flat field-related header keywords
dead_search_kw = 'BPFDEAD'
low_qe_search_kw = 'BPFLOWQE'
dead_search_type_kw = 'BPFSCHTP'
mean_sig_threshold_kw = 'BPFSIGMA'
norm_method_kw = 'BPFNORM'
smooth_box_width_kw = 'BPFSMOTH'
smoothing_type_kw = 'BPFSMTYP'
dead_sig_thresh_kw = 'BPFDEDSG'
dead_zero_sig_frac_kw = 'BPFZEROF'
dead_flux_check_kw = 'BPFFXCHK'
#dead_flux_file_kw = 'BPFFXFIL'
max_dead_sig_kw = 'BPFMXDED'
manual_flag_kw = 'BPFMANFL'
flat_do_not_use_kw = 'BPFDONOT'
max_low_qe_kw = 'BPFMXLQE'
max_open_adj_kw = 'BPFMXOAD'

# Dark current-related header keywords
bad_from_dark_kw = 'BPDSERCH'
dark_clip_sigma_kw = 'BPDCLPSG'
dark_clip_iters_kw = 'BPDCLPIT'
dark_noisy_thresh_kw = 'BPDNSETH'
max_sat_frac_kw = 'BPDMXSAT'
jump_limit_kw = 'BPDMXJMP'
jump_ratio_thresh_kw = 'BPDJPRAT'
cutoff_frac_kw = 'BPDCUTFC'
pedestal_sig_thresh_kw = 'BPDPEDTH'
rc_frac_thresh_kw = 'BPDRCTH'
low_ped_frac_kw = 'BPDLOPFC'
high_cr_frac_kw = 'BPDCRFC'
dark_do_not_use_kw = 'BPDDONOT'
flag_mapping_kw = 'BPDMAPS'


def bad_pixels(flat_slope_files=None, dead_search=True, low_qe_and_open_search=True,
               dead_search_type='sigma_rate', flat_mean_sigma_threshold=3, flat_mean_normalization_method='smoothed',
               smoothing_box_width=15, smoothing_type='Box2D', dead_sigma_threshold=5.,  max_dead_norm_signal=None,
               run_dead_flux_check=False, dead_flux_check_files=None, flux_check=45000, max_low_qe_norm_signal=0.5,
               max_open_adj_norm_signal=1.05, manual_flag_file='default', flat_do_not_use=[],
               dark_slope_files=None, dark_uncal_files=None, dark_jump_files=None, dark_fitopt_files=None,
               dark_stdev_clipping_sigma=5., dark_max_clipping_iters=5,
               dark_noisy_threshold=5, max_saturated_fraction=0.5, max_jump_limit=10, jump_ratio_threshold=5,
               early_cutoff_fraction=0.25, pedestal_sigma_threshold=5, rc_fraction_threshold=0.8, low_pedestal_fraction=0.8,
               high_cr_fraction=0.8,
               flag_values={'hot': ['HOT'], 'rc': ['RC'], 'low_pedestal': ['OTHER_BAD_PIXEL'], 'high_cr': ["TELEGRAPH"]},
               dark_do_not_use=['hot', 'rc', 'low_pedestal', 'high_cr'], plot=False,
               output_file=None, author='jwst_reffiles', description='A bad pix mask',
               pedigree='GROUND', useafter='2019-04-01 00:00:00', history='', quality_check=False):
    """
    Wrapper that calls the two modules for finding bad pixels from input flat
    field files, and bad pixels from dark current files.

    Parameters
    ----------
    flat_slope_files : list
        List of flat field slope files to be used for the dead pixel search.
        If None, the search is skipped.

    dead_search : bool
        Whether or not to search for DEAD pixels using the flat field files

    low_qe_and_open_search : bool
        Whether or not to search for LOW_QE, OPEN, and ADJ_OPEN pixels using
        the flat field files

    dead_search_type : str
        Type of search to use when looking for dead pixels. Options are:
        ``sigma_rate``: Using a normalized signal rate image, dead pixels
                        are defined as those with a rate smaller than
                        ``dead_sigma_threshold`` standard deviations below
                        the mean.
        ``absolute_rate``: Using a normalized signal rate image, dead pixels
                           are defined as those with a rate less than
                           ``max_dead_norm_signal``.

    flat_mean_sigma_threshold : float
        Number of standard deviations to use when sigma-clipping to
        calculate the mean slope image or the mean across the detector

    flat_mean_normalization_method : str
        Specify how the mean image is normalized prior to searching for
        bad pixels. Options are:
        'smoothed': Mean image will be smoothed using a
                    ``smoothing_box_width`` x ``smoothing_box_width``
                    box kernel. The mean image is then normalized by
                    this smoothed image.
        'none': No normalization is done. Mean slope image is used as is
        'mean': Mean image is normalized by its sigma-clipped mean
        'fit2d': Mean image will be normalized by a fit of 2-D surface to
                 mean image. The degree of the fit is controlled by the
                 ``fit_degree`` parameters

    smoothing_box_width : float
        Width in pixels of the box kernel to use to compute the smoothed
        mean image

    smoothing_typ : string
        Type of smoothing to do ``Box2D `` or ``median`` filtering

    smoothing_sigma : float
        Number of standard deviations to use when smoothing in a box defined
        by smoothing_box_width.

    dead_sigma_threshold : float
        Number of standard deviations below the mean at which a pixel is
        considered dead.

    max_dead_norm_signal : float
        Maximum normalized signal rate of a pixel that is considered dead

    run_dead_flux_check : bool
        Whether or not to check for dead pixels using an absolute flux value

    dead_flux_check_files : list
        List of ramp (uncalibrated) files to use to check the flux of average
        of last 4 groups. If None then the uncalibration files are not read in
        and no flux_check is done.

    flux_check: float
        Tolerance on average signal in last 4 groups. If dead_flux_check_files is
        a list of uncalibrated files, then the average of the last four groups
        for all the integrations is determined. If this average > flux_check
        then this pixel is not a dead pixel.

    max_low_qe_norm_signal: float
        The maximum normalized signal a pixel can have and be considered
        low QE.

    max_open_adj_norm_signal : float
        The maximum normalized signal a pixel adjacent to a low QE pixel can have
        in order for the low QE pixel to be reclassified as OPEN

    manual_flag_file : str
        Name of the ascii file containing a list of pixels to be added manually
        to the output bad pixel mask file. Default is 'default', in which case
        the file contained in the ``bad_pixel_mask`` directory of the repo will
        be used.

    flat_do_not_use : list
        List of bad pixel types where the DO_NOT_USE flag should also be
        applied (e.g. ['DEAD', 'LOW_QE'])

    dark_slope_files : list
        List of dark current slope files to be used for the noisy pixel search.
        If None, the search is skipped.

    dark_uncal_files : list
        List of uncalibrated dark current ramp files. These should correspond
        1-to-1 with the files listed in ``dark_slope_files``. If None,
        the code assumes the files are in the same location as the slope
        files and have names ending in uncal.fits

    dark_jump_files : list
        List of dark current ramp files output from the jump step of the pipeline.
        These should correspond 1-to-1 with the files listed in ``dark_slope_files``.
        If None, the code assumes the files are in the same location as the slope
        files and have names ending in jump.fits

    dark_fitopt_files : list
        List of optional output files produced by the ramp-fitting step of the
        pipeline. These should correspond 1-to-1 with the files listed in
        ``dark_slope_files``. If None, the code assumes the files are in the
        same location as the slope files and have names ending in fitopt.fits

    dark_stdev_clipping_sigma : int
        Number of sigma to use when sigma-clipping the 2D array of
        standard deviation values. The sigma-clipped mean and standard
        deviation are used to locate noisy pixels.

    dark_max_clipping_iters : int
        Maximum number of iterations to use when sigma clipping to find
        the mean and standard deviation values that are used when
        locating noisy pixels.

    dark_noisy_threshold : int
        Number of sigma above the mean noise (associated with the slope)
        to use as a threshold for identifying noisy pixels.

    max_saturated_fraction : float
        When identifying pixels that are fully saturated (in all groups
        of an integration), this is the fraction of integrations within
        which a pixel must be fully saturated before flagging it as HOT

    max_jump_limit : int
        The maximum number of jumps a pixel can have in an integration
        before it is flagged as a ``high jump`` pixel (which may be
        flagged as noisy later)

    jump_ratio_threshold : int
        Cutoff for the ratio of jumps early in the ramp to jumps later in
        the ramp. Pixels with a ratio greater than this value (and which
        also have a high total number of jumps) will be flagged as
        potential (I)RC pixels.

    early_cutoff_fraction : float
        Fraction of the integration to use when comparing the jump rate
        early in the integration to that across the entire integration.
        Must be <= 0.5

    pedestal_sigma_threshold : int
        Used when searching for RC pixels via the pedestal image. Pixels
        with pedestal values more than ``pedestal_sigma_threshold`` above
        the mean are flagged as potential RC pixels

    rc_fraction_threshold : float
        Used when searching for RC pixels. This is the fraction of input
        files within which the pixel must be identified as an RC pixel
        before it will be flagged as a permanent RC pixel

    low_pedestal_fraction : float
        This is the fraction of input files within which a pixel must be
        identified as a low pedestal pixel before it will be flagged as
        a permanent low pedestal pixel

    high_cr_fraction : float
        This is the fraction of input files within which a pixel must be
        flagged as having a high number of jumps before it will be flagged
        as permanently noisy

    flag_values : dict
        This dictionary maps the types of bad pixels searched for to the
        flag mnemonics to use when creating the bad pixel file. Keys are
        the types of bad pixels searched for, and values are lists that
        include mnemonics recognized by the jwst calibration pipeline
        e.g. {'hot': ['HOT'], 'rc': ['RC'], 'low_pedestal': ['OTHER_BAD_PIXEL'], 'high_cr': ["TELEGRAPH"]}

    dark_do_not_use : list
        List of bad pixel types to be flagged as DO_NOT_USE
        e.g. ['hot', 'rc', 'low_pedestal', 'high_cr']

    plot : bool
        If True, produce and save intermediate results from noisy pixel search

    output_file : str
        Name of the CRDS-formatted bad pixel reference file to save the final
        bad pixel map into

    author : str
        CRDS-required name of the reference file author, to be placed in the
        referece file header

    description : str
        CRDS-required description of the reference file, to be placed in the
        reference file header

    pedigree : str
        CRDS-required pedigree of the data used to create the reference file

    useafter : str
        CRDS-required date of earliest data with which this referece file
        should be used. (e.g. '2019-04-01 00:00:00')

    history : str
        Text to be added to the HISOTRY section of the output bad pixel file

    quality_check : bool
        If True, the pipeline is run using the output reference file to be
        sure the pipeline doens't crash
    """
    instrument = None
    detector = None
    all_files = []
    history = [history]
    hdu = fits.PrimaryHDU()

    if flat_slope_files is not None:
        all_files = copy.deepcopy(flat_slope_files)
        instrument, detector = instrument_info(flat_slope_files[0])

        # Get output filenames
        if output_file is None:
            output_file = create_output_filename(instrument, detector)

        flat_output_file = output_file.replace('.fits', '_from_flats.fits')

        # Get bad pixels from the flats
        flatmask = badpix_from_flats.find_bad_pix(flat_slope_files, dead_search=dead_search,
                                                  low_qe_and_open_search=low_qe_and_open_search,
                                                  dead_search_type=dead_search_type,
                                                  sigma_threshold=flat_mean_sigma_threshold,
                                                  normalization_method=flat_mean_normalization_method,
                                                  smoothing_type=smoothing_type,
                                                  smoothing_box_width=smoothing_box_width,
                                                  dead_sigma_threshold=dead_sigma_threshold,
                                                  #dead_zero_signal_fraction=dead_zero_signal_fraction,
                                                  run_dead_flux_check=run_dead_flux_check,
                                                  dead_flux_check_files=dead_flux_check_files,
                                                  max_dead_norm_signal=max_dead_norm_signal,
                                                  manual_flag_file=manual_flag_file,
                                                  max_low_qe_norm_signal=max_low_qe_norm_signal,
                                                  max_open_adj_norm_signal=max_open_adj_norm_signal,
                                                  do_not_use=flat_do_not_use,
                                                  output_file=flat_output_file,
                                                  author=author,
                                                  description=description,
                                                  pedigree=pedigree,
                                                  useafter=useafter,
                                                  history=history[0],
                                                  quality_check=quality_check)

        # Convert the do not use list to a string to add to the header
        if len(flat_do_not_use) > 0:
            flat_do_not_use_string = ', '.join(flat_do_not_use)
        else:
            flat_do_not_use_string = 'None'
        flat_do_not_use_string = '{} {}'.format('Bad pixel types from flat to which DO_NOT_USE is applied: ', flat_do_not_use_string)

        # Add the do not use string to the list of history entries to add,
        # since it may end up being longer than 8 characters
        history.append(flat_do_not_use_string)

        # Define the non-standard fits header keywords by placing them in a
        # fits HDU List
        hdu.header[dead_search_kw] = dead_search
        hdu.header[low_qe_search_kw] = low_qe_and_open_search
        hdu.header[dead_search_type_kw] = dead_search_type
        hdu.header[mean_sig_threshold_kw] = flat_mean_sigma_threshold
        hdu.header[norm_method_kw] = flat_mean_normalization_method
        hdu.header[smooth_box_width_kw] = smoothing_box_width
        hdu.header[dead_sig_thresh_kw] = dead_sigma_threshold
        #hdu.header[dead_zero_sig_frac_kw] = dead_zero_signal_fraction
        hdu.header[dead_flux_check_kw] = run_dead_flux_check
        #hdu.header[dead_flux_file_kw] = dead_flux_check_files
        hdu.header[max_dead_sig_kw] = max_dead_norm_signal
        hdu.header[manual_flag_kw] = manual_flag_file
        hdu.header[max_low_qe_kw] = max_low_qe_norm_signal
        hdu.header[max_open_adj_kw] = max_open_adj_norm_signal

    else:
        flatmask = 0
        hdu.header[dead_search_kw] = False
        hdu.header[low_qe_search_kw] = False

    if dark_slope_files is not None:
        if len(all_files) == 0:
            all_files = copy.deepcopy(dark_slope_files)
            instrument, detector = instrument_info(dark_slope_files[0])
        else:
            all_files = all_files + dark_slope_files

        # Get output filenames
        if output_file is None:
            output_file = create_output_filename(instrument, detector)

        dark_output_file = output_file.replace('.fits', '_from_darks.fits')

        # Get bad pixels from the darks
        darkmask = badpix_from_darks.find_bad_pix(dark_slope_files, uncal_filenames=dark_uncal_files,
                                                  jump_filenames=dark_jump_files,
                                                  fitopt_filenames=dark_fitopt_files,
                                                  clipping_sigma=dark_stdev_clipping_sigma,
                                                  max_clipping_iters=dark_max_clipping_iters,
                                                  noisy_threshold=dark_noisy_threshold,
                                                  max_saturated_fraction=max_saturated_fraction,
                                                  max_jump_limit=max_jump_limit,
                                                  jump_ratio_threshold=jump_ratio_threshold,
                                                  early_cutoff_fraction=early_cutoff_fraction,
                                                  pedestal_sigma_threshold=pedestal_sigma_threshold,
                                                  rc_fraction_threshold=rc_fraction_threshold,
                                                  low_pedestal_fraction=low_pedestal_fraction,
                                                  high_cr_fraction=high_cr_fraction,
                                                  flag_values=flag_values,
                                                  do_not_use=dark_do_not_use,
                                                  outfile=dark_output_file, plot=False)

        # Convert the do not use list to a string to add to the header
        if len(dark_do_not_use) > 0:
            dark_do_not_use_string = ', '.join(dark_do_not_use)
        else:
            dark_do_not_use_string = 'None'
        dark_do_not_use_string = '{} {}'.format('Bad pixel types from dark to which DO_NOT_USE is applied: ', dark_do_not_use_string)

        # Add the do not use string to the list of history entries to add,
        # since it may end up being longer than 8 characters
        history.append(dark_do_not_use_string)

        # Convert the bad pixel type mapping into a string so it can be
        # added to the output header
        if len(flag_values) > 0:
            mapping_str = ''
            for key in flag_values:
                substr = '{}: {}, '.format(key, flag_values[key])
                mapping_str = mapping_str + substr
        else:
            mapping_str = 'None'
        mapping_str = '{} {}'.format('Mapping of jwst_reffiles bad pixel types to jwst cal bad pixel flags: ', mapping_str)

        # Add the do not use string to the list of history entries to add,
        # since it may end up being longer than 8 characters
        history.append(mapping_str)

        # Define the non-standard fits header keywords by placing them in a
        # fits HDU List
        hdu.header[bad_from_dark_kw] = True
        hdu.header[dark_clip_sigma_kw] = dark_stdev_clipping_sigma
        hdu.header[dark_clip_iters_kw] = dark_max_clipping_iters
        hdu.header[dark_noisy_thresh_kw] = dark_noisy_threshold
        hdu.header[max_sat_frac_kw] = max_saturated_fraction
        hdu.header[jump_limit_kw] = max_jump_limit
        hdu.header[jump_ratio_thresh_kw] = jump_ratio_threshold
        hdu.header[cutoff_frac_kw] = early_cutoff_fraction
        hdu.header[pedestal_sig_thresh_kw ] = pedestal_sigma_threshold
        hdu.header[rc_frac_thresh_kw ] = rc_fraction_threshold
        hdu.header[low_ped_frac_kw] = low_pedestal_fraction
        hdu.header[high_cr_frac_kw] = high_cr_fraction

    else:
        darkmask = 0.
        hdu.header[bad_from_dark_kw] = False

    # Combine the two masks
    final_mask = np.bitwise_or(flatmask, darkmask)

    # Some pixels that are saturated in all groups may be flagged as hot and dead,
    # because the slope in the flat ramp appears to be zero, while the more in-depth
    # checking of the darks shows that the pixel is saturated the entire time. For any
    # pixels flagged as both, keep only the hot flag and throw out the dead flag.
    hot_and_or_dead = (final_mask & dqflags.pixel['HOT']) + (final_mask & dqflags.pixel['DEAD'])
    hot_and_dead = hot_and_or_dead == (dqflags.pixel['HOT'] + dqflags.pixel['DEAD'])
    final_mask[hot_and_dead] = np.bitwise_xor(final_mask[hot_and_dead], dqflags.pixel['DEAD'])

    # Save mask in reference file
    hdu_list = fits.HDUList([hdu])
    save_final_map(final_mask, instrument.upper(), detector.upper(), hdu_list, all_files, author, description,
                   pedigree, useafter, history, output_file)


def create_output_filename(inst_name, det_name):
    """Create a default output filename for the bad pixel mask given
    instrument and detector names

    Parameters
    ----------
    inst_name : str
        Instrument name

    det_name : str
        Detector name

    Returns
    -------
    outfile : str
        Default bad pixel mask filename
    """
    # Add in timestamp as a way to prevent overwriting past runs
    current_time = datetime.datetime.now()

    # Use the current working directory
    outfile = '{}_{}_{}_badpix_mask.fits'.format(inst_name, det_name, current_time)
    outfile = os.path.join(os.getcwd(), outfile)
    return outfile


def instrument_info(filename):
    """Get the instrument and detector name from the header of the
    input file

    Parameters
    ----------
    filename : str
        Name of fits file

    Returns
    -------
    inst : str
        Instrument name

    det : str
        Detector name
    """
    with fits.open(filename) as hdulist:
        try:
            inst = hdulist[0].header['INSTRUME'].lower()
        except KeyError:
            raise KeyError("ERROR: expecting instrument name in main header of {}".format(filename))

        try:
            det = hdulist[0].header['DETECTOR'].lower()
        except KeyError:
            raise KeyError("ERROR: expecting detector name in main header of {}".format(filename))
    return inst, det


def save_final_map(bad_pix_map, instrument, detector, hdulist, files, author, description, pedigree, useafter,
                   history_text, outfile):
    """Save a bad pixel map into a CRDS-formatted reference file

    Parameters
    ----------
    bad_pix_map : numpy.ndarray
        2D bad pixel array

    instrument : str
        Name of instrument associated with the bad pixel array

    detector : str
        Name of detector associated with the bad pixel array

    hdulist : astropy.fits.HDUList
        HDUList containing "extra" fits keywords

    files : list
        List of files used to create ``bad_pix_map``

    author : str
        Author of the bad pixel mask reference file

    description : str
        CRDS description to use in the final bad pixel file

    pedigree : str
        CRDS pedigree to use in the final bad pixel file

    useafter : str
        CRDS useafter string for the bad pixel file

    history_text : list
        List of strings to add as HISTORY entries to the bad pixel file

    outfile : str
        Name of the output bad pixel file
    """
    yd, xd = bad_pix_map.shape

    # Initialize the MaskModel using the hdu_list, so the new keywords will
    # be populated
    model = MaskModel(hdulist)
    model.dq = bad_pix_map

    # Create dq_def data
    dq_def = badpix_from_flats.create_dqdef()
    model.dq_def = dq_def
    model.meta.reftype = 'MASK'
    model.meta.subarray.name = 'FULL'
    model.meta.subarray.xstart = 1
    model.meta.subarray.xsize = xd
    model.meta.subarray.ystart = 1
    model.meta.subarray.ysize = yd
    model.meta.instrument.name = instrument.upper()
    model.meta.instrument.detector = detector

    # Get the fast and slow axis directions from one of the input files
    fastaxis, slowaxis = badpix_from_flats.get_fastaxis(files[0])
    model.meta.subarray.fastaxis = fastaxis
    model.meta.subarray.slowaxis = slowaxis

    model.meta.author = author
    model.meta.description = description
    model.meta.pedigree = pedigree
    model.meta.useafter = useafter

    # Add information about parameters used
    # Parameters from badpix_from_flats
    package_note = ('This file was created using the bad_pixel_mask.py module within the '
                    'jwst_reffiles package.')

    software_dict = {'name': 'jwst_reffiles.bad_pixel_mask.bad_pixel_mask.py', 'author': 'STScI',
                     'homepage': 'https://github.com/spacetelescope/jwst_reffiles',
                     'version': '0.0.0'}
    entry = util.create_history_entry(package_note, software=software_dict)
    model.history.append(entry)

    model.history.append(util.create_history_entry('Parameter values and descriptions:'))
    dead_search_descrip = ('dead_search: Boolean, whether or not to run the dead pixel search '
                           'using flat field files. The value is stored in the {} keyword.'.format(dead_search_kw))
    model.history.append(util.create_history_entry(dead_search_descrip))

    low_qe_search_descrip = ('low_qe_and_open_search: Boolean, whether or not to run the low QE '
                             'and open pixel search using flat field files. The value is stored in the {} '
                             'keyword.'.format(low_qe_search_kw))
    model.history.append(util.create_history_entry(low_qe_search_descrip))

    dead_type_descrip = ('dead_search_type: Method used to identify dead pixels. The value is stored in the '
                         '{} keyword.'.format(dead_search_type_kw))
    model.history.append(util.create_history_entry(dead_type_descrip))

    sigma_descrip = ('flat_mean_sigma_threshold: Number of standard deviations to use when sigma-clipping to '
                     'calculate the mean slope image or the mean across the detector. The value '
                     'used is stored in the {} keyword.'.format(mean_sig_threshold_kw))
    model.history.append(util.create_history_entry(sigma_descrip))

    norm_descrip = ('flat_mean_normalization_method: Specify how the mean image is normalized prior to searching '
        'for bad pixels. The value used is stored in the {} keyword.'.format(norm_method_kw))
    model.history.append(util.create_history_entry(norm_descrip))

    smooth_descrip = ('smoothing_box_width: Width in pixels of the box kernel to use to compute the '
                      'smoothed mean image. The value used is stored in the {} keyword.'.format(smooth_box_width_kw))
    model.history.append(util.create_history_entry(smooth_descrip))

    smooth_type_descrip = ('smoothing_type: Type of smoothing to do: Box2D or median filtering. The value used '
                           'is stored in the {} keyword.'.format(smoothing_type_kw))
    model.history.append(util.create_history_entry(smooth_type_descrip))

    dead_sig_descrip = ('Number of standard deviations below the mean at which a pixel is considered dead. '
                        'The value used is stored in the {} keyword.'.format(dead_sig_thresh_kw))
    model.history.append(util.create_history_entry(dead_sig_descrip))

    max_dead_descrip = ('Maximum normalized signal rate of a pixel that is considered dead. The value '
                        'used is stored in the {} keyword.'.format(max_dead_sig_kw))
    model.history.append(util.create_history_entry(max_dead_descrip))

    run_dead_flux_descrip = ('run_dead_flux_check: Boolean, if True, search for pixels erroneously flagged '
                             'as dead because they are saturated in all groups. The value used is stored '
                             'in the {} keyword.'.format(dead_flux_check_kw))
    model.history.append(util.create_history_entry(run_dead_flux_descrip))

    dead_flux_limit_descrip = ('Signal limit in raw data above which the pixel is considered not dead. The '
                               'value used is stored in the {} keyword.'.format(max_dead_sig_kw))
    model.history.append(util.create_history_entry(dead_flux_limit_descrip))

    max_low_qe_descrip = ('The maximum normalized signal a pixel can have and be considered low QE. The '
                          'value used is stored in the {} keyword.'.format(max_low_qe_kw))
    model.history.append(util.create_history_entry(max_low_qe_descrip))

    max_open_adj_descrip = ('The maximum normalized signal a pixel adjacent to a low QE pixel can have '
                            'in order for the low QE pixel to be reclassified as OPEN. The value used '
                            'is stored in the {} keyword.'.format(max_open_adj_kw))
    model.history.append(util.create_history_entry(max_open_adj_descrip))

    flat_do_not_use_descrip = ('List of bad pixel types (from flats) where the DO_NOT_USE flag is also applied. '
                          'The values used are stored in the {} keyword.'.format(flat_do_not_use_kw))
    model.history.append(util.create_history_entry(flat_do_not_use_descrip))

    manual_file_descrip = ('Name of the ascii file containing a list of pixels to be added manually. The '
                           'value used is stored in the {} keyword.'.format(manual_flag_kw))
    model.history.append(util.create_history_entry(manual_file_descrip))

    # Parameters from badpix_from_darks
    bad_from_dark_descrip = ('badpix_from_dark: Boolean, whether or not the bad pixel from dark search  '
                             'has been run. The value is stored in the {} keyword.'.format(bad_from_dark_kw))
    model.history.append(util.create_history_entry(bad_from_dark_descrip))

    dark_clip_sig_descrip = ('Number of sigma to use when sigma-clipping 2D stdev image. The value used '
                             'is stored in the {} keyword.'.format(dark_clip_sigma_kw))
    model.history.append(util.create_history_entry(dark_clip_sig_descrip))

    dark_clip_iter_descrip = ('Max number of iterations to use when sigma clipping mean and stdev values. '
                             'The value used is stored in the {} keyword.'.format(dark_clip_iters_kw))
    model.history.append(util.create_history_entry(dark_clip_iter_descrip))

    dark_noisy_thresh_descrip = ('Number of sigma above mean noise for noisy pix threshold. The value '
                                 'used is stored in the {} keyword.'.format(dark_noisy_thresh_kw))
    model.history.append(util.create_history_entry(dark_noisy_thresh_descrip))

    max_sat_frac_descrip = ('Fraction of integrations within which a pixel must be fully saturated before '
                            'flagging it as HOT. The value used is stored in the {} keyword.'.format(max_sat_frac_kw))
    model.history.append(util.create_history_entry(max_sat_frac_descrip))

    jump_limit_descrip = ('Maximum number of jumps a pixel can have in an integration before it is flagged as a '
                          '"high jump" pixel. The value used is stored in the {} keyword.'.format(jump_limit_kw))
    model.history.append(util.create_history_entry(jump_limit_descrip))

    jump_ratio_descrip = ('Cutoff for the ratio of jumps early in the ramp to jumps later in the ramp when '
                          'looking for RC pixels. The value used is stored in the {} keyword.'.format(jump_ratio_thresh_kw))
    model.history.append(util.create_history_entry(jump_ratio_descrip))

    cutoff_frac_descrip = ('Fraction of the integration to use when comparing the jump rate early in the integration to '
                           'that across the entire integration. The value used is stored in the {} keyword.'.format(cutoff_frac_kw))
    model.history.append(util.create_history_entry(cutoff_frac_descrip))

    ped_sigma_descrip = ('Pixels with pedestal values more than this limit above the mean are flagged as RC. '
                         'The value used is stored in the {} keyword.'.format(pedestal_sig_thresh_kw))
    model.history.append(util.create_history_entry(ped_sigma_descrip))

    rc_thresh_descrip = ('Fraction of input files within which a pixel must be identified as an RC pixel before '
                         'it will be flagged as a permanent RC pixel. The value used is stored in the {} '
                         'keyword.'.format(rc_frac_thresh_kw))
    model.history.append(util.create_history_entry(rc_thresh_descrip))

    low_ped_descrip = ('Fraction of input files within which a pixel must be identified as a low pedestal '
                       'pixel before it will be flagged as a permanent low pedestal pixel. The value used '
                       'is stored in the {} keyword.'.format(low_ped_frac_kw))
    model.history.append(util.create_history_entry(low_ped_descrip))

    high_cr_descrip = ('Fraction of input files within which a pixel must be flagged as having a high number '
                       'of jumps before it will be flagged as permanently noisy. The value used '
                       'is stored in the {} keyword.'.format(high_cr_frac_kw))
    dark_do_not_use_descrip = ('List of bad pixel types (from darks) where the DO_NOT_USE flag is also applied. '
                          'The values used are stored in the {} keyword.'.format(dark_do_not_use_kw))
    model.history.append(util.create_history_entry(dark_do_not_use_descrip))

    # Add the list of input files used to create the map
    model.history.append('DATA USED:')
    for file in files:
        totlen = len(file)
        div = np.arange(0, totlen, 60)
        for val in div:
            if totlen > (val+60):
                model.history.append(util.create_history_entry(file[val:val+60]))
            else:
                model.history.append(util.create_history_entry(file[val:]))

    # Add the do not use lists, pixel flag mappings, and user-provided
    # history text
    for history_entry in history_text:
        if history_entry != '':
            model.history.append(util.create_history_entry(history_entry))

    model.save(outfile, overwrite=True)
    print('Final bad pixel mask reference file save to: {}'.format(outfile))

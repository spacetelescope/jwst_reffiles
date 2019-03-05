#! /usr/bin/env python

"""This module can be used to create JWST-format bad pixel mask reference
files for use in the ``dq_init`` step of the JWST calibration pipeline.

Author
------
     - Bryan Hilbert
 Use
---
     This module can be imported and used as such:
     ::
         from jwst_reffiles.bad_pixel_mask import bad_pixel_mask
         bad_pixel_mask.something(arguments)

         or

         command line call here

Notes
-----
    This algorithm is used to identify types of bad pixels that are
    flagged in the bad pixel mask. This includes the following types:
    DEAD, LOW_QE, OPEN, ADJ_OPEN

    This is based on discussions within the JWST Reference File Generation
    Working Group in February 2019. The goal is to produce an algorithm
    that can be used by all JWST instruments to produce bad pixel mask
    reference files.

    Overview:
    Inputs: A set of flatfield countrate images (NIRCam, NIRISS, NIRSpec)
            A set of dark current exposures i.e. Ramps (MIRI)

    Algorithm:
        0. For MIRI, extract and use only the 10th group from each exposure
        1. Create average image from the set of input images, using sigma-clipping
        2. Create smoothed version of the average image (15x15 smoothing)
        3. Divide average image by smoothed image
        4. Calculate sigma-clipped mean and stdev in normalized image
        5. Find bad pix:
            NIRCam, NIRISS: DEAD+DO_NOT_USE if signal < (mean-N*stdev)
            MIRI: DEAD+DO_NOT_USE if slope = 0 in >90% of the input images
            NIRSPec: DEAD+DO_NOT_USE if signal < 0.05
                     LOW_QE if 0.05 < signal < 0.5 and signal in 4 adjacent pixels
                     are below ADJ_OPEN threshold (1.05)

References
----------
    Working group discussion and algorithm details are presented at:
    https://outerspace.stsci.edu/display/JWSTCC/Algorithm+details%3A+DQ+Init

Dependencies
------------
    - jwst
    - astropy
    - numpy
"""
import argparse
import numpy as np
from scipy.stats import sigmaclip

from astropy.convolution import convolve, Box2DKernel
from astropy.io import fits
from jwst.datamodels import dqflags, util, MaskModel, Level1bModel
from jwst.dq_init import DQInitStep


def find_bad_pix(input_files, sigma_threshold=3, smoothing_box_width=15, dead_sigma_threshold=5.,
                 max_dead_norm_signal=0.05, max_low_qe_norm_signal=0.5, max_open_adj_norm_signal=1.05,
                 do_not_use=[], output_file=None, author='jwst_reffiles', description='A bad pix mask',
                 pedigree='GROUND', useafter='2019-04-01 00:00:00', quality_check=True):
    """MAIN SCRIPT: Given a set of input files, create dead, low QE, and
    open pixel maps

    Parameters
    ----------
    input_files : list
        Fits files to be used in the search for dead/open pixels

    sigma_threshold : float
        Number of standard deviations to use when sigma-clipping to
        calculate the mean slope image or the mean across the detector

    smoothing_box_width : float
        Width in pixels of the box kernel to use to compute the smoothed
        mean image

    dead_sigma_threshold : float
        Number of standard deviations below the mean at which a pixel is
        considered dead.

    max_dead_norm_signal : float
        Maximum normalized signal rate of a pixel that is considered dead

    max_low_qe_norm_signal: float
        The maximum normalized signal a pixel can have and be considered
        low QE.

    max_open_adj_norm_signal : float
        The maximum normalized signal a pixel adjacent to a low QE pixel can have
        in order for the low QE pixel to be reclassified as OPEN

    do_not_use : list
        List of bad pixel types where the DO_NOT_USE flag should also be
        applied (e.g. ['DEAD', 'LOW_QE'])

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

    quality_check : bool
        If True, the pipeline is run using the output reference file to be
        sure the pipeline doens't crash

    Returns
    -------
    dead_map : numpy.ndarray
        Map of dead pixel locations. 1 = dead, 0 = not dead

    lowqe_map : numpy.ndarray
        Map of low QE pixel locations. 1 = low QE, 0 = not low QE

    open_map : numpy.ndarray
        Map of open and adjacent to open pixel locations. 1 = open,
        2 = adjacent to open, 0 = good.
    """
    crds_input_checks(author, description, pedigree, useafter)

    # Read in input files
    input_exposures, instrument, detector = read_files(input_files)

    if instrument == 'MIRI':
        dead_map = find_miri_dead_pixels(input_exposures)

    # Create mean and stdev images
    mean_img, stdev_img = mean_stdev_images(input_exposures, sigma=sigma_threshold)

    # Create smoothed version of mean image
    smoothed_image = smooth(mean_img, box_width=smoothing_box_width)

    # Normalize
    normalized = mean_img / smoothed_image

    # Sigma-clipped mean and stdev
    norm_mean, norm_dev = image_stats(normalized, sigma=sigma_threshold)

    if instrument != 'MIRI':
        # Find dead pixels
        dead_map = find_dead_pixels(normalized, norm_mean, norm_dev, sigma=dead_sigma_threshold,
                                    max_dead_signal=max_dead_norm_signal, max_low_qe=max_low_eq_norm_signal)

    # Find open pixels
    lowqe_map, open_map, adjacent_to_open_map = find_open_and_low_qe_pixels(normalized, max_dead_signal=max_dead_norm_signal,
                                                                            max_low_qe=max_low_qe_norm_signal,
                                                                            max_adj_open=max_open_adj_norm_signal)

    # Flag MIRI's bad columns
    bad_col_map = miri_bad_columns(dead_map.shape)

    # Create a map showing locations of reference pixels
    reference_pix = reference_pixel_map(dead_map.shape)

    # Wrap up all of the individual bad pixel maps into a dictionary
    stack_of_maps = {'DEAD': dead_map, 'LOW_QE': lowqe_map, 'OPEN': open_map,
                     'ADJ_OPEN': adjacent_to_open_map, 'REFERENCE_PIXEL': reference_pix}  # '???': bad_col_map,

    # Combine the individual maps into a final bad pixel map
    do_not_use = conform_check(do_not_use)
    final_map = combine_individual_maps(stack_of_maps, do_not_use)

    # Save the bad pixel map in DMS format
    if output_file is None:
        output_file = os.path.join(os.getcwd(), '{}_{}_mask.fits'.format(instrument, detector))
    save_final_map(final_map, input_files, sigma_threshold, smoothing_box_width, dead_sigma_threshold,
                   max_dead_norm_signal, max_low_qe_norm_signal, max_open_adj_norm_signal,
                   output_file)

    # Basic compatibility check: run the pipeline on a dummy file using
    # this reference file. Make sure the pipeline doesn't crash
    if quality_check:
        pipeline_check(output_file)


def combine_individual_maps(bad_maps, do_not_use_flags):
    """Given multiple separate bad pixel maps, combine using bit
    values matching those defined in the JWST calibraiton pipeline

    Parameters
    ----------
    bad_maps : dict
        Dictionary containing the various individual bad pixel maps to be
        combined. Format is (e.g.) {'DEAD': dead_map, 'LOW_QE': low_qe_map}

    Returns
    -------
    final_map : numpy.ndarray
        2D array containing the combined bad pixel map


    """
    # The official bit definitions for bad pixel flavors
    dq_defs = dqflags.pixel

    # Convert each type of bad pixel input to have the proper value,
    # and add to the final map
    final_map = np.zeros((bad_maps[].shape))
    for pix_type in bad_maps:
        if pix_type in dq_defs.keys():
            map = bad_maps[pix_type]
            flagged = map != 0
            value = dq_defs[pix_type]
            if pix_type in do_not_use_flags:
                value += dq_defs['DO_NOT_USE']
            final_map[flagged] += value
        else:
            raise ValueError('Unrecognized bad pixel type: {}'.format(pix_type))

    return final_map


def conform_check(flags):
    """Check that the flags in the input list are valid JWST bad pixel types

    Parameters
    ----------
    flag : list
        List of flag types (e.g. ['DEAD', 'LOW_QE'])

    Returns
    -------
    flag : list
        Modified list (all flags made uppercase)
    """
    possible_flags = dqflags.pixel
    flags = [entry.upper() for entry in flags]
    for flag in flags:
        if flag not in possible_flags:
            raise ValueError('Unrecognized flag type: {}'.format(flag))
    return flags


def crds_input_checks(author_name, descrip, pedigree, use_after):
    """Perform some basic checks on the input values to be placed in the
    header of the reference file. These check against CRDS requirements

    Parameters
    ----------

    Results
    -------
    """
    if len(descrip) > 65:
        raise ValueError('Description cannot exceed 65 characters.')

    allowed_pedigrees = ['GROUND']
    if pedigree not in allowed_pedigrees:
        raise ValueError('Pedigree must be one of: {}'.format(allowed_pedigrees))


def create_dummy_hdu_list(sigma_threshold, smoothing_width, dead_sigma_threshold, dead_max_rate,
                          low_qe_max_rate, open_adj_max_rate):
    """Create a fits HDU List to contain the "extra" fits header keywords
    that specific to this module and not included in the datamodels. This
    will be used to initialize the datamodel when creating the reference
    file

    Parameters
    ----------
    sigma_threshold : float
        Number of standard deviations to use when sigma-clipping to
        calculate the mean slope image or the mean across the detector

    smoothing_width : float
        Width in pixels of the box kernel to use to compute the smoothed
        mean image

    dead_sigma_threshold : float
        Number of standard deviations below the mean at which a pixel is
        considered dead.

    dead_max_rate : float
        Maximum normalized signal rate of a pixel that is considered dead

    low_qe_max_rate: float
        The maximum normalized signal a pixel can have and be considered
        low QE.

    open_adj_max_rate: float
        The maximum normalized signal a pixel adjacent to a low QE pixel can have
        in order for the low QE pixel to be reclassified as OPEN


    Returns
    -------
    hdu_list : astropy.io.fits.HDUList
        HDU List containing a header with just the extra keywords specific
        to this module
    """
    sig_thresh_keyword = 'BPMSIGMA'
    smooth_keyword = 'BPMSMOTH'
    dead_sigma_keyword = 'BPMDEDSG'
    max_dead_keyword = 'BPMMXDED'
    max_low_qe_keyword = 'BPMMXLQE'
    max_open_adj_keyword = 'BPMMXOAD'
    hdu = fits.PrimaryHDU()
    hdu[0].header[sig_thresh_keyword] = sigma_threshold
    hdu[0].header[smooth_keyword] = smoothing_width
    hdu[0].header[dead_sigma_keyword] = dead_sigma_threshold
    hdu[0].header[max_dead_keyword] = dead_mx_rate
    hdu[0].header[max_low_qe_keyword] = low_qe_max_rate
    hdu[0].header[max_open_adj_keyword] = open_adj_max_rate
    hdu_list = fits.HDUList([hdu])
    return hdu_list


def extract_10th_group(data):
    """Keep only the 10th group from each integration

    Parameters
    ----------
    data : numpy.ndarray
        3D or 4D array of data

    Returns
    -------
    data : numpy.ndarray
        3D array (10th group only)
    """
    dims = data.shape
    if len(dims) == 4:
        group10 = data[:, 9, :, :]
    elif len(dims) == 3:
        group10 = np.expand_dims(data[9, :, :], axis=0)
    return group10


def find_dead_pixels(rate_image, mean_rate, stdev_rate, instrument, sigma=5., max_dead_signal=0.05):
    """Create a map of dead pixels given a normalized rate image

    Parameters
    ----------
    rate_image : numpy.ndarray
        2D normalized rate image

    mean_rate : float
        Sigma-clipped mean value of the ``rate_image``

    stdev_rate : float
        Sigma-clipped standard deviation of the ``rate_image``

    instrument : str
        Name of the instrument used in the observations

    sigma : float
        Number of standard deviations below the mean at which a pixel is
        considered dead.

    max_dead_signal : float
        Maximum normalized signal rate of a pixel that is considered dead

    Returns
    -------
    dead_pix_map : numpy.ndarray
        2D map showing DEAD pixels. Good pixels have a value of 0.

    low_qe_map : numpy.ndarray
        2D map showing LOW QE pixels. Good pixels have a value of 0.
    """
    instrument = instrument.lower()
    if instrument in ['nircam', 'niriss']:
        dead_pix_map = (rate_image < (mean_rate - sigma * stdev_rate)).astype(np.int)
    elif instrument == 'nirspec':
        dead_pix_map = (rate_image < max_dead_signal).astype(np.int)
    return dead_pix_map


def find_miri_dead_pixels(groups, min_zero_signal_fraction=0.9):
    """Create a map of dead pixels given a set of individual groups
    from multiple integrations

    Parameters
    ----------
    groups : numpy.ndarray
        3D stack of group images

    min_zero_signal_fraction : float
        Threshold for the fraction of groups in which a pixel can have
        zero signal and not be considered dead.

    Returns
    -------
    dead_pix_map : numpy.ndarray
        2D map showing DEAD pixels. Good pixels have a value of 0.

    low_qe_map : numpy.ndarray
        2D map showing LOW QE pixels. Good pixels have a value of 0.
    """
    num_groups, ydim, xdim = groups.shape
    zero_signal = (groups == 0.).astype(np.int)
    total_zeros = np.sum(zero_signal, axis=0)
    total_fraction = total_zeros / num_groups
    dead_pix_map = (total_fraction >= min_zero_signal_fraction).astype(np.int)
    low_qe_map = np.zeros((ydim, xdim))
    return dead_pix_map, low_qe_map


def find_open_and_low_qepixels(rate_image, max_dead_signal=0.05, max_low_qe=0.5, max_adj_open=1.05):
    """Create maps of open (and adjacent to open) and low QE pixels given a
    normalized rate image

    Parameters
    ----------
    rate_image : numpy.ndarray
        2D normalized rate image

    mean_rate : float
        Sigma-clipped mean value of the ``rate_image``

    stdev_rate : float
        Sigma-clipped standard deviation of the ``rate_image``

    max_dead_signal : float
        The maximum normalized signal a pixel can have and be considered dead

    max_low_qe : float
        The maximum normalized signal a pixel can have and be considered
        low QE.

    max_adj_open : float
        The maximum normalized signal a pixel adjacent to a low QE pixel can have
        in order for the low QE pixel to be reclassified as OPEN

    Returns
    -------
    open_pix_map : numpy.ndarray
        2D map showing OPEN (1) and ADJ_OPEN (2) pixels. Good pixels have a
        value of 0.

    low_qe_map : numpy.ndarray
        2D map showing LOW QE pixels. Good pixels have a value of 0.
    """
    low_qe_map = np.zeros(rate_image.shape)
    open_pix_map = np.zeros(rate_image.shape)
    adj_pix_map = np.zeros(rate_image.shape)
    low_sig_x, low_sig_y = np.where((rate_image >= max_dead_signal) & (rate_image < max_low_qe))
    for x, y in zip(low_sig_x, low_sig_y):
        adj_pix_x = np.array([x, x+1, x, x-1])
        adj_pix_y = np.array([y+1, y, y-1, y])
        adj_pix = rate_image[adj_pix_y, adj_pix_x]
        adj_check = adj_pix > max_adj_open
        if all(adj_check):
            adj_pix_map[y-1:y+2, x-1:x+2] = 2
            open_pix_map[y, x] = 1
        else:
            low_qe_map[y, x] = 1
    return low_qe_map, open_pix_map, adj_pix_map


def get_fastaxis(filename):
    """Get the fastaxis and slowaxis values from the input file

    Parameters
    ----------
    filename : str
        Name of fits file to get values from

    Returns
    -------
    fastaxis : int
        Value of FASTAXIOS keyword

    slowaxis : int
        Value of SLOWAXIS keyword
    """
    with fits.open(filename) ash hdu_list:
        try:
            fastaxis = hdu_list[0].header['FASTAXIS']
            slowaxis = hdu_list[0].header['SLOWAXIS']
        except (KeyError, FileNotFoundError) as e:
            print(e)
    return fastaxis, slowaxis


def image_stats(image, sigma=3.):
    """Calculate the sigma-clipped mean and standard deviation across the
    given image

    Parameters
    ----------
    image : numpy.ndarray
        2D image

    sigma : float
        Sigma threshold to use for sigma-clipping

    Returns
    -------
    mean_val : float
        Sigma-clipped mean value

    stdev_val : float
        Sigma-clipped standard deviation
    """
    clipped_image = sigma_clip(image, sigma=sigma, masked=False)
    mean_val = np.nanmean(clipped_image)
    stdev_val = np.nanstd(clipped_image)
    return mean_val, stdev_val


def mean_stdev_images(data, sigma=3.):
    """Calculate the sigma-clipped mean and stdev images from a stack of
    input images

    Parameters
    ----------
    data : numpy.ndarray
        Stack of 2D images (i.e 3D array)

    sigma : float
        Threshold value to use for sigma clipping.

    Returns
    -------
    mean_image : numpy.ndarray
        2D array of the mean image

    stdev_image : numpy.ndarray
        2D array of the stdev image
    """
    clipped_cube = sigma_clip(data, sigma=sigma, axis=0, masked=False)
    mean_image = np.nanmean(clipped_cube, axis=0)
    stdev_image = np.nanstd(clipped_cube, axis=0)
    return mean_image, stdev_image


def miri_bad_columns(dimensions):
    """Create a map that flags the shorted columns on the MIRI detector

    Parameters
    ----------
    dimensions : tup
        (y, x) dimensions, in pixels, of the map to create

    Returns
    -------
    shorted_map : numpy.ndarray
        2D map showing the locations of the shorted columns (1)
    """
    shorted_map = np.zeros(dimensions).astype(np.int)
    shorted_map[:, 385:387] = 1
    return shorted_map


def pipeline_check(filename):
    """Run the dq_init step using the provided file, to check that the
    pipeline runs successfully

    Parameters
    ----------
    filename : str
        Name of the bad pixel mask reference filename
    """
    dummy_model = Level1bModel()
    dummy_model.data = np.zeros(())
    DQInitStep.call(dummy_model, override_mask=filename)


def read_files(filenames):
    """Read the data in the input files. Perform basic sanity checks
    to be sure data are consistent (same instrument, aperture, etc)

    Parameters
    ----------
    filenames : list
        List of fits filenames to be opened

    Returns
    -------
    data : numpy.ndarray
        3D stack of data
    """
    for i, filename in enumerate(filenames):
        with fits.open(filename) as hdu_list:
            exposure = hdu_list['SCI'].data
            insturment = hdu_list[0].header['INSTRUME'].upper()
            detector = hdu_list[0].header['DETECTOR'].upper()
            aperture = hdu_list[0].header['SUBARRAY'].upper()

        if instrument == 'MIRI':
            exposure = extract_10th_group(exposure)
        else:
            # Make 3D if it's not already
            dims = exposure.shape
            ndims = len(dims)
            if ndims == 2:
                exposure = np.expand_dims(exposure, axis=0)

        # Create the comparison cases and output array if we are
        # working on the first file
        if i == 0:
            comparison_instrument = copy.deepcopy(instrument)
            comparison_detector = copy.deepcopy(detector)
            comparison_aperture = copy.deepcopy(aperture)
            integrations = copy.deepcopy(exposure)

        # Consistency checks
        if instrument != comparison_instrument:
            raise ValueError('Inconsistent instruments in input data!')
        if detector != comparison_detector:
            raise ValueError('Inconsistent detectors in input data!')
        if aperture != comparison_aperture:
            raise ValueError('Inconsistent apertures in input data!')

        # Stack the new integrations onto the outuput
        integrations = np.concatenate((integrations, exposure), axis=0)
    return integrations, comparison_instrument, comparison_detector


def reference_pixel_map(dimensions):
    """Create a map that flags all reference pixels as such

    Parameters
    ----------
    dimensions : tup
        (y, x) dimensions, in pixels, of the map to create

    Returns
    -------
    ref_map : numpy.ndarray
        2D map showing the locations of reference pixels (1)

    """
    ref_map = np.zeros(dimensions).astype(np.int)
    ref_map[0:4, :] += 1
    ref_map[:, 0:4] += 1
    ref_map[yd-4:yd, :] += 1
    ref_map[:, xd-4:xd] += 1
    return ref_map


def save_final_map(bad_pix_map, instrument, detector, files, author, description, pedigree, useafter,
                   sigma_thresh, smooth_width, dead_sigma_thresh, max_dead_rate, max_low_qe_rate,
                   max_open_adj_rate, outfile):
    """Save a bad pixel map into a CRDS-formatted reference file

    Parameters
    ----------

    """
    # Define the non-standard fits header keywords by placing them in a
    # fits HDU List
    sig_thresh_keyword = 'BPMSIGMA'
    smooth_keyword = 'BPMSMOTH'
    dead_sigma_keyword = 'BPMDEDSG'
    max_dead_keyword = 'BPMMXDED'
    max_low_qe_keyword = 'BPMMXLQE'
    max_open_adj_keyword = 'BPMMXOAD'
    hdu = fits.PrimaryHDU()
    hdu[0].header[sig_thresh_keyword] = sigma_threshold
    hdu[0].header[smooth_keyword] = smooth_width
    hdu[0].header[dead_sigma_keyword] = dead_sigma_thresh
    hdu[0].header[max_dead_keyword] = max_dead_rate
    hdu[0].header[max_low_qe_keyword] = max_low_qe_rate
    hdu[0].header[max_open_adj_keyword] = max_open_adj_rate
    hdu_list = fits.HDUList([hdu])

    yd, xd = bad_pix_map.shape

    # Initialize the MaskModel using the hdu_list, so the new keywords will
    # be populated
    model = MaskModel(hdu_list)
    model.dq = np.uint8(bad_pix_map)
    model.dq_def = dqflags.pixel

    model.meta.reffile.type = 'MASK'
    model.meta.subarray.name = 'FULL'
    model.meta.subarray.xstart = 1
    model.meta.subarray.xsize = xd
    model.meta.subarray.ystart = 1
    model.meta.subarray.ysize = yd
    model.meta.instrument.name = instrument.upper()
    model.meta.instrument.detector = detector

    # Get the fast and slow axis directions from one of the input files
    fastaxis, slowaxis = get_fastaxis(files[0])
    model.meta.subarray.fastaxis = fastaxis
    model.meta.subarray.slowaxis = slowaxis

    model.meta.reffile.author = author
    model.meta.reffile.description = description
    model.meta.reffile.pedigree = pedigree
    model.meta.reffile.useafter = useafter

    # Populate "extra" header keywords that will contain parameters used
    # in this module

    package_note = ('This file was created using the bad_pixel_mask.py module within the '
                    'jwst_reffiles package.')

    software_dict = {'name': 'jwst_reffiles.bad_pixel_mask.py', 'author': 'STScI',
                     'homepage': 'https://github.com/spacetelescope/jwst_reffiles',
                     'version': 0.0.0}
    entry = util.create_history_entry(package_note, software=software_dict)
    model.history.append(entry)

    model.history.append(util.create_history_entry('Parameter values and descriptions:'))
    sigma_descrip = ('sigma_thresh: Number of standard deviations to use when sigma-clipping to '
                     'calculate the mean slope image or the mean across the detector. The value '
                     'used is stored in the {} keyword.'.format(sig_thresh_keyword))
    model.history.append(util.create_history_entry(sigma_descrip))

    smooth_descrip = ('smoothing_box_width: Width in pixels of the box kernel to use to compute the '
                      'smoothed mean image. The value used is stored in the {} keyword.'.format(smooth_keyword))
    model.history.append(util.create_history_entry(smooth_descrip))

    dead_sig_descrip = ('Number of standard deviations below the mean at which a pixel is considered dead. '
                        'The value used is stored in the {} keyword.'.format(dead_sigma_keyword))
    model.history.append(util.create_history_entry(dead_sig_descrip))

    max_dead_descrip = ('Maximum normalized signal rate of a pixel that is considered dead. The value '
                        'used is stored in the {} keyword.'.format(max_dead_keyword))
    model.history.append(util.create_history_entry(mx_dead_descrip))

    max_low_qe_descrip = ('The maximum normalized signal a pixel can have and be considered low QE. The '
                          'value used is stored in the {} keyword.'.format(max_low_qe_keyword))
    model.history.append(util.create_history_entry(max_low_qe_descrip))

    max_open_adj_descrip = ('The maximum normalized signal a pixel adjacent to a low QE pixel can have '
                            'in order for the low QE pixel to be reclassified as OPEN. The value used '
                            'is stored in the {} keyword.'.format(max_open_adj_keyword))
    model.history.append(util.create_history_entry(max_open_adj_descrip))

    do_not_use_descrip = ('List of bad pixel types where the DO_NOT_USE flag should also be applied. '
                          'Values used are: {}'.format(do_not_use_list))
    model.history.append(util.create_history_entry(do_not_use_descrip))

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

    if history_text is not None:
        model.history.append(util.create_history_entry(history_text))

    model.save(outfile, overwrite=True)
    print('Final bad pixel mask reference file save to: {}'.format(outfile))


def smooth(data, box_width=15):
    """Create a smoothed version of the 2D input data

    Parameters
    ----------
    data : numpy.ndarray
        2D array of data to be smoothed

    box_width : int
        Width of the smoothing box, in pixels

    Returns
    -------
    smoothed : numpy.ndarray
        A smoothed version of ``data``
    """
    smoothing_kernel = Box2DKernel(box_width)
    smoothed = convolve(data, smoothing_kernel)
    return smoothed


if __name__ == '__main__':
    usagestring = ('USAGE: bad_pixel_mask.py input_files sigma_threshold=3 smoothing_box_width=15 '
                   'dead_sigma_threshold=5. max_dead_norm_signal=0.05, max_low_qe_norm_signal=0.5, '
                   'max_open_adj_norm_signal=1.05, do_not_use=[], output_file=None')

    gainim = gainimclass()
    parser = gainim.add_options(usage=usagestring)
    gainim.options, args = parser.parse_args()
    gainim.doit()







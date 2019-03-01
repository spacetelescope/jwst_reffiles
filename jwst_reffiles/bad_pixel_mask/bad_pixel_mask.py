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
import numpy as np

from astropy.convolution import convolve, Box2DKernel
from astropy.io import fits
from scipy.stats import sigmaclip


def find_bad_pix(input_files, sigma_threshold=3, smoothing_box_width=15, dead_sigma_threshold=5.,
                 max_dead_norm_signal=0.05, max_low_qe_norm_signal=0.5, max_open_adj_norm_signal=1.05):
    """MAIN SCRIPT: Given a set of input files, create deadd, low QE, and open pixel maps

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
    lowqe_map, open_map = find_open_and_low_qe_pixels(normalized, max_dead_signal=max_dead_norm_signal,
                                                      max_low_qe=max_low_qe_norm_signal,
                                                      max_adj_open=max_open_adj_norm_signal)

    return dead_map, lowqe_map, open_map


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
    low_sig_x, low_sig_y = np.where((rate_image >= max_dead_signal) & (rate_image < max_low_qe))
    for x, y in zip(low_sig_x, low_sig_y):
        adj_pix_x = np.array([x, x+1, x, x-1])
        adj_pix_y = np.array([y+1, y, y-1, y])
        adj_pix = rate_image[adj_pix_y, adj_pix_x]
        adj_check = adj_pix > max_adj_open
        if all(adj_check):
            open_pix_map[y-1:y+2, x-1:x+2] = 2
            open_pix_map[y, x] = 1
        else:
            low_qe_map[y, x] = 1
    return low_qe_map, open_pix_map


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

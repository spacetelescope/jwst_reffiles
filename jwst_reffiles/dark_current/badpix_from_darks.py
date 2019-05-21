#! /usr/bin/env python

"""Find bad pixels from dark current files


Start with a stack of dark ramps and slope images?
Or maybe a stack of ramps that have been processed through the jump step
and then ramp-fitting is performed here?

Input for the dark current reference file step is going to be a stack
of ramps. So maybe use that and ramp-fit here.


inputs:
1. list of dark current ramps that have been run through jump step
2. list of same exposures after ramp-fitting has been done

Plot summary:

0. Check to see if IPC correction has been run
1. Look through stack of slopes images, get mean and rms per pixel
   (do we sigma-clip the inputs or not?)
2. Potential bad pixels are those with noise values above some threshold


NOTE: when producing slope images of these data, make sure to save the
optional output parameters into the *fitopt.fits files.
https://jwst-pipeline.readthedocs.io/en/stable/jwst/ramp_fitting/main.html?highlight=intercept
"A third, optional output product is also available and is produced only when the step parameter ‘save_opt’ is True"

"""
from astropy.io import fits
import copy
from jwst.datamodels import dqflags
import numpy as np
from scipy.stats import sigmaclip


def get_cr_flags(dq_array):
    """Return a map of cosmic ray flags given a data quality array

    Parameters
    ----------
    dq_array : numpy.ndarray
        GROUP_DQ extension of an integration

    Returns
    -------
    cr_map : numpy.ndarray
        GROUP_DQ extension with all flags other than jumps removed
    """
    cr_map = (dq_array & dqflags.pixel['JUMP'] > 0)
    return cr_map


def find_pix_with_many_jumps(jump_map, sigma_threshold=5, early_cutoff_fraction=0.25):
    """Identify pixels that have an abnormal number of flagged jumps

    Parameters
    ----------
    jump_map : numpy.ndarray
        Map of jump flags for all pixels (e.g. output from ``get_cr_flags``)
        Assume we are working one integration at a time.

    sigma_threshold : int
        Number of sigma from the mean that defines the minimum number of
        jumps to classify a pixel as having an abnormally high number of
        jumps

    early_cutoff_fraction : float
        Fraction of the integration to use when comparing the jump rate
        early in the integration to that across the entire integration

    Returns
    -------

    """
    jump_map = jump_map.astype(np.int)
    number_of_jumps = np.sum(jump_map, axis=0)
    avg_number_of_jumps = np.median(number_of_jumps)
    stdev_number_of_jumps = np.std(number_of_jumps)

    high_jumps = number_of_jumps >= (avg_number_of_jumps + stdev_number_of_jumps*sigma_threshold)

    # Repeat using only the first N groups, and compare to the results
    # above as a way of finding pixels with many early jumps, which is
    # a sign of an RC or IRC pixel
    number_of_groups = jump_map.shape[0]
    early_cutoff = np.int(early_cutoff_fraction * number_of_groups)
    number_of_early_jumps = np.sum(jump_map[early_cutoff], axis=0)
    avg_number_of_early_jumps = np.median(number_of_early_jumps)
    stdev_number_of_early_jumps = np.std(number_of_early_jumps)

    # Compare avg_number_of_early_jumps with avg_number_of_jumps?

    # Compare the ratio of jump rate early in the ramp to that across
    # the entire ramp
    jump_ratio = number_of_early_jumps / early_cutoff_fraction / number_of_jumps
    potential_rc = jump_ratio >= (1./(2.*early_cutoff_fraction))
    print('there is probably a better way to look for these')

    return high_jumps, potential_rc


def read_pedestal_data(filename):
    """Read in the PEDESTAL values from a *fitopt.fits file

    Parameters
    ----------
    filename : str
        Name of output file to check. This should be a *fitopt.fits file.

    Returns
    -------
    pedestal : numpy.ndarray
        3D array of pedestal values (signal extrapolated to time=0)
    """
    with fits.open(filename) as hdulist:
        pedestal = hdulist['PEDESTAL'].data
    return pedestal


def saturated_in_all_groups(pedestal_array):
    """Generate a list of pixels that are saturated in all groups

    Parameters
    ----------
    pedestal_array : numpy.ndarray
        3D array of pedestal values (signal extrapolated to time=0)

    Returns
    -------
    full_saturation : tup
        Tuple of (y, x) coordinate lists (output from np.where)
    """
    full_saturation = pedestal_array == 0
    return full_saturation


def pedestal_stats(pedestal_array, threshold=5):
    """Get statsitics on the pedestal array corresponding to an
    integration

    """
    median_pedestal = np.median(pedestal_array)
    stdev_pedestal = np.std(pedestal_array)

    suspicious_pedestal = ((pedestal_array > (median_pedestal + stdev_pedestal*threshold)) or
                           (pedestal_array < (median_pedestal - stdev_pedestal*threshold)))
    return suspicious_pedestal


def check_metadata(hdr, comp):
    """Basic metadata check for consistency from one file to another
    """
    keywords = ['EFFEXPTM', 'TFRAME', 'TGROUP', 'NFRAMES', 'NGROUPS', 'SUBARRAY']
    file_name = hdr['FILENAME']
    compare_name = comp['FILENAME']
    for key in keywords:
        value = hdr[key]
        compare_value = comp[key]
        if isinstance(value, str) or isinstance(value, int):
            if value != compare_value:
                raise ValueError('Inconsistent input files. {} is different between {} and {}.'
                                 .format(key, file_name, compare_name))
        elif isinstance(value, float):
            if not np.isclose(value, compare_value, rtol=0, atol=0.001):
                raise ValueError('Inconsistent input files. {} is different between {} and {}.'
                                 .format(key, file_name, compare_name))


def read_slope_files(filenames):
    """Read in the science extension from a group of slope images
    """
    for i, filename in enumerate(filenames):
        # Read all of the slope data into an array
        slope_file = filename.replace('jump.fits', 'rate.fits')
        with fits.open(slope_file) as hdulist:
            slope_img = hdulist['SCI'].data
            header = hdulist['SCI'].header
        if i == 0:
            ydim, xdim = slope_img.shape
            slope_data = np.zeros((len(filenames), ydim, xdim))
            header_to_compare = copy.deepcopy(header)
        else:
            # Check to be sure the input files are consistent
            check_metadata(header, header_to_compare)

        slope_data[i, :, :] = slope_img
    return slope_data


def find_bad_pix(filenames, noisy_threshold=5, saturation_threshold=0.5):
    """MAIN FUNCTION
    Assume that pipeline products all exist.
    Assume input filenames are those of the ramp files
    (which have been processed
    through the JUMP step)
    """
    # Read in the slope data
    slopes = read_slope_files(filenames)

    # Calculate the mean and standard deviation through the stack for
    # each pixel. Assuming that we are looking for noisy pixels, we don't
    # want to do any sigma clipping on the inputs here, right?
    mean_slope = np.mean(slopes, axis=0)
    std_slope = np.std(slopes, axis=0)
    avg_of_std = np.mean(std_slope)
    std_of_std = np.std(std_slope)

    # Identify noisy pixels as those with noise values more than
    # noisy_threshold*sigma above the average noise level
    noisy = std_slope > (avg_of_std + std_of_std*noisy_threshold)

    # Read in the optional outputs from the ramp-fitting step, so that
    # we can look at the y-intercepts and the jump flags
    saturated = np.zeros(slopes.shape)
    rc_from_pedestal = np.zeros(slopes.shape)
    low_pedestal = np.zeros(slopes.shape)
    high_cr_rate = np.zeros(slopes.shape)
    rc_from_flags = np.zeros(slopes.shape)

    for i, filename in enumerate(filenames):
        # Read in the fitops file associated with the exposure and get
        # the pedestal array
        pedestal_file = filename.replace('rate.fits', 'fitops.fits')
        pedestal = read_pedestal_data(pedestal_file)

        clipped_pedestal = sigmaclip(pedestal, low=3., hgih=3.)
        mean_pedestal = np.mean(clipped_pedestal)
        std_pedestal = np.std(clipped_pedestal)
        rc_from_pedestal[i, :, :] = pedestal > (mean_pedestal + std_pedestal*pedestal_threshold)

        # Are these useful to track?
        low_pedestal[i, :, :] = pedestal < (mean_pedestal + std_pedestal*pedestal_threshold)

        # Find pixels that are saturated in all groups. These will have
        # a pedestal value of 0 (according to the pipeline documentation).
        # These should end up flagged as HOT and DO_NOT_USE
        saturated[i, :, :] = saturated_in_all_groups(pedestal)

        # Read in the ramp and get the data and dq arrays
        ramp_file = filename.replace('rate.fits', 'jump.fits')
        with fits.open(ramp_file) as hdulist:
            data = hdulist['SCI'].data
            groupdq = hdulist['GROUPDQ'].data

        # Generate a map of JUMP flags in the ramp
        cr_map = get_cr_flags(groupdq)

        # Find pixels that have an abnormally high number of jumps, as
        # well as those that have most of their jumps concentrated in the
        # early part of the integration. The latter are possibly RC or IRC
        # pixels
        many_jumps, rc_candidates = find_pix_with_many_jumps(cr_map, sigma_threshold=5,
                                                             early_cutoff_fraction=0.25)
        high_cr_rate[i, :, :] = copy.deepcopy(many_jumps)
        rc_from_flags[i, :, :] = copy.deepcopy(rc_candidates)

    # Look through the stack of saturated pixels and keep those saturated
    # more than N% of the time
    fully_saturated = np.sum(saturated, axis=0) / len(filenames)
    fully_saturated[full_saturated < saturation_threshold] = 0
    fully_saturated = np.ceil(fully_saturated).astype(np.int)

    # How do we want to combine these to identify RC pixels?
    rc_pedestal = np.sum(rc_from_pedestal, axis=0) / len(filenames)
    rc_flags = np.sum(rc_from_flags, axis=0) / len(filenames)
    rc = (rc_pedestal > rc_threshold) or (rc_flags > rc_threshold)

    # Pixels with lots of CR flags should be added to the list of noisy pixels?
    high_cr = np.sum(high_cr_rate, axis=0) / len(filenames)
    noisy_second_pass = high_cr > high_cr_threshold
    combined_noisy = np.bitwise_or(noisy, noisy_second_pass)










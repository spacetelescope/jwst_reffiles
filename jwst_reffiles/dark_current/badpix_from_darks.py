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
from astropy.stats import sigma_clip
import copy
import os
from jwst.datamodels import dqflags
import numpy as np
from os import path
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip
import matplotlib.cm as cm
from jwst_reffiles.bad_pixel_mask.bad_pixel_mask import create_dqdef
from jwst_reffiles.utils import dq_flags
from jwst_reffiles.utils.constants import RATE_FILE_SUFFIXES


def find_bad_pix(filenames, clipping_sigma=5., max_clipping_iters=5, noisy_threshold=5,
                 max_saturated_fraction=0.5,
                 max_jump_limit=10, jump_ratio_threshold=5, early_cutoff_fraction=0.25,
                 pedestal_sigma_threshold=5, rc_fraction_threshold=0.8, low_pedestal_fraction=0.8,
                 high_cr_fraction=0.8,
                 flag_values={'hot': ['HOT'], 'rc': ['RC'], 'low_pedestal': ['OTHER_BAD_PIXEL'], 'high_cr': ["TELEGRAPH"]},
                 do_not_use=['hot', 'rc', 'low_pedestal', 'high_cr'], outfile=None, plot=False):
    """MAIN FUNCTION

    Parameters
    ----------
    filenames : list
        List of dark current files. These should be slope images.

    clipping_sigma : int
        Number of sigma to use when sigma-clipping the 2D array of
        standard deviation values. The sigma-clipped mean and standard
        deviation are used to locate noisy pixels.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma clipping to find
        the mean and standard deviation values that are used when
        locating noisy pixels.

    noisy_threshold : int
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

    do_not_use : list
        List of bad pixel types to be flagged as DO_NOT_USE
        e.g. ['hot', 'rc', 'low_pedestal', 'high_cr']

    outfile : str
        Name of fits file to save the resulting bad pixel mask to
    """
    # Add DO_NOT_USE to all requested types of bad pixels
    do_not_use = [element.lower() for element in do_not_use]
    for key in flag_values:
        if key.lower() in do_not_use:
            flag_values[key].append('DO_NOT_USE')

    # Form the outfile and outdir
    if outfile is None:
        outfile = 'badpixels_from_darks.fits'

    outdir = os.path.dirname(outfile)
    if not outdir:
        outdir = '.'

    # Read in the slope data. Strip off reference pixels.
    # Return a 3D array of slopes and a 3D array mapping where the
    # science pixels are.
    print('Reading slope files...')
#    instrument,slopes, refpix_additions = read_slope_files(filenames)

    instrument, slopes, indexes, refpix_additions = read_slope_integrations(filenames)

    shape_slope = slopes.shape
    print('Number of integrations used to flag bad pixels', shape_slope[0])
    # Calculate the mean and standard deviation through the stack for
    # each pixel. Assuming that we are looking for noisy pixels, we don't
    # want to do any sigma clipping on the inputs here, right?
    mean_slope = np.mean(slopes, axis=0)
    std_slope = np.std(slopes, axis=0)
    hdout = fits.PrimaryHDU(mean_slope)
    hdout.writeto('average_of_slopes.fits', overwrite=True)
    hdout = fits.PrimaryHDU(std_slope)
    hdout.writeto('sigma_of_slopes.fits', overwrite=True)

    # Use sigma-cliping when calculating the mean and standard deviation
    # of the standard deviations
    clipped_stdevs, cliplow, cliphigh = sigma_clip(std_slope, sigma=clipping_sigma,
                                                   maxiters=max_clipping_iters,
                                                   masked=False, return_bounds=True)

    avg_of_std = np.mean(clipped_stdevs)
    std_of_std = np.std(clipped_stdevs)
    cut_limit = avg_of_std + std_of_std*noisy_threshold

    # Identify noisy pixels as those with noise values more than
    # noisy_threshold*sigma above the average noise level
    # noisy = std_slope > cut_limit # not a good stat we need to remove slopes with cr hits
    # Plot histogram to later compare with better std_slope only containing
    # slopes with no jumps detected.
    if plot:
        xhigh = avg_of_std + std_of_std*noisy_threshold
        plot_image(std_slope, xhigh, outdir,
                   "Pixel Standard devations", "pixel_std_withjumps.png")

        nbins = 5000
        titleplot = 'Histogram of Pixel Slope STD with cosmic ray jumps: Clipped Ave ' + \
            '{:6.4f}'.format(avg_of_std) + '  Std ' + '{:6.4f}'.format(std_of_std)

        plot_histogram_stats(std_slope, cut_limit, nbins,
                             outdir, titleplot,
                             "histo_std_withjumps.png", xaxis_log=True)

    # Read in the optional outputs from the ramp-fitting step, so that
    # we can look at the y-intercepts and the jump flags

    saturated = np.zeros(slopes.shape)
    rc_from_pedestal = np.zeros(slopes.shape)
    low_pedestal = np.zeros(slopes.shape)
    high_cr_rate = np.zeros(slopes.shape)
    rc_from_flags = np.zeros(slopes.shape)
    slope_stack = []
    islope_stack = []

    total_ints = 0
    counter = 0

    for i, filename in enumerate(filenames):

        # Read in the ramp and get the data and dq arrays
        jump_file = None
        for suffix in RATE_FILE_SUFFIXES:
            if suffix in filename:
                slope_suffix = '{}.fits'.format(suffix)
                jump_file = filename.replace(slope_suffix, '_jump.fits')
                break
        if jump_file is None:
            raise ValueError("ERROR: Unrecognized slope filename suffix.")
        if not os.path.isfile(jump_file):
            raise FileNotFoundError("ERROR: Jump file {} not found.".format(jump_file))

        print('Opening Jump File {}'.format(jump_file))
        groupdq = dq_flags.get_groupdq(jump_file, refpix_additions)
        cr_map = dq_flags.flag_map(groupdq, 'JUMP_DET')

        # Get slope data corresponding to this file by extracting the
        # appropriate frames from the ``slopes`` stack
        slope = slopes[indexes[i]: indexes[i+1], :, :]

        # Read in the fitops file associated with the exposure and get
        # the pedestal array (y-intercept)
        pedestal_file = filename.replace(slope_suffix, '_fitopt.fits')
        if not os.path.isfile(pedestal_file):
            raise FileNotFoundError("ERROR: Pedestal file {} not found.".format(pedestal_file))
        print('Opening Pedestal File {}'.format(pedestal_file))
        pedestal = read_pedestal_data(pedestal_file, refpix_additions)

        # for MIRI the zero point of the ramp drifts with time. Adjust the
        # pedestal to be a relative pedestal wrt to group 2
        if instrument == 'MIRI':
            uncal_file = filename.replace(slope_suffix, '_uncal.fits')
            if not os.path.isfile(uncal_file):
                raise FileNotFoundError("ERROR: Uncal file {} not found.".format(uncal_file))
            group2 = extract_group2(uncal_file, refpix_additions)
            pedestal_org = copy.deepcopy(pedestal)
            pedestal = np.fabs(group2 - pedestal)

        # Work one integration at a time
        for int_num in range(pedestal.shape[0]):

            # pull out the DQ of the first group. This will be use to remove
            # Low pedestal values that have a pedestal of 0 because they are
            # saturated on group 1.
            first_group = groupdq[int_num, 0, :, :]
            pedestal_int = pedestal[int_num, :, :]
            slope_int = slope[int_num, :, :]

            clipped_pedestal, cliplow, cliphigh = sigmaclip(pedestal_int, low=3., high=3.)
            mean_pedestal = np.mean(clipped_pedestal)
            std_pedestal = np.std(clipped_pedestal)

            rc_from_pedestal[counter, :, :] += pedestal_int > (mean_pedestal + std_pedestal * pedestal_sigma_threshold)

            # Pixels with abnormally low pedestal values
            pedestal_low = pedestal_int < (mean_pedestal - std_pedestal * pedestal_sigma_threshold)
            first_group_sat = np.bitwise_and(first_group, dqflags.pixel['SATURATED'])

            # do not allow pixels saturated on group 1 to be marked as low pedestal
            pedestal_results = np.logical_and(pedestal_low, (first_group_sat == 0))
            low_pedestal[counter, :, :] += pedestal_results

            # Find pixels that are saturated in all groups. These will have
            # a pedestal value of 0 (according to the pipeline documentation).
            # These should end up flagged as HOT and DO_NOT_USE
            # Remove all the cases where ped = 0, but group 1 is not saturated
            # This can be dead pixels

            if instrument == 'MIRI':
                pedestal_int = pedestal_org[int_num, :, :]

            saturated[counter, :, :] += saturated_in_all_groups(pedestal_int, first_group_sat)

            # Find pixels that have an abnormally high number of jumps, as
            # well as those that have most of their jumps concentrated in the
            # early part of the integration. The latter are possibly RC or IRC
            # pixels
            many_jumps, rc_candidates, number_of_jumps =\
                find_pix_with_many_jumps(cr_map[int_num, :, :, :], max_jump_limit=10,
                                         jump_ratio_threshold=5,
                                         early_cutoff_fraction=0.25)

            high_cr_rate[counter, :, :] += many_jumps
            rc_from_flags[counter, :, :] += rc_candidates

            # using the number_of_jumps (a per integration value) create a clean set of
            # pixel slopes with no cosmic rays
            clean_slopes, iclean_slopes = slopes_not_cr(slope_int, number_of_jumps)
            slope_stack.append(clean_slopes)
            islope_stack.append(iclean_slopes)

            total_ints += 1
        counter += 1

    # now find the mean and standard deviation of the "clean" pixel slopes
    clean_mean_slope, clean_std_slope, num_good = combine_clean_slopes(slope_stack, islope_stack)
    hdout = fits.PrimaryHDU(clean_mean_slope)
    hdout.writeto('average_of_slopes_nojumps.fits', overwrite=True)
    hdout = fits.PrimaryHDU(clean_std_slope)
    hdout.writeto('sigma_of_slopes_nojumps.fits', overwrite=True)
    num_good_slopes = num_good.astype(np.int16)
    hdout = fits.PrimaryHDU(num_good_slopes)
    hdout.writeto('number_of_slopes_nojumps.fits', overwrite=True)

    # Use sigma-cliping to remove large outliers to have clean stats to flag
    # noisy pixels.
    # removing nans from clean_std_slope because it causes warning messages to be print
    clean_std_slope_nonan = clean_std_slope[np.isfinite(clean_std_slope)]

    clipped_stdevs, cliplow, cliphigh = sigma_clip(clean_std_slope_nonan, sigma=clipping_sigma,
                                                   maxiters=max_clipping_iters,
                                                   masked=False, return_bounds=True)

    avg_of_std = np.mean(clipped_stdevs)
    std_of_std = np.std(clipped_stdevs)
    cut_limit = avg_of_std + std_of_std*noisy_threshold

    # assigning nans from clean_std_slope to very large values that will be cut
    # because it causes warning messages to be print
    values_nan = np.isnan(clean_std_slope)
    clean_std_slope[values_nan] = avg_of_std + std_of_std*50

    noisy = clean_std_slope > cut_limit
    num_noisy = len(np.where(noisy)[0])

    if plot:
        # plot the number of good slopes per pixel
        max_values = np.amax(num_good)
        plot_image(num_good, max_values, outdir,
                   "Number of Good slopes/pixel ",
                   "clean_pixel_number.png")

        # plot the standard deviation of pixels slope after eliminating
        # values having jumps detectect in ramp
        xhigh = avg_of_std + std_of_std
        plot_image(clean_std_slope, xhigh, outdir,
                   "Clean Pixel Standard devations",
                   "clean_pixel_std.png")

        # plot the histogram before the clipping
        nbins = 5000
        titleplot = 'Histogram of Clean Pixel Slope STD  Average ' + \
            '{:6.4f}'.format(avg_of_std) + '  Std ' + '{:6.4f}'.format(std_of_std)

        plot_histogram_stats(clean_std_slope, cut_limit, nbins, outdir,
                             titleplot, "histo_clean_std.png", xaxis_log=True)

    # Look through the stack of saturated pixels and keep those saturated
    # more than N% of the time

    fully_saturated = np.sum(saturated, axis=0) / total_ints
    fully_saturated[fully_saturated < max_saturated_fraction] = 0
    fully_saturated = np.ceil(fully_saturated).astype(np.int)

    fully_saturated = apply_flags(fully_saturated, flag_values['hot'])
    num_saturated = len(np.where(fully_saturated != 0)[0])
    print('\n\nFound {} fully saturated pixels.'.format(num_saturated))

    # How do we want to combine these to identify RC pixels?
    rc_pedestal = np.sum(rc_from_pedestal, axis=0) / total_ints
    rc_flags = np.sum(rc_from_flags, axis=0) / total_ints

    rc_from_pedestal_only = (rc_pedestal > rc_fraction_threshold).astype(np.int)
    rc_from_jumps_only = (rc_flags > rc_fraction_threshold).astype(np.int)
    num_rc_ped = len(np.where(rc_from_pedestal_only != 0)[0])
    num_rc_jump = len(np.where(rc_from_jumps_only != 0)[0])
    print("Found {} RC pixels from pedestal search".format(num_rc_ped))
    print("Found {} RC pixels from Jump search".format(num_rc_jump))

    rc = ((rc_pedestal > rc_fraction_threshold) | (rc_flags > rc_fraction_threshold))
    rc = apply_flags(rc.astype(np.int), flag_values['rc'])
    num_rc = len(np.where(rc != 0)[0])
    print('Found {} RC pixels.'.format(num_rc))

    # Low pedestal pixels
    low_pedestal_vals = np.sum(low_pedestal, axis=0) / total_ints
    low_ped = low_pedestal_vals > low_pedestal_fraction

    # Pixels that are saturated on the first group will have a PEDESTAL value
    # of 0. Pull these out of this set (these are hot pixels)
    low_ped = apply_flags(low_ped.astype(np.int), flag_values['low_pedestal'])
    num_low_ped = len(np.where(low_ped != 0)[0])
    print('Found {} low pedestal pixels.'.format(num_low_ped))

    # Pixels with lots of CR flags should be added to the list of noisy pixels?
    high_cr = np.sum(high_cr_rate, axis=0) / total_ints
    noisy_second_pass = high_cr > high_cr_fraction
    combined_noisy = np.bitwise_or(noisy, noisy_second_pass)
    combined_noisy = apply_flags(combined_noisy.astype(np.int), flag_values['high_cr'])

    num_high_cr = len(np.where(noisy_second_pass != 0)[0])
    print('Found {} pixels with a high number of jumps.'.format(num_high_cr))
    print('Found {} pixels with noise above the threshold.'.format(num_noisy))
    num_combined_noisy = len(np.where(combined_noisy != 0)[0])
    print('Combining noisy and high jump pixels, found {} noisy pixels.'.format(num_combined_noisy))

    # Combine the various flavors of bad pixels into a final DQ map
    bad_pixels = combine_bad_pixel_types(fully_saturated, rc, low_ped, combined_noisy)

    # Add the reference pixels back into the bad pixel map
    bad_pixels = add_refpix(bad_pixels, refpix_additions)

    # Create DQ definitions to be saved with the output file
    dq_def = create_dqdef()

    # Save the bad pixel mask to a fits file
    # Eventually this routine will be called as part of the dark current reference file
    # generator, and the bad pixel mask will be saved in the DQ extension of the
    # reference file

    h0 = fits.PrimaryHDU(fully_saturated)
    h0.header['EXTNAME'] = 'SATURATED'
    h1a = fits.ImageHDU(rc_from_pedestal_only)
    h1a.header['EXTNAME'] = 'RC_FROM_PED'
    h1b = fits.ImageHDU(rc_from_jumps_only)
    h1b.header['EXTNAME'] = 'RC_FROM_JUMPS'
    h1 = fits.ImageHDU(rc)
    h1.header['EXTNAME'] = 'RC'
    h2 = fits.ImageHDU(low_ped)
    h2.header['EXTNAME'] = 'LOW_PEDESTAL'
    h3 = fits.ImageHDU(noisy.astype(np.int))
    h3.header['EXTNAME'] = 'NOISY'
    h4 = fits.ImageHDU(noisy_second_pass.astype(np.int))
    h4.header['EXTNAME'] = 'MANY_CRS'
    h5 = fits.ImageHDU(combined_noisy)
    h5.header['EXTNAME'] = 'NOISY_AND_CRS'
    hlist = fits.HDUList([h0, h1a, h1b, h1, h2, h3, h4, h5])
    hlist.writeto(outfile, overwrite=True)
    print('Multi-extension file with individual types of bad pixels saved to:')
    print(outfile)


def add_refpix(array, to_add):
    """Place ``map`` within a larger array that contains the reference
    pixels.

    Parameters
    ----------
    array : numpy.ndarray
        2D array of bad pixels that does not contain reference pixels

    to_add : tup
        4-element tuple containing the number of rows and columns to
        add around the outside of the science pixels.
        (left cols, right cols, bottom rows, top rows)

    Returns
    -------
    array : numpy.ndarray
        2D array with rows and columns added
    """
    left_cols, right_cols, bottom_rows, top_rows = to_add
    y_array, x_array = array.shape
    xdim = x_array + left_cols + right_cols
    ydim = y_array + bottom_rows + top_rows
    full_array = np.zeros((ydim, xdim))

    full_array[bottom_rows: bottom_rows+y_array, left_cols: left_cols+x_array] = array
    return array


def apply_flags(pixmap, flag_list):
    """Beginning with a map indicating locations of a particular type of
    bad pixel, apply the bits specified in ``flag_list`` to come up with
    the ``jwst`` bad pixel value.

    Parameters
    ----------
    pixmap : numpy.ndarray
        2D array indicating bad pixels. 1 for a bad pixel and 0 for a
        good pixel

    flag_list : list
        List of bad pixel mnemonics to be applied. These mnemonics must
        be in the dictionary of bad pixel types recognized by the jwst
        calibration pipeline

    Returns
    -------
    pixmap : numpy.ndarray
        Map updated with proper bad pixel values
    """
    for mnemonic in flag_list:
        if mnemonic in dqflags.pixel.keys():
            value = dqflags.pixel[mnemonic.upper()]
            pixmap[pixmap != 0] += value
        else:
            raise ValueError("ERROR: unrecognized DQ mnemonic: {}".format(mnemonic))
    pixmap = pixmap.astype(np.uint8)
    return pixmap


def check_metadata(hdr, comp):
    """Basic metadata check for consistency from one file to another

    Parameters
    ----------
    hdr : astropy.fits.header
        Header read in from fits file

    comp : astropy.fits.header
        Baseline header to which the comparison is done
    """
    keywords = ['EFFEXPTM', 'TFRAME', 'TGROUP', 'NFRAMES', 'NGROUPS', 'SUBARRAY']

    for key in hdr:
        print(key, hdr[key])

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


def combine_bad_pixel_types(sat_map, rc_map, low_pedestal_map, high_cr_map):
    """Copmbine individual maps of bad pixel types into a final bad pixel
    map, using flag values defined in ``dq_flags``.

    Parameters
    ----------
    sat_map : numpy.ndarray
        2D array giving the location of pixels saturated all the way up
        the ramp

    rc_map : numpy.ndarray
        2D array giving the location of RC pixels

    low_pededtal_map : numpy.ndarray
        2D array giving the location of pixels with abnormally low
        pedestal values

    high_cr_map : numpy.ndarray
        2D array giving the location of pixels with abnormally high
        numbers of jumps

    Returns
    -------
    final_map : numpy.ndarray
        2D array containing the bitwise combined bad pixel maps
    """

    sat_and_rc = np.bitwise_or(sat_map, rc_map)
    add_pedestal = np.bitwise_or(sat_and_rc, low_pedestal_map)
    final_map = np.bitwise_or(add_pedestal, high_cr_map)
    return final_map


def find_pix_with_many_jumps(jump_map, max_jump_limit=10, jump_ratio_threshold=5, early_cutoff_fraction=0.25):
    """Identify pixels that have an abnormal number of flagged jumps. Do
    this by finding the jump rate early in the ramp versus that later in
    the ramp.

    Parameters
    ----------
    jump_map : numpy.ndarray
        Map of jump flags for all pixels (e.g. output from ``dqflags.flag_map``)
        Assume we are working one integration at a time.

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

    Returns
    -------
    high_jumps : numpy.ndarray
        Map of pixels that have more than ``max_jump_limit`` jumps.

    potential_rc : numpy.ndarray
        Map of pixels which have: 1) a large number of jumps, and
        2) a higher rate of jumps early in the ramp than later in
        the ramp.
    """
    if early_cutoff_fraction > 0.5:
        raise ValueError("ERROR: early_cutoff_fraction must be less than or equal to 0.5")

    # First look across the entire ramp for pixels that have a large number
    # of jumps. Those with more than max_jump_limit will be flagged
    jump_map = jump_map.astype(np.int)
    number_of_jumps = np.sum(jump_map, axis=0)
    high_jumps = number_of_jumps >= max_jump_limit

    # Next compare the number of jumps early in the ramp to the number
    # later in the ramp. This is a way of finding pixels with many early
    # jumps, which is a sign of an RC or IRC pixel
    number_of_groups = jump_map.shape[0]
    early_cutoff = np.int(early_cutoff_fraction * number_of_groups)
    early_jump_rate = np.sum(jump_map[0:early_cutoff, :, :], axis=0) / early_cutoff

    # When looking later in the ramp, use the same number of groups
    late_cutoff = np.int(number_of_groups - early_cutoff_fraction * number_of_groups)
    late_jump_rate = np.sum(jump_map[late_cutoff:, :, :], axis=0) / (number_of_groups - late_cutoff)

    # Pixels with no CRs in the late groups have their rate set to a small
    # positive number to avoid divide by zero errors
    late_jump_rate[late_jump_rate == 0] = 1e-6

    # A potential RC or IRC pixel will have a higher jump rate early in
    # the ramp compared to later, and will have an overall high number
    # of jumps.
    jump_ratio = early_jump_rate / late_jump_rate
    potential_rc = ((jump_ratio >= jump_ratio_threshold) & (high_jumps == 1))
#    print('Number of potential_rc pixels based on Jumps: ', len(np.where(potential_rc == 1)[0]))

    return high_jumps, potential_rc, number_of_jumps


def slopes_not_cr(slope, number_of_jumps):
    """ Create an array of pixel slopes which are clean and have not detected cosmic rays

    Parameters
    ----------
    slope : numpy.ndarray
        Array of pixel slopes for integration

    number_of_jumps : numpy.ndarray
        Array holding the number of jumps detected for each pixel ramp

    Returns
    -------
    clean_slope : numpy.ndarray
        array of slopes for an integration containing no cosmic rays
         slopes for a pixel ramp that have a cosmic ray detected are set to nan
    iclean_slope: numpy.ndarray
         an array of that holds if a good slope was detected. A value of 1 is
         assigned to array element
    """

    good = number_of_jumps == 0
    bad = number_of_jumps != 0

    clean_slope = np.zeros(slope.shape, dtype=np.float)
    iclean_slope = np.zeros(slope.shape, dtype=np.int)
    clean_slope[good] = slope[good]
    clean_slope[bad] = np.nan
    iclean_slope[good] = 1
    return clean_slope, iclean_slope


def combine_clean_slopes(slope_stack, islope_stack):
    """ Combine the stack of slopes and form the stanard deviation of the pixel slopes

    Parameters
    ----------
    slope_stack : list
        A list of slopes for full array stacked for each integration
    islope_stack: list
        A list of of 1 or 0 for each integration. A 1 is a good slope and 0 is slope
        with cosmic ray
    """

    slopes = np.array(slope_stack)
    islopes = np.array(islope_stack)

    mean_slope = np.nanmean(slopes, axis=0)
    std_slope = np.nanstd(slopes, axis=0)
    num_good_array = np.sum(islopes, axis=0)
    # picked the value of 5 at random - should this be a parameter to program ?
    few_values = num_good_array < 5
    nfew_values = np.where(few_values)
    print('Number pixels with less than 5 pixel slopes to determine standard deviation',
          len(nfew_values[0]))
    std_slope[few_values] = np.nan

    return mean_slope, std_slope, num_good_array


def pedestal_stats(pedestal_array, threshold=5):
    """Get statsitics on the pedestal array corresponding to an
    integration

    Parameters
    ----------
    pedestal_array : numpy.ndarray
        Array of pedestal values

    threshold : int
        Number of sigma above or below the mean at which a pedestal
        value is considered bad

    Returns
    -------
    suspicioius_pedestal : numpy.ndarray
        Map of bad pedestal values\
    """
    median_pedestal = np.median(pedestal_array)
    stdev_pedestal = np.std(pedestal_array)

    suspicious_pedestal = ((pedestal_array > (median_pedestal + stdev_pedestal*threshold)) or
                           (pedestal_array < (median_pedestal - stdev_pedestal*threshold)))
    return suspicious_pedestal


def read_pedestal_data(filename, refpix):
    """Read in the PEDESTAL values from a *fitopt.fits file

    Parameters
    ----------
    filename : str
        Name of output file to check. This should be a *fitopt.fits file.

    refpix : tup
        4-element tuple listing the number of outer rows and columns that
        are reference pixels

    Returns
    -------
    pedestal : numpy.ndarray
        3D array of pedestal values (signal extrapolated to time=0)
    """
    with fits.open(filename) as hdulist:
        pedestal = hdulist['PEDESTAL'].data

    # Crop the reference pixels
    left, right, bottom, top = refpix

    if len(pedestal.shape) == 2:
        ydim, xdim = pedestal.shape
        pedestal = pedestal[bottom: ydim-top, left: xdim-right]
    elif len(pedestal.shape) == 3:
        nint, ydim, xdim = pedestal.shape
        pedestal = pedestal[:, bottom: ydim-top, left: xdim-right]

    return pedestal


def extract_group2(filename, refpix):
    """Read in the PEDESTAL values from a *fitopt.fits file

    Parameters
    ----------
    filename : str
        Name of uncalibrated file. This should be a *uncal.fits file.

    refpix : tup
        4-element tuple listing the number of outer rows and columns that
        are reference pixels

    Returns
    -------
    group2 : numpy.ndarray
        3D array of group 2
    """
    with fits.open(filename) as hdulist:
        dims = hdulist['SCI'].data.shape
        if len(dims) == 4:
            group2 = hdulist['SCI'].data[:, 1, :, :]
        elif len(dims) == 3:
            group2 = np.expand_dim(hdulist['SCI'].data[1, :, :], axis=0)
        nint, ydim, xdim = group2.shape
    # Crop the reference pixels
        left, right, bottom, top = refpix
        group2 = group2[:, bottom:  ydim-top, left: xdim-right]

    return group2


def read_slope_integrations(filenames):
    """Read in the science extension from a group of slope images

    Parameters
    ----------
    filenames : list
        List of fits files containing slope values

    Returns
    -------
    instrument : str
        Name of instrument associated with the data

    slope_data : numpy.ndarray
        3D array containing slope values for science pixels only.
        Reference pixels have been stripped off.

    starting_indexes : list
        List of numbers corresponding to the index numbers within slope_data
        where each file's data begins. For example, if slope_data is an array
        of size (10, 2048, 2048), and starting_indexes = [0, 5, 7, 10] then
        we can pull apart slope_data into its constituent exposures using
        slope_data[starting_indexes[0]: starting_indexes[1]], etc.

    left_cols : int
        Number of columns of reference pixels on the left side of the array

    right_cols : int
        Number of columns of reference pixels on the right side of the array

    bottom_rows : int
        Number of rows of reference pixels on the bottom of the array

    top_rows : int
        Number of rows of reference pixels on the top of the array
    """
    slope_stack = []
    starting_indexes = []
    for i, slope_file in enumerate(filenames):
        if not os.path.isfile(slope_file):
            raise FileNotFoundError('ERROR: Input slope file {} does not exist'.format(slope_file))

        with fits.open(slope_file) as hdulist:
            slope_img = hdulist['SCI'].data
            dq_int = hdulist['DQ'].data
            header = hdulist[0].header
            instrument = header['INSTRUME']
            slope_shape = slope_img.shape
            if len(slope_shape) == 2:
                dq_img = (dq_int[:, :] & dqflags.pixel['REFERENCE_PIXEL'] == 0)
            elif len(slope_shape) == 3:
                dq_img = (dq_int[0, :, :] & dqflags.pixel['REFERENCE_PIXEL'] == 0)
            else:
                raise ValueError("Slope image should be either 2D or 3D.")

        # Create a mask where 1 indicates a science pixel and 0 indicates
        # a reference pixel

        science = np.where(dq_img == 1)
        left_edge = np.min(science[1])
        right_edge = np.max(science[1]) + 1
        bottom_edge = np.min(science[0])
        top_edge = np.max(science[0]) + 1

        left_cols = left_edge
        right_cols = dq_img.shape[1] - right_edge
        bottom_rows = bottom_edge
        top_rows = dq_img.shape[0] - top_edge

        # Add to the list of starting indexes
        starting_indexes.append(len(slope_stack))

        # loop over integrations and pull out slope for int
        # Crop the reference pixels from the array.

        if len(slope_shape) == 2:
            slopes = slope_img[bottom_edge:top_edge, left_edge:right_edge]
            slope_stack.append(slopes)
        elif len(slope_shape) == 3:
            num_int = slope_shape[0]
            for i in range(num_int):
                slopes = slope_img[i, bottom_edge:top_edge, left_edge:right_edge]
                slope_stack.append(slopes)
    starting_indexes.append(len(slope_stack))
    slope_data = np.array(slope_stack)
    return instrument, slope_data, starting_indexes, (left_cols, right_cols, bottom_rows, top_rows)


def read_slope_data(filename, refpix):
    """Read in the science extension from a group of slope images

    Parameters
    ----------
    filenames : list
        List of fits files containing slope values

    refpix : tup
        4-element tuple listing the number of outer rows and columns that
        are reference pixels

    Returns
    -------
    slope : numpy.ndarray
        2D or 3D array of slope values (signal extrapolated to time=0)
        Reference pixels have been stripped off.

    """

    left, right, bottom, top = refpix
    if not os.path.isfile(filename):
        raise FileNotFoundError('ERROR: Slope file {} does not exist.'.format(filename))

    with fits.open(filename) as hdulist:
        slope_img = hdulist['SCI'].data
        # dq_int = hdulist['DQ'].data  # should we check DQ array to toss out any data
        slope_shape = slope_img.shape

        # Crop the reference pixels from the array.
        if len(slope_shape) == 2:
            ydim, xdim = slope_img.shape
            slope = slope_img[bottom: ydim-top, left: xdim-right]
            slopes = np.expand_dims(slope, axis=0)
        elif len(slope_shape) == 3:
            nint, ydim, xdim = slope_img.shape
            slopes = slope_img[:, bottom: ydim-top, left: xdim-right]
    return slopes


def read_slope_files(filenames):
    """Read in the science extension from a group of slope images

    Parameters
    ----------
    filenames : list
        List of fits files containing slope values

    Returns
    -------
    slope_data : numpy.ndarray
        3D array containing slope values for science pixels only.
        Reference pixels have been stripped off.

    left_cols : int
        Number of columns of reference pixels on the left side of the array

    right_cols : int
        Number of columns of reference pixels on the right side of the array

    bottom_rows : int
        Number of rows of reference pixels on the bottom of the array

    top_rows : int
        Number of rows of reference pixels on the top of the array
    """
    print('METADATA check turned off for testing with old NIRCAM data that is missing keywords')
    for i, filename in enumerate(filenames):
        # Read all of the slope data into an array
        slope_file = filename.replace('jump.fits', 'rateint.fits')
        with fits.open(slope_file) as hdulist:
            slope_img = hdulist['SCI'].data
            dq_img = hdulist['DQ'].data
            header = hdulist[0].header
            instrument = header['INSTRUME']

        # Create a mask where 1 indicates a science pixel and 0 indicates
        # a reference pixel

        dq_img = (dq_img & dqflags.pixel['REFERENCE_PIXEL'] == 0)

        science = np.where(dq_img == 1)
        left_edge = np.min(science[1])
        right_edge = np.max(science[1]) + 1
        bottom_edge = np.min(science[0])
        top_edge = np.max(science[0]) + 1

        left_cols = left_edge
        right_cols = dq_img.shape[1] - right_edge
        bottom_rows = bottom_edge
        top_rows = dq_img.shape[0] - top_edge

        # Crop the reference pixels from the array. Make sure we can handle
        # exposures with multiple integrations
        slope_shape = slope_img.shape
        if len(slope_shape) == 2:
            slopes = np.expand_dims(slope_img[bottom_edge:top_edge, left_edge:right_edge], axis=0)
        elif len(slope_shape) == 3:
            slopes = slope_img[:, bottom_edge:top_edge, left_edge:right_edge]
        else:
            raise ValueError("Slope image should be either 2D or 3D.")

        if i == 0:
            slope_data = copy.deepcopy(slopes)
            # if len(slope_shape) == 2:
            #    scipix = np.expand_dims(dq_img, axis=0)
            # elif len(slope_shape) == 3:
            #    scipix = copy.deepcopy(dq_img)
            # header_to_compare = copy.deepcopy(header)
        else:
            # Check to be sure the input files are consistent
            # check_metadata(header, header_to_compare)

            slope_data = np.vstack([slope_data, slopes])
            # scipix = np.vstack([scipix, np.expand_dims(dq_img, axis=0)])

    return instrument, slope_data, (left_cols, right_cols, bottom_rows, top_rows)


def saturated_in_all_groups(pedestal_array, first_group_sat):
    """Generate a list of pixels that are saturated in all groups

    Parameters
    ----------
    pedestal_array : numpy.ndarray
        3D array of pedestal values (signal extrapolated to time=0)
    first_group_sat: numpy.ndarray
        2D array of the first group DQ containing either 0 = not saturated or 2 = saturated.

    Returns
    -------
    full_saturation : tup
        Tuple of (y, x) coordinate lists (output from np.where)
    """
    full_saturation_ped0 = pedestal_array == 0
    # to be marked as saturated first_group_sat = 2 (saturated) and ped = 0
    full_saturation = np.logical_and(full_saturation_ped0, (first_group_sat == 2))
    return full_saturation.astype(int)


def plot_image(image, image_max, outdir, titleplot, fileout):
    """ Plot an Image

    Parameters
    ----------
    image : numpy.ndarray
         2D image to plot
    image_max :  float
         maximum of image to use for scaling the image
    titleplot : string
         title of the plot
    fileout : string
         output file of the plot

    Returns
    -------
      prints the plot to disk

    """
    fig = plt.figure(figsize=(9, 9))
    ax1 = fig.add_subplot(1, 1, 1)
    ysize = image.shape[0]
    xsize = image.shape[1]
    im = ax1.imshow(image, extent=[0, xsize, 0, ysize], interpolation='None',
                    cmap=cm.RdYlGn, origin='lower', vmin=0, vmax=image_max)
    plt.colorbar(im)

    ax1.set_title(titleplot)
    fig.tight_layout()
    fileout = outdir + '/' + fileout
    plt.savefig(fileout, bbox_inches='tight')
    # plt.show(block=False)
    # input('Press Enter to continue')
    plt.close()


def plot_histogram_stats(data_array, cut_limit, nbins, outdir,
                         titleplot, fileout,
                         xaxis_log=False):
    """ Plot a histogram of stats and over the upper limit cut off

    Parameters
    ----------
    data_array : numpy.ndarray
         2D data to make a histogram from
    sigma_threshold :  float
         used to plotting sigma clip line on plot
    nbins : integer
         number of bins in creating histogram
    titleplot : string
         title of the plot
    fileout : string
         output file of the plot

    Returns
    -------
      prints the plot to disk

    """

    # plot histogram
    data = data_array.flatten()
    data_good = np.isfinite(data)
    data = data[data_good]

    fig = plt.figure(figsize=(9, 9))
    ax1 = fig.add_subplot(1, 1, 1)
    h = np.histogram(data, bins=nbins)
    if not xaxis_log:
        ax1.hist(data, bins=nbins)
        ymax = np.amax(h[0])
        x = np.array([cut_limit, cut_limit])
        y = np.array([0, ymax])
        ax1.plot(x, y)
        ax1.set_xlabel(' Pixel Slope Standard Deviation')
    else:
        xh = h[1]
        data_small = np.logical_and(data > 0, data < 1)
        xsmall = np.amin(data[data_small])
        logbins = np.logspace(np.log(xsmall), np.log10(xh[-1]), len(xh))
        hlog = np.histogram(data, bins=logbins)

        ymax = np.amax(hlog[0])
        ax1.hist(data, bins=logbins)

        x = np.array([cut_limit, cut_limit])
        y = np.array([0, ymax])
        ax1.plot(x, y)
        # print('min and max histogram',np.amin(data),np.amax(data))
        ax1.set_xlim(0.01, 10)
        ax1.set_xscale('log')
        ax1.set_xlabel(' Log Pixel Slope Standard Deviation')
    # noisy flag set based on stats from clipped array

    ax1.set_ylabel(' Number of Pixels')
    num_above = len(np.where(data > cut_limit)[0])
    # print('number beyond cut',num_above)

    titleplot = titleplot + ' # beyond limit' + '{:6d}'.format(num_above)
    ax1.set_title(titleplot)
    fig.tight_layout()

    # plt.show(block=False)
    # cont = input('Press Enter to continue')
    fileout = outdir + '/' + fileout
    plt.savefig(fileout, bbox_inches='tight')
    plt.close()

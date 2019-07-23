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

from jwst_reffiles.bad_pixel_mask.bad_pixel_mask import create_dqdef


def find_bad_pix(filenames, noisy_threshold=5, max_saturated_fraction=0.5, max_jump_limit=10,
                 jump_ratio_threshold=5, early_cutoff_fraction=0.25, pedestal_sigma_threshold=5,
                 rc_fraction_threshold=0.8, low_pedestal_fraction=0.8, high_cr_fraction=0.8,
                 flag_values={'hot': ['HOT'], 'rc': ['RC'], 'low_pedestal': ['OTHER_BAD_PIXEL'], 'high_cr': ["TELEGRAPH"]},
                 do_not_use=['hot', 'rc', 'low_pedestal', 'high_cr'], outfile=None):
    """MAIN FUNCTION

    Parameters
    ----------
    filenames : list
        List of dark current files. These should be slope images.

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

    # Read in the slope data. Strip off reference pixels.
    # Return a 3D array of slopes and a 3D array mapping where the
    # science pixels are.
    print('Reading slope files...')
    slopes, refpix_additions = read_slope_files(filenames)
    print('Searching for noisy pixels')

    # Calculate the mean and standard deviation through the stack for
    # each pixel. Assuming that we are looking for noisy pixels, we don't
    # want to do any sigma clipping on the inputs here, right?
    mean_slope = np.mean(slopes, axis=0)
    std_slope = np.std(slopes, axis=0)

    # Use sigma-cliping when calculating the mean and standard deviation
    # of the standard deviations
    clipped_stdevs, cliplow, cliphigh = sigmaclip(std_slope, low=clipping_sigma, high=clipping_sigma)
    avg_of_std = np.mean(clipped_stdevs)
    std_of_std = np.std(clipped_stdevs)

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
    counter = 0
    print('Working on jump step outputs.')
    for i, filename in enumerate(filenames):
        print('Opening {}'.format(filename))
        # Read in the ramp and get the data and dq arrays
        ramp_file = filename.replace('_ramp_fit_0.fits', '_jump.fits')
        groupdq = get_jump_dq_values(ramp_file, refpix_additions)

        # Generate a map of JUMP flags in the ramp
        cr_map = get_cr_flags(groupdq)

        # Read in the fitops file associated with the exposure and get
        # the pedestal array
        pedestal_file = filename.replace('_ramp_fit_0.fits', '_fitopt.fits')
        pedestal = read_pedestal_data(pedestal_file, refpix_additions)

        # Work one integration at a time
        for int_num in range(pedestal.shape[0]):
            pedestal_int = pedestal[int_num, :, :]
            clipped_pedestal, cliplow, cliphigh = sigmaclip(pedestal_int, low=3., high=3.)
            mean_pedestal = np.mean(clipped_pedestal)
            std_pedestal = np.std(clipped_pedestal)

            rc_from_pedestal[counter, :, :] = pedestal_int > (mean_pedestal + std_pedestal * pedestal_sigma_threshold)

            # Pixels with abnormally low pedestal values
            low_pedestal[counter, :, :] = pedestal_int < (mean_pedestal - std_pedestal * pedestal_sigma_threshold)

            # Find pixels that are saturated in all groups. These will have
            # a pedestal value of 0 (according to the pipeline documentation).
            # These should end up flagged as HOT and DO_NOT_USE
            saturated[counter, :, :] = saturated_in_all_groups(pedestal_int)

            # Find pixels that have an abnormally high number of jumps, as
            # well as those that have most of their jumps concentrated in the
            # early part of the integration. The latter are possibly RC or IRC
            # pixels
            many_jumps, rc_candidates = find_pix_with_many_jumps(cr_map[int_num, :, :, :], max_jump_limit=10,
                                                                 jump_ratio_threshold=5,
                                                                 early_cutoff_fraction=0.25)
            high_cr_rate[counter, :, :] = copy.deepcopy(many_jumps)
            rc_from_flags[counter, :, :] = copy.deepcopy(rc_candidates)
            counter += 1

    # Look through the stack of saturated pixels and keep those saturated
    # more than N% of the time
    fully_saturated = np.sum(saturated, axis=0) / len(filenames)
    fully_saturated[fully_saturated < max_saturated_fraction] = 0
    fully_saturated = np.ceil(fully_saturated).astype(np.int)
    fully_saturated = apply_flags(fully_saturated, flag_values['hot'])
    num_saturated = len(np.where(fully_saturated != 0)[0])
    print('\n\nFound {} fully saturated pixels.'.format(num_saturated))

    # How do we want to combine these to identify RC pixels?
    rc_pedestal = np.sum(rc_from_pedestal, axis=0) / len(filenames)
    rc_flags = np.sum(rc_from_flags, axis=0) / len(filenames)

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
    low_pedestal_vals = np.sum(low_pedestal, axis=0) / len(filenames)
    low_ped = low_pedestal_vals > low_pedestal_fraction
    low_ped = apply_flags(low_ped.astype(np.int), flag_values['low_pedestal'])
    num_low_ped = len(np.where(low_ped != 0)[0])
    print('Found {} low pedestal pixels.'.format(num_low_ped))

    # Pixels with lots of CR flags should be added to the list of noisy pixels?
    high_cr = np.sum(high_cr_rate, axis=0) / len(filenames)
    noisy_second_pass = high_cr > high_cr_fraction
    combined_noisy = np.bitwise_or(noisy, noisy_second_pass)
    combined_noisy = apply_flags(combined_noisy.astype(np.int), flag_values['high_cr'])
    num_noisy = len(np.where(noisy != 0)[0])
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
    if outfile is None:
        outfile = 'badpixels_from_darks.fits'
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
        Map of jump flags for all pixels (e.g. output from ``get_cr_flags``)
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
    print('Number of potential_rc pixels based on Jumps: ', len(np.where(potential_rc == 1)[0]))

    return high_jumps, potential_rc


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
    cr_map = (dq_array & dqflags.pixel['JUMP_DET'] > 0)
    return cr_map


def get_jump_dq_values(filename, refpix):
    """Read in the GROUPDQ extension of a file and crop reference pixels

    Parameters
    ----------
    filename : str
        fits file to be read in

    refpix : tup
        4-element tuple giving the left, right, bottom, and top number
        of rows and columns of reference pixels

    Returns
    -------
    groupdq : numpy.ndarray
        4D array of DQ values
    """
    with fits.open(filename) as hdulist:
        groupdq = hdulist['GROUPDQ'].data

    if len(groupdq.shape) != 4:
        raise ValueError("ERROR: expecting groupdq to be a 4D array")

    nint, ngroup, ydim, xdim = groupdq.shape
    left, right, bottom, top = refpix
    return groupdq[:, :, bottom: ydim-top, left: xdim-right]


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
        slope_file = filename.replace('jump.fits', 'rate.fits')
        with fits.open(slope_file) as hdulist:
            slope_img = hdulist['SCI'].data
            dq_img = hdulist['DQ'].data
            header = hdulist[0].header
        # Create a mask where 1 indicates a science pixel and 0 indicates
        # a reference pixel
        dq_img = (dq_img & dqflags.pixel['REFERENCE_PIXEL'] == 0)

        science = np.where(dq_img == 1)
        left_edge = np.min(science[1])
        right_edge = np.max(science[1]) + 1
        bottom_edge = np.min(science[0])
        top_edge = np.max(science[1]) + 1

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
            #if len(slope_shape) == 2:
            #    scipix = np.expand_dims(dq_img, axis=0)
            #elif len(slope_shape) == 3:
            #    scipix = copy.deepcopy(dq_img)
            header_to_compare = copy.deepcopy(header)
        else:
            # Check to be sure the input files are consistent
            #check_metadata(header, header_to_compare)

            slope_data = np.vstack([slope_data, slopes])
            #scipix = np.vstack([scipix, np.expand_dims(dq_img, axis=0)])
    return slope_data, (left_cols, right_cols, bottom_rows, top_rows)


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

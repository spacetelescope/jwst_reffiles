#! /usr/bin/env python

"""This module creates a superbias reference file that can be used
in the JWST calibration pipeline.

Author
------
     - Ben Sunnquist

Use
---
    This module can be imported and used as such:
    ::
        from jwst_reffiles.superbias import superbias
        superbias.make_superbias(arguments)

Notes
-----
    Overview:
    Inputs: A list of dark current ramps

    Algorithm: The superbias is the sigma-clipped mean of all 0th frames.
               The superbias error is the sigma-clipped stddev of all 0th frames.
               Pixels with high superbias error values are flagged as
               UNRELIABLE_BIAS and DO_NOT_USE in the superbias DQ array.
"""

import os
from copy import deepcopy

from astropy.io import fits
from astropy.stats import sigma_clip
import datetime as dt
import json
from jwst.datamodels import dqflags, SuperBiasModel
from sklearn.mixture import GaussianMixture
from stdatamodels import util
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from scipy.ndimage import median_filter


def bias_trend_plot(array, filename_integ_list, mjd_list, aperture, output_file):
    """Calculate the median bias for each integration in the stack, and plot

    Parameters
    ----------
    array : numpy.ndarray
        3D array of bias frames

    filename_integ_list : list
        List of strings such as "jw01068001001_02101_00001_nrcalong/3" that provide the filename and
        integration number associated with each value in ''array''

    mjd_list : list
        List of MJD associated with each value in ''array''

    aperture : str
        Aperture name, for the plot title

    output_file : str
        Name of the superbias output file. The trending plot will be saved under a similar filename

    Returns
    -------
    trending_file : str
        Name of file to which the trending plot has been saved
    """
    nframe, ny, nx = array.shape

    hw = 10  # pixels
    ystart = ny // 2 - hw
    yend = ystart + (2 * hw)
    xstart = nx // 2 - hw
    xend = xstart + (2 * hw)
    bias = []
    for frame in range(nframe):
        bias.append(np.nanmedian(array[frame, ystart:yend, xstart:xend]))

    # Convert MJD values to datetime values that will be more intuitive on the plot
    dates = [dt.datetime(1970, 1, 1) + dt.timedelta(days=m - 40587)
         for m in mjd_list]

    detector = filename_integ_list[0].split('_')[-1].upper().split('/')[0]
    trending_file = f"{output_file.split('.')[0]}_trending.png"

    f, a = plt.subplots(figsize=(10, 8))
    a.scatter(dates, bias, marker='o')
    a.set_xlabel('Date')
    a.set_ylabel('Median bias level (DN)')
    a.set_title(f'{aperture} Bias Trending')
    a.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %Y'))
    f.autofmt_xdate()  # rotates labels nicely
    f.savefig(trending_file)
    print(f'Save superbias trending plot to {trending_file}')
    print('To do: maybe turn this into a bokeh figure where we can hover over points and get filename/int info')


def calculate_mean(stack, clipping_sigma=3.0, max_clipping_iters=3, mean_stat='mean'):
    """Calculates the sigma-clipped mean through a stack of images.

    Parameters
    ----------
    stack : numpy.ndarray
        A 3D stack of images.

    clipping_sigma : float
        Number of sigmas to use when sigma-clipping the input stack.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping the input
        stack.

    mean_stat : string
        Can be 'mean' or 'median'. Controls which statistic is used on the clipped data

    Returns
    -------
    mean_image : numpy.ndarray
        2D image of the sigma-clipped mean through the input stack.
    """

    clipped = sigma_clip(stack, sigma=clipping_sigma,
                         maxiters=max_clipping_iters, axis=0)
    if mean_stat == 'mean':
        mean_image = np.nanmean(clipped, axis=0)
    elif mean_stat == 'median':
        mean_image = np.nanmedian(clipped, axis=0)

    return mean_image


def calculate_stddev(stack, clipping_sigma=3.0, max_clipping_iters=3):
    """Calculates the sigma-clipped standard deviation through a stack
    of images.

    Parameters
    ----------
    stack : numpy.ndarray
        A 3D stack of images.

    clipping_sigma : float
        Number of sigmas to use when sigma-clipping the input stack.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping the input
        stack.

    Returns
    -------
    stddev_image : numpy.ndarray
        2D image of the sigma-clipped standard deviation through the
        input stack.
    """

    clipped = sigma_clip(stack, sigma=clipping_sigma,
                         maxiters=max_clipping_iters, axis=0)
    stddev_image = np.nanstd(clipped, axis=0)

    return stddev_image


def find_low_bias(array):
    """Use a Gaussian Mixture Model to locate the integrations within the exposure that exhibit low
    bias signal
    """
    array = array.reshape(-1, 1)
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(array)

    labels = gmm.predict(array)
    means = gmm.means_.flatten()

    bad_cluster = np.argmin(means)
    bad_indices = np.where(labels == bad_cluster)[0]

    #probs = gmm.predict_proba(array)
    #bad_prob = probs[:, bad_cluster]

    # Sanity check: the number of low bias entries shouldn't be more than 1 of every 4
    # integrations. If more were found, then we consider this test a failure and include
    # all integrations
    if len(bad_indices) > 0.25 * len(array):
        print('Too many low bias integrations found. Keeping them all.')
        bad_indices = np.array([], dtype=np.int64)

    print(bad_indices)
    return bad_indices


def find_unreliable_bias_pixels(error, bad_clipping_sigma=3.0,
                                bad_max_clipping_iters=3,
                                bad_sigma_threshold=5.0, plot=False):
    """Identifies pixels with high superbias error values (i.e.
    UNRELIABLE_BIAS pixels).

    Parameters
    ----------
    error : numpy.ndarray
        A 2D superbias error image.

    bad_clipping_sigma : float
        Number of sigma to use when sigma-clipping the superbias error
        image to find UNRELIABLE_BIAS pixels.

    bad_max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping the superbias
        error image to find UNRELIABLE_BIAS pixels.

    bad_sigma_threshold : float
        Number of standard deviations above the mean of the superbias error
        image at which a pixel is flagged as UNRELIABLE_BIAS.

    plot : bool
        Option to save a histogram of the superbias error image values, with
        the threshold used for flagging UNRELIABLE_BIAS pixels also shown.

    Returns
    -------
    dq : numpy.ndarray
        A 2D DQ image that flags pixels with high superbias error values as
        UNRELIABLE_BIAS (2) and DO_NOT_USE (1).
    """

    # Make initial DQ array and flag pixels with nan/inf superbias error values
    dq = np.zeros(error.shape, dtype=int)
    dq[~np.isfinite(error)] = dqflags.pixel['UNRELIABLE_BIAS']

    # Flag pixels X sigma greater than the mean in the error image
    clipped = sigma_clip(error, sigma=bad_clipping_sigma,
                         maxiters=bad_max_clipping_iters)
    mean = np.mean(clipped)
    std = np.std(clipped)
    thresh = mean + bad_sigma_threshold*std
    dq[error > thresh] = 3

    # Plot a histogram of the error image values
    if plot:
        error_no_nan = error.flatten()
        error_no_nan = error_no_nan[np.isfinite(error_no_nan)]
        plt.hist(error_no_nan, range=(0, thresh*3), bins=80,
                 color='steelblue')
        plt.axvline(thresh, ls='--', color='red')
        plt.yscale('log')
        plt.ylabel('Number of Pixels')
        plt.xlabel('Superbias Error [DN]')
        plt.title('# of UNRELIABLE_BIAS pixels: {}'.format(len(dq[dq!=0])))
        plt.savefig('superbias_error_hist.png', dpi=100, bbox_inches='tight')
        plt.clf()

    return dq


def get_crds_info(filename, subarray, readpatt):
    """Get CRDS-required header information to put into the final superbias
    reference file.

    Parameters
    ----------
    filename : str
        Path to one of the files being used to generate the superbias
        reference file. Will be used to find default values.

    subarray : str
        CRDS-required subarray for which to use this reference file for.

    readpatt : str
        CRDS-required read pattern for which to use this reference file for.

    Returns
    -------
    instrument : str
        CRDS-required instrument for which to use this reference file for.

    detector : str
        CRDS-required detector for which to use this reference file for.

    subarray : str
        CRDS-required subarray for which to use this reference file for.

    readpatt : str
        CRDS-required read pattern for which to use this reference file for.

    fastaxis : int
        CRDS-required fastaxis of the reference file.

    slowaxis : int
        CRDS-required slowaxis of the reference file.

    substrt1 : int
        CRDS-required starting pixel in axis 1 direction.

    substrt2 : int
        CRDS-required starting pixel in axis 2 direction.
    """

    header = fits.getheader(filename)
    instrument = header['INSTRUME']
    detector = header['DETECTOR']
    fastaxis = header['FASTAXIS']
    slowaxis = header['SLOWAXIS']
    substrt1 = header['SUBSTRT1']
    substrt2 = header['SUBSTRT2']

    if subarray is None:
        subarray = header['SUBARRAY']
    if readpatt is None:
        readpatt = 'N/A'

    return instrument, detector, subarray, readpatt, fastaxis, slowaxis, substrt1, substrt2


def strip_filename_suffixes(filename):
    """Strip away arbitrary suffixes after the detector in a filename
    """
    detectors = ['nrcalong', 'nrcblong', 'nrca1', 'nrca2', 'nrca3', 'nrca4', 'nrcb1', 'nrcb2', 'nrcb3', 'nrcb4']
    pattern = re.compile(rf'^(.*?_({"|".join(detectors)}))')

    match = pattern.match(filename)
    if not match:
        raise ValueError(f"No detector found in filename: {filename}")
    return match.group(1)


def get_dark_data(files):
    """Read 0th groups from dark current data into an array
    """
    # Create dictionary to hold the list of good integrations for each file
    good_idx = {}

    try:
        stack = fits.getdata(files[0], 'ZEROFRAME').astype(float)
    except:
        stack = fits.getdata(files[0], 'SCI')[:, 0, :, :].astype(float)

    # Create a list of filename/integ_number for all integrations in the file
    hdr = fits.getheader(files[0])
    nints = hdr['NINTS']
    filenames = [get_file_base(hdr['FILENAME'])] * nints
    intnums = np.arange(nints)
    fileint = [f'{filename}/{intnum}' for filename, intnum in zip(filenames, intnums)]

    # Get the observation date time
    obsdate = hdr['EXPSTART']
    tint = hdr['EFFINTTM'] / 3600. / 24.
    obstimes = [obsdate + tint * inum for inum in range(nints)]

    # Strip away the integrations with the low bias due to cycling
    stack, good = remove_bias_cycling(stack)
    fileint = fileint[good]
    obstimes = obstimes[good]

    # Use the filename base as a key. Remove extra suffixes in order to make it easier
    # for other reffile generation scripts to use
    file_base = get_file_base(files[0])
    good_idx[file_base] = good

    for filename in files[1:]:
        try:
            data = fits.getdata(filename, 'ZEROFRAME')[:, 0, :, :].astype(float)
        except:
            data = fits.getdata(filename, 'SCI')[:, 0, :, :].astype(float)

        # Create a list of filename/integ_number for all integrations in the file
        hdr = fits.getheader(filename)
        nints = hdr['NINTS']
        filenames = [get_file_base(hdr['FILENAME'])] * nints
        intnums = np.arange(nints)
        fileintstr = [f'{filename}/{intnum}' for filename, intnum in zip(filenames, intnums)]

        # Get the observation date time
        obsdate = hdr['EXPSTART']
        tint = hdr['EFFINTTM'] / 3600. / 24.
        exp_obstimes = [obsdate + tint * inum for inum in range(nints)]

        # Strip away the integrations with the low bias due to cycling
        data, good = remove_bias_cycling(data)
        fileintstr = np.array(fileintstr)[good]
        fileint.extend(fileintstr)
        exp_obstimes = np.array(exp_obstimes)[good]
        obstimes.extend(exp_obstimes)

        file_base = get_file_base(filename)
        good_idx[file_base] = good
        stack = np.concatenate([stack, data])
    return stack, good_idx, fileint, obstimes


def get_and_filter_dark_data(files, good_ints):
    """Read in 0th frames of all integrations for the given file list. Then using good_ints,
    which lists the integration numbers that we want to keep, filter out the "bad" integrations.
    The motivation here is to filter out all of the integrations where the bias level is different
    because the fast reset intersected the subarray

    Paramters
    ---------
    files : list
        List of filenames

    good_ints : dict
        Dictionary with keys that are file bases (i.e. no suffix), and the values are lists of
        good integration numbers

    Returns
    -------
    stack : numpy.ndarray
        3D stack of 0th frames
    """
    for i, filename in enumerate(files):
        try:
            data = fits.getdata(filename, 'ZEROFRAME').astype(float)
        except:
            data = fits.getdata(filename, 'SCI')[:, 0, :, :].astype(float)

        file_base = get_file_base(filename)

        if file_base in good_ints:
            idx = np.array(good_ints[file_base])
            data = data[idx, :, :]

        if i == 0:
            stack = deepcopy(data)
        else:
            stack = np.concatenate([stack, data])
    return stack


def get_file_base(filename):
    """Get the filename base. i.e. remove suffixes
    """
    filename = os.path.basename(filename)
    nrc_index = filename.find('nrc')
    if nrc_index != -1:
        underscore_idx = filename.find('_', nrc_index)
        if underscore_idx != -1:
            file_base = filename[0:underscore_idx]
            return file_base
        else:
            return filename
    else:
        return filename


def get_pedestal_data(files):
    """Read pedestal data from fitops file into an array
    """
    good_idx = {}
    fileint = []
    obstimes = []
    for i, filename in enumerate(files):
        print(f'{filename}:')
        array = fits.getdata(filename, 'PEDESTAL')

        # Create a list of filename/integ_number for all integrations in the file
        hdr = fits.getheader(filename)
        nints = hdr['NINTS']
        filenames = [get_file_base(hdr['FILENAME'])] * nints
        intnums = np.arange(nints)
        exp_fileint = [f'{filename}/{intnum}' for filename, intnum in zip(filenames, intnums)]

        # Get the observation date time
        obsdate = hdr['EXPSTART']
        tint = hdr['EFFINTTM'] / 3600. / 24.
        exp_obstimes = [obsdate + tint * inum for inum in range(nints)]

        if i == 0:
            # Identify and remove integrations where the bias cycling has resulted in a low bias.
            # See https://jira.stsci.edu/browse/JP-2601 for details
            stack, good = remove_bias_cycling(array)
        else:
            arr, good = remove_bias_cycling(array)
            stack = np.concatenate([stack, arr])

        # Append to the list of filenames/ints and obstimes
        exp_fileint = np.array(exp_fileint)[good]
        exp_obstimes = np.array(exp_obstimes)[good]
        fileint.extend(exp_fileint)
        obstimes.extend(exp_obstimes)

        file_base = get_file_base(filename)
        good_idx[file_base] = good
    return stack, good_idx, fileint, obstimes

def get_all_good_pedestal_info(files):
    """Create the outputs from get_pedestal_data, excepting the actual pixel data
    """
    good_idx = {}
    fileint = []
    obstimes = []
    for i, filename in enumerate(files):
        print(f'{filename}:')
        #array = fits.getdata(filename, 'PEDESTAL')

        # Create a list of filename/integ_number for all integrations in the file
        hdr = fits.getheader(filename)
        nints = hdr['NINTS']
        filenames = [get_file_base(hdr['FILENAME'])] * nints
        intnums = np.arange(nints)
        exp_fileint = [f'{filename}/{intnum}' for filename, intnum in zip(filenames, intnums)]

        # Get the observation date time
        obsdate = hdr['EXPSTART']
        tint = hdr['EFFINTTM'] / 3600. / 24.
        exp_obstimes = [obsdate + tint * inum for inum in range(nints)]

        #if i == 0:
        #    # Identify and remove integrations where the bias cycling has resulted in a low bias.
        #    # See https://jira.stsci.edu/browse/JP-2601 for details
        #    stack, good = remove_bias_cycling(array)
        #else:
        #    arr, good = remove_bias_cycling(array)
        #    stack = np.concatenate([stack, arr])

        # Append to the list of filenames/ints and obstimes
        exp_fileint = np.array(exp_fileint)
        exp_obstimes = np.array(exp_obstimes)
        fileint.extend(exp_fileint)
        obstimes.extend(exp_obstimes)

        file_base = get_file_base(filename)
        good_idx[file_base] = np.arange(nints)
    return good_idx, fileint, obstimes


def make_filled_superbias(superbias, refpix_map, bad_clipping_sigma=3.0,
                          bad_max_clipping_iters=3, bad_sigma_threshold=5.0,
                          box_width=3, outfile='superbias_jwst_reffiles.fits'):
    """Create a version of the superbias image where the values of high
    superbias pixels are replaced by their smoothed values.

    Parameters
    ----------
    superbias : numpy.ndarray
        A 2D superbias image.

    refpix_map : numpy.ndarray
        A 2D image where refpix have a value of 1 and all other pixels have
        a value of 0.

    bad_clipping_sigma : float
        Number of sigma to use when sigma-clipping to find UNRELIABLE_BIAS
        pixels.

    bad_max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping to find
        UNRELIABLE_BIAS pixels.

    bad_sigma_threshold : float
        Number of standard deviations above the mean at which a pixel is
        flagged as UNRELIABLE_BIAS.

    box_width : int
        Width in pixels of the box kernel to use to compute the smoothed
        superbias image.

    outfile : str
        Name of the CRDS-formatted superbias reference file to save the final
        superbias map to. This is only used in this function to name the
        filled output file.
    """

    # Trim the refpix from the input superbias
    scipix = np.where(refpix_map==0)
    ymin = np.min(scipix[0])
    xmin = np.min(scipix[1])
    ymax = np.max(scipix[0]) + 1
    xmax = np.max(scipix[1]) + 1
    superbias_no_refpix = superbias[ymin:ymax, xmin:xmax].copy()

    # Make initial DQ arrays, with and without refpix, to store pixels
    # with high superbias values.
    dq = np.zeros(superbias.shape, dtype=int)
    dq_no_refpix = np.zeros(superbias_no_refpix.shape, dtype=int)

    # Median-filter the superbias image; replace any nans/infs first as
    # these mess up the smoothing. Then subtract this smoothed image from
    # the original superbias image.
    superbias_no_refpix[~np.isfinite(superbias_no_refpix)] = \
        np.nanmedian(superbias_no_refpix)
    smoothed = median_filter(superbias_no_refpix, box_width)
    diff = superbias_no_refpix - smoothed

    # Flag pixels X sigma greater than the mean in the difference image
    clipped = sigma_clip(diff, sigma=bad_clipping_sigma,
                         maxiters=bad_max_clipping_iters)
    mean = np.mean(clipped)
    std = np.std(clipped)
    thresh = mean + bad_sigma_threshold*std
    dq_no_refpix[diff > thresh] = 3

    # Put DQ array back into original dimensions
    dq[ymin:ymax, xmin:xmax] = dq_no_refpix

    # Save an image where high superbias pixels are replaced by their smoothed
    # values.
    outname = outfile.replace('.fits', '_filled.fits')
    superbias_filled_no_refpix = superbias[ymin:ymax, xmin:xmax].copy()
    superbias_filled_no_refpix[dq_no_refpix != 0] = \
        smoothed[dq_no_refpix != 0]
    superbias_filled = superbias.copy()
    superbias_filled[ymin:ymax, xmin:xmax] = superbias_filled_no_refpix
    fits.writeto(outname, superbias_filled, overwrite=True)


def make_refpix_map(filename):
    """Creates a map of the reference pixels based on the input fits file.

    Parameters
    ----------
    filename : str
        The fits image file to create the refpix map from. The PIXELDQ
        extension is used to find the refpix.

    Returns
    -------
    refpix_map : numpy.ndarray
        A map of the reference pixels (same dimensions as the input fits file
        image extensions), where refpix are 1 and all other pixels are 0.
        An array of all zeros is returned if no refpix are found or if no
        PIXELDQ extension exists.
    """

    try:
        dq = fits.getdata(filename, 'PIXELDQ')
        refpix_map = (dq & dqflags.pixel['REFERENCE_PIXEL'] > 0).astype('int')
        n_refpix = len(refpix_map[refpix_map != 0])
        if n_refpix == 0:
            print('No reference pixels found in PIXELDQ extension; all pixels '
                  'will be used when smoothing the superbias to create the '
                  'filled superbias image.')
    except KeyError:
        print('Warning: No PIXELDQ extension found in input file; all pixels '
              'will be used when smoothing the superbias to create the filled '
              'superbias image.')
        data = fits.getdata(filename, 'SCI')[0, 0, :, :]
        refpix_map = np.zeros(data.shape).astype('int')

    return refpix_map


def make_superbias(ramp_filenames=None, fitopt_filenames=None, extname='SCI', clipping_sigma=3.0, max_clipping_iters=3,
                   bad_clipping_sigma=3.0, bad_max_clipping_iters=3,
                   bad_sigma_threshold=5.0, plot=False, save_filled=False,
                   box_width=3, save_reffile=True, mean_stat='mean',
                   outfile='superbias_jwst_reffiles.fits',
                   author='jwst_reffiles', description='Super Bias Image',
                   pedigree='GROUND', useafter='2000-01-01T00:00:00',
                   history='', subarray=None, readpatt=None):
    """The main function. Creates a superbias reference file using the input
    dark current ramps. See module docstring for more details.

    Parameters
    ----------
    filenames : list
        List of dark current ramp files. The data shape in these images is
        assumed to be a 4D array in DMS format (integration, group, y, x).

    extname : str
        Extension name from which the input data are retrieved. When working with
        dark current images, use 'SCI'. With fitopt files, use 'PEDESTAL'

    clipping_sigma : float
        Number of sigma to use when sigma-clipping through the 0th frame
        stack to create the superbias image.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping through
        the 0th frame stack to create the superbias error image.

    bad_clipping_sigma : float
        Number of sigma to use when sigma-clipping the superbias error image
        to find UNRELIABLE_BIAS pixels.

    bad_max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping the superbias
        error image to find UNRELIABLE_BIAS pixels.

    bad_sigma_threshold : float
        Number of standard deviations above the mean of the superbias error
        image at which a pixel is flagged as UNRELIABLE_BIAS.

    plot : bool
        Option to save a histogram of the superbias error image values, with
        the threshold used for flagging UNRELIABLE_BIAS pixels also shown.

    save_filled : bool
        Option to save an image where high superbias pixels are replaced by
        their smoothed values.

    box_width : int
        Width in pixels of the box kernel to use to compute the smoothed
        superbias image. Only relevant if save_filled=True.

    save_reffile : bool
        Option to save the generated superbias map into a CRDS-formatted
        reference file.

    mean_stat : str
        Can be 'mean' or 'median'. Controls which statistic is used through the
        stack of input data.

    outfile : str
        Name of the CRDS-formatted superbias reference file to save the final
        superbias map to.

    author : str
        CRDS-required name of the reference file author, to be placed in the
        referece file header.

    description : str
        CRDS-required description of the reference file, to be placed in the
        reference file header.

    pedigree : str
        CRDS-required pedigree of the data used to create the reference file.

    useafter : str
        CRDS-required date of earliest data with which this referece file
        should be used (e.g. '2019-04-01T00:00:00').

    history : str
        CRDS-required history section to place in the reference file header.

    subarray : str
        CRDS-required subarray for which to use this reference file for.

    readpatt : str
        CRDS-required read pattern for which to use this reference file for.
    """

    # Make a stack of all 0th frames
    #if extname == 'SCI':
    #    stack, good_ints, file_integ_list, integ_mjd_list = get_dark_data(filenames)
    #elif extname == 'PEDESTAL':
    #    stack, good_ints, file_integ_list, integ_mjd_list = get_pedestal_data(filenames)

    # Turns out it's much easier to locate the high/low bias integrations in the pedestal
    # data than it is when using 0th reads. So, when working with subarrays, let's try
    # reading in the pedestal data, finding the high/low bias, and then, in cases where
    # the 0th read is to be used, read in the 0th group data, and apply the masks from
    # the pedestal data
    if ramp_filenames is not None:
        check_filename = ramp_filenames[0]
        filenames = ramp_filenames
    else:
        check_filename = fitopt_filenames[0]
        filenames = fitopt_filenames

    subarray_name = fits.getheader(check_filename)['SUBARRAY']
    aperture = fits.getheader(check_filename)['APERNAME']
    if (extname == 'PEDESTAL' or 'FULL' not in subarray_name):
        stack, good_ints, file_integ_list, integ_mjd_list = get_pedestal_data(fitopt_filenames)
    elif 'FULL' in subarray_name:
        # In this case, we just build up the lists of good integrations, etc to include all ints
        good_ints, file_integ_list, integ_mjd_list = get_all_good_pedestal_info(ramp_filenames)
    if extname == 'SCI':
        # If we want to use 0th frames, we overwrite what's in stack with that data
        stack = get_and_filter_dark_data(ramp_filenames, good_ints)

    # Save the information about which integrations are good, so that it can be used
    # in the construction of other reference files
    outdir = os.path.dirname(outfile)
    if 'FULL' not in aperture:
        integ_filename = os.path.join(outdir, f'{aperture}_good_bias_integrations.json')

        with open(integ_filename, 'w') as f:
            json.dump(good_ints, f, indent=4)
        print(f"Dictionary of filename and good integrations saved to {integ_filename}")
    else:
        print(f'Full frame data. No list of high/low bias integrations saved to file.')
        integ_filename = None

    # Create a plot of the bias level within each of the good input integrations
    trending_plot = bias_trend_plot(stack, file_integ_list, integ_mjd_list, aperture, os.path.join(outdir, outfile))

    num_ints = stack.shape[0]

    # Calculate the superbias as the sigma-clipped mean through the stack
    superbias = calculate_mean(stack, clipping_sigma=clipping_sigma,
                               max_clipping_iters=max_clipping_iters,
                               mean_stat=mean_stat)

    # Convert masked array to normal numpy array and check for any missing
    # data.
    #superbias = superbias.fill(fill_value=np.nan)
    n_bad = len(superbias[~np.isfinite(superbias)])
    print(f'Found {n_bad} NaN pixels')
    if n_bad > 0:
        print('Warning: Superbias file has {} nan/inf pixels.'.format(n_bad))

    if save_reffile:

        # Inputs listed as None in the config file are read in as strings.
        # Change these to NoneType objects.
        subarray = none_check(subarray)
        readpatt = none_check(readpatt)

        # CRDS doesnt allow descriptions over 65 characters
        #if len(description) > 65:
        #    raise ValueError('Description cannot exceed 65 characters.')

        # Get info needed by CRDS to put into the final superbias reference file
        instrument, detector, subarray, readpatt, fastaxis, slowaxis, substrt1, substrt2 = \
            get_crds_info(check_filename, subarray, readpatt)

        # Calculate the superbias error as the sigma-clipped stddev through the stack
        error = calculate_stddev(stack, clipping_sigma=clipping_sigma,
                                 max_clipping_iters=max_clipping_iters)

        # Convert masked array to normal numpy array and check for any missing data
        error = error.filled(fill_value=np.nan)
        n_bad = len(error[~np.isfinite(error)])
        if n_bad > 0:
            print('Warning: Superbias Error array has {} nan/inf pixels.'
                  .format(n_bad))

        # Flag pixels with high superbias error values as UNRELIABLE_BIAS
        dq = find_unreliable_bias_pixels(error,
                                         bad_clipping_sigma=bad_clipping_sigma,
                                         bad_max_clipping_iters=bad_max_clipping_iters,
                                         bad_sigma_threshold=bad_sigma_threshold,
                                         plot=plot)

        # Save a version of the superbias where high superbias values are
        # replaced by their smoothed values.
        if save_filled:
            refpix_map = make_refpix_map(check_filename)
            make_filled_superbias(superbias, refpix_map,
                                  bad_clipping_sigma=bad_clipping_sigma,
                                  bad_max_clipping_iters=bad_max_clipping_iters,
                                  bad_sigma_threshold=bad_sigma_threshold,
                                  box_width=box_width, outfile=outfile)

        # Save CRDS-formatted superbias reference file
        save_superbias(superbias, error, dq, instrument=instrument,
                       detector=detector, subarray=subarray, readpatt=readpatt,
                       outfile=outfile, author=author, description=description,
                       pedigree=pedigree, useafter=useafter, history=history,
                       fastaxis=fastaxis, slowaxis=slowaxis, substrt1=substrt1,
                       substrt2=substrt2, filenames=filenames, num_ints=num_ints)

    return superbias, integ_filename, file_integ_list, integ_mjd_list


def none_check(value):
    """If value is a string containing 'none', then change it to a
    NoneType object.

    Parameters
    ----------
    value : str or NoneType

    Returns
    -------
    new_value : NoneType
    """

    if isinstance(value, str):
        if 'none' in value.lower():
            new_value = None
        else:
            new_value = value
    else:
        new_value = value

    return new_value


def remove_bias_cycling(array):
    """Remove the integrations where the detector reset passes through the subarray. This leads to a bias
    signal that is much lower than the case where the reset does not pass through the subarray. Details are
    given in https://jira.stsci.edu/browse/JP-2601
    """
    nints, ny, nx = array.shape
    if ((ny < 2048) & (nx < 2048)):
        hw = 10
        r1 = array[:, ny//2-hw:ny//2+hw, 5:5+2*hw]
        r1stats = np.nanmedian(r1, axis=(1,2))
        r1_low_idx = find_low_bias(r1stats)

        r2 = array[:, ny-4-2*hw:ny-4, nx//2-hw:nx//2+hw]
        r2stats = np.nanmedian(r2, axis=(1,2))
        r2_low_idx = find_low_bias(r2stats)

        r3 = array[:, ny//2-hw:ny//2+hw, nx-4-2*hw:nx-4]
        r3stats = np.nanmedian(r3, axis=(1,2))
        r3_low_idx = find_low_bias(r3stats)

        r4 = array[:, 5:5+2*hw, nx//2-hw:nx//2+hw]
        r4stats = np.nanmedian(r4, axis=(1,2))
        r4_low_idx = find_low_bias(r4stats)

        allidx = np.concatenate([r1_low_idx, r2_low_idx, r3_low_idx, r4_low_idx])
        unique_values, counts = np.unique(allidx, return_counts=True)

        # In order to be considered a bad bias value, at least 2 of the values
        # must be low
        bad_idx = []
        for val, ct in zip(unique_values, counts):
            if ct >= 2:
                bad_idx.append(val)
        bad_idx = np.array(bad_idx)

        print('Integrations with low bias: ', bad_idx)
        print(r1stats)
        print(r2stats)

        # Create a boolean mask initialized to True (all are good initially)
        boolean_mask = np.ones(nints, dtype=bool)

        # Set the "bad" indices to False
        if len(bad_idx) > 0:
            boolean_mask[bad_idx] = False
    else:
        boolean_mask = np.ones(nints, dtype=bool)

    # Use the boolean mask to select only the "good" frames
    array = array[boolean_mask]
    good_idx = list(np.where(boolean_mask)[0])
    good_idx = [int(e) for e in good_idx]
    return array, good_idx


def save_superbias(superbias, error, dq, instrument='', detector='',
                   subarray='GENERIC', readpatt='ANY',
                   outfile='superbias_jwst_reffiles.fits',
                   author='jwst_reffiles', description='Super Bias Image',
                   pedigree='GROUND', useafter='2000-01-01T00:00:00',
                   history='', fastaxis=-1, slowaxis=2, substrt1=1, substrt2=1,
                   filenames=[], num_ints=0):
    """Saves a CRDS-formatted superbias reference file.

    Parameters
    ----------
    superbias : numpy.ndarray
        The 2D superbias image.

    error : numpy.ndarray
        The 2D superbias error image.

    dq : numpy.ndarray
        The 2D superbias data quality image.

    instrument : str
        CRDS-required instrument for which to use this reference file for.

    detector : str
        CRDS-required detector for which to use this reference file for.

    subarray : str
        CRDS-required subarray for which to use this reference file for.

    readpatt : str
        CRDS-required read pattern for which to use this reference file for.

    outfile : str
        Name of the CRDS-formatted superbias reference file to save the final
        superbias map to.

    author : str
        CRDS-required name of the reference file author, to be placed in the
        referece file header.

    description : str
        CRDS-required description of the reference file, to be placed in the
        reference file header.

    pedigree : str
        CRDS-required pedigree of the data used to create the reference file.

    useafter : str
        CRDS-required date of earliest data with which this referece file
        should be used. (e.g. '2019-04-01T00:00:00').

    history : str
        CRDS-required history section to place in the reference file header.

    fastaxis : int
        CRDS-required fastaxis of the reference file.

    slowaxis : int
        CRDS-required slowaxis of the reference file.

    substrt1 : int
        CRDS-required starting pixel in axis 1 direction.

    substrt2 : int
        CRDS-required starting pixel in axis 2 direction.

    filenames : list
        List of dark current files that were used to generate the reference
        file.

    num_ints : int
        Number of integrations used to construct the superbias array
    """

    s = SuperBiasModel()

    s.data = superbias
    s.err = error
    s.dq = dq
    s.dq_def = [(0, 0, 'GOOD', ''),
                (0, 1, 'DO_NOT_USE', ''),
                (1, dqflags.pixel['UNRELIABLE_BIAS'], 'UNRELIABLE_BIAS', '')]

    s.meta.instrument.name = instrument
    s.meta.instrument.detector = detector
    s.meta.subarray.name = subarray

    if subarray == 'FULL':
        s.meta.subarray.p_subarr = 'FULLP|'


    s.meta.exposure.readpatt = readpatt
    s.meta.author = author
    s.meta.description = description
    s.meta.pedigree = pedigree
    s.meta.useafter = useafter
    s.meta.subarray.fastaxis = fastaxis
    s.meta.subarray.slowaxis = slowaxis
    s.meta.reftype = 'SUPERBIAS'

    yd, xd = superbias.shape
    s.meta.subarray.xstart = substrt1
    s.meta.subarray.xsize = xd
    s.meta.subarray.ystart = substrt2
    s.meta.subarray.ysize = yd

    package_note = ('This file was created using the superbias.py module '
                    'within the jwst_reffiles package.')
    software_dict = {'name': 'jwst_reffiles.superbias.py', 'author': 'STScI',
                     'homepage': 'https://github.com/spacetelescope/jwst_reffiles',
                     'version': '1.0.2'}
    entry = util.create_history_entry(package_note, software=software_dict)
    s.history.append(entry)

    #if isinstance(history, str):
    #    history = [history]

    print('history type is ', type(history))
    print(history)
    print('\n\n')

    #for ele in history:
    entry = util.create_history_entry(history)
    s.history.append(entry)

    # Add the list of input files used to create the superbias reference file
    s.history.append(f'DATA USED (total of {num_ints} integrations):')
    for f in filenames:
        f = os.path.basename(f)
        totlen = len(f)
        div = np.arange(0, totlen, 60)
        for val in div:
            if totlen > (val+60):
                s.history.append(util.create_history_entry(f[val:val+60]))
            else:
                s.history.append(util.create_history_entry(f[val:]))

    #if history != '':
    #    s.history.append(util.create_history_entry(history))

    s.save(outfile, overwrite=True)
    print('Final CRDS-formatted superbias map saved to {}'.format(outfile))

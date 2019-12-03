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
               Pixels with high superbias values are flagged as UNRELIABLE_BIAS
               and DO_NOT_USE in the superbias DQ array.
"""

import os

from astropy.io import fits
from astropy.stats import sigma_clip
from jwst.datamodels import dqflags, SuperBiasModel, util
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import median_filter

def calculate_mean(stack, clipping_sigma=3.0, max_clipping_iters=3):
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

    Returns
    -------
    mean_image : numpy.ndarray
        2D image of the sigma-clipped mean through the input stack.
    """

    clipped = sigma_clip(stack, sigma=clipping_sigma, 
                         maxiters=max_clipping_iters, axis=0)
    mean_image = np.mean(clipped, axis=0)

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
    stddev_image = np.std(clipped, axis=0)

    return stddev_image

def find_unreliable_bias_pixels(superbias, refpix_map, bad_clipping_sigma=3.0, 
                                bad_max_clipping_iters=3, bad_sigma_threshold=5.0, 
                                box_width=3, plot=False, save_tmp=False, 
                                save_filled=False, outfile='superbias_jwst_reffiles.fits'):
    """Identifies pixels with high superbias values (i.e. UNRELIABLE_BIAS 
    pixels).

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

    plot : bool
        Option to save a histogram of the superbias image minus the smoothed 
        superbias image values, with the threshold used for flagging 
        UNRELIABLE_BIAS pixels also shown.

    save_tmp : bool
        Option to save intermediate images that were used to identify 
        UNRELIABLE_BIAS pixels (e.g. the smoothed superbias image).

    save_filled : bool
        Option to save an image where UNRELIABLE_BIAS pixels are replaced by 
        their smoothed values.

    outfile : str
        Name of the CRDS-formatted superbias reference file to save the final
        superbias map to. This is only used in this function to name the  
        optional output files.

    Returns
    -------
    dq : numpy.ndarray
        A 2D DQ image that flags pixels with high superbias values as 
        UNRELIABLE_BIAS (2) and DO_NOT_USE (1).
    """

    # Trim the refpix from the input superbias
    scipix = np.where(refpix_map==0)
    ymin = np.min(scipix[0])
    xmin = np.min(scipix[1])
    ymax = np.max(scipix[0]) + 1
    xmax = np.max(scipix[1]) + 1
    superbias_no_refpix = superbias[ymin:ymax, xmin:xmax].copy()

    # Make initial DQ arrays, with and without refpix
    dq = np.zeros(superbias.shape, dtype=int)
    dq_no_refpix = np.zeros(superbias_no_refpix.shape, dtype=int)

    # Flag pixels with nan/inf superbias values
    dq_no_refpix[~np.isfinite(superbias_no_refpix)] = 3

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

    # Plot a histogram of the difference image values
    if plot:
        plt.hist(diff.flatten(), range=(-thresh*4, thresh*4), bins=80, 
                 color='steelblue')
        plt.axvline(thresh, ls='--', color='red')
        plt.yscale('log')
        plt.ylabel('Number of Pixels')
        plt.xlabel('Difference (superbias - smoothed superbias) [DN]')
        plt.title('# of UNRELIABLE_BIAS pixels: {}'.format(len(dq[dq!=0])))
        plt.savefig('diff_hist.png', dpi=100, bbox_inches='tight')
        plt.clf()

    # Save images used to identify UNRELIABLE_BIAS pixels
    if save_tmp:
        fits.writeto('tmp_superbias_original.fits', 
                     superbias_no_refpix, overwrite=True)
        fits.writeto('tmp_superbias_smoothed.fits', smoothed, 
                     overwrite=True)
        fits.writeto('tmp_superbias_original_minus_smoothed.fits', 
                     diff, overwrite=True)

    # Save an image where bad superbias pixels are replaced by their smoothed
    # values.
    if save_filled:
        outname = outfile.replace('.fits', '_filled.fits')
        superbias_filled_no_refpix = superbias[ymin:ymax, xmin:xmax].copy()
        superbias_filled_no_refpix[dq_no_refpix != 0] = \
            smoothed[dq_no_refpix != 0]
        superbias_filled = superbias.copy()
        superbias_filled[ymin:ymax, xmin:xmax] = superbias_filled_no_refpix
        fits.writeto(outname, superbias_filled, overwrite=True)

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
        readpatt = 'ANY'

    return instrument, detector, subarray, readpatt, fastaxis, slowaxis, substrt1, substrt2

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
                  'will be used when identifying UNRELIABLE_BIAS pixels.')
    except KeyError:
        print('Warning: No PIXELDQ extension found in input file. Assuming '
              'no reference pixels exist (if this isnt true, i.e. if refpix '
              'actually do exist in this file, it may affect the search for '
              'UNRELIABLE_BIAS pixels).')
        data = fits.getdata(filename, 'SCI')[0, 0, :, :]
        refpix_map = np.zeros(data.shape).astype('int')

    return refpix_map

def make_superbias(filenames, clipping_sigma=3.0, max_clipping_iters=3, 
                   bad_clipping_sigma=3.0, bad_max_clipping_iters=3, 
                   bad_sigma_threshold=5.0, box_width=3, plot=False, 
                   save_tmp=False, save_filled=False, save_reffile=True, 
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

    clipping_sigma : float
        Number of sigma to use when sigma-clipping through the 0th frame  
        stack to create the superbias image.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping through   
        the 0th frame stack to create the superbias error image.

    bad_clipping_sigma : float
        Number of sigma to use when sigma-clipping the superbias image to 
        find UNRELIABLE_BIAS pixels.

    bad_max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping the superbias 
        image to find UNRELIABLE_BIAS pixels.

    bad_sigma_threshold : float
        Number of standard deviations above the mean of the superbias image 
        at which a pixel is flagged as UNRELIABLE_BIAS.

    box_width : int
        The box size to use in pixels when smoothing the superbias image 
        to find UNRELIABLE_BIAS pixels.

    plot : bool
        Option to save a histogram of the superbias image minus the smoothed 
        superbias image values (i.e. the image used to identify UNRELIABLE_BIAS
        pixels), with the threshold used for flagging UNRELIABLE_BIAS 
        pixels also shown.

    save_tmp : bool
        Option to save intermediate images that were used to identify 
        UNRELIABLE_BIAS pixels (e.g. the smoothed superbias image).

    save_filled : bool
        Option to save an image where UNRELIABLE_BIAS pixels are replaced by 
        their smoothed values.

    save_reffile : bool
        Option to save the generated superbias map into a CRDS-formatted 
        reference file.

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
    n_ints, n_groups, n_y, n_x = fits.getdata(filenames[0], 'SCI').shape
    stack = np.zeros((len(filenames), n_y, n_x))
    for i,filename in enumerate(filenames):
        stack[i] = fits.getdata(filename, 'SCI')[0, 0, :, :].astype(float)

    # Calculate the superbias as the sigma-clipped mean through the stack
    superbias = calculate_mean(stack, clipping_sigma=clipping_sigma, 
                               max_clipping_iters=max_clipping_iters)

    # Convert masked array to normal numpy array and check for any missing 
    # data.
    superbias = superbias.filled(fill_value=np.nan)
    n_bad = len(superbias[~np.isfinite(superbias)])
    if n_bad > 0:
        print('Warning: Superbias file has {} nan/inf pixels.'.format(n_bad))

    if save_reffile:

        # Get a map of the refpix pixels, so they can be excluded when 
        # creating the superbias DQ extension
        refpix_map = make_refpix_map(filenames[0])

        # Inputs listed as None in the config file are read in as strings.
        # Change these to NoneType objects.
        subarray = none_check(subarray)
        readpatt = none_check(readpatt)

        # CRDS doesnt allow descriptions over 65 characters
        if len(description) > 65:
            raise ValueError('Description cannot exceed 65 characters.')

        # Get info needed by CRDS to put into the final superbias reference file
        instrument, detector, subarray, readpatt, fastaxis, slowaxis, substrt1, substrt2 = \
            get_crds_info(filenames[0], subarray, readpatt)

        # Calculate the superbias error as the sigma-clipped stddev through 
        # the stack.
        error = calculate_stddev(stack, clipping_sigma=clipping_sigma, 
                                 max_clipping_iters=max_clipping_iters)

        # Convert masked array to normal numpy array and check for any missing 
        # data.
        error = error.filled(fill_value=np.nan)
        n_bad = len(error[~np.isfinite(error)])
        if n_bad > 0:
            print('Warning: Superbias Error array has {} nan/inf pixels.'
                  .format(n_bad))

        # Flag pixels with high superbias values as UNRELIABLE_BIAS
        dq = find_unreliable_bias_pixels(superbias, refpix_map, 
                                         bad_clipping_sigma=bad_clipping_sigma,
                                         bad_max_clipping_iters=bad_max_clipping_iters,
                                         bad_sigma_threshold=bad_sigma_threshold,
                                         box_width=box_width,
                                         plot=plot, save_tmp=save_tmp, 
                                         save_filled=save_filled, outfile=outfile)

        # Save CRDS-formatted superbias reference file
        save_superbias(superbias, error, dq, instrument=instrument, 
                       detector=detector, subarray=subarray, readpatt=readpatt, 
                       outfile=outfile, author=author, description=description, 
                       pedigree=pedigree, useafter=useafter, history=history, 
                       fastaxis=fastaxis, slowaxis=slowaxis, substrt1=substrt1, 
                       substrt2=substrt2, filenames=filenames)

    return superbias

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

def save_superbias(superbias, error, dq, instrument='', detector='', 
                   subarray='GENERIC', readpatt='ANY', 
                   outfile='superbias_jwst_reffiles.fits', 
                   author='jwst_reffiles', description='Super Bias Image', 
                   pedigree='GROUND', useafter='2000-01-01T00:00:00', 
                   history='', fastaxis=-1, slowaxis=2, substrt1=1, substrt2=1, 
                   filenames=[]):
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
    """

    s = SuperBiasModel()
    
    s.data = superbias
    s.err = error
    s.dq = dq
    s.dq_def = [(0, 0, 'GOOD', ''),
                (0, 1, 'DO_NOT_USE', ''),
                (1, 2, 'UNRELIABLE_BIAS', '')]

    s.meta.instrument.name = instrument
    s.meta.instrument.detector = detector
    s.meta.subarray.name = subarray
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
                     'version': '0.0.0'}
    entry = util.create_history_entry(package_note, software=software_dict)
    s.history.append(entry)

    # Add the list of input files used to create the superbias reference file
    s.history.append('DATA USED:')
    for f in filenames:
        f = os.path.basename(f)
        totlen = len(f)
        div = np.arange(0, totlen, 60)
        for val in div:
            if totlen > (val+60):
                s.history.append(util.create_history_entry(f[val:val+60]))
            else:
                s.history.append(util.create_history_entry(f[val:]))
    
    if history != '':
        s.history.append(util.create_history_entry(history))
    
    s.save(outfile, overwrite=True)
    print('Final CRDS-formatted superbias map saved to {}'.format(outfile))

#! /usr/bin/env python

"""This module creates a readnoise reference file that can be used 
in the JWST calibration pipeline.

Author
------
     - Ben Sunnquist

Use
---
    This module can be imported and used as such:
    ::
        from jwst_reffiles.readnoise import readnoise
        readnoise.make_readnoise(arguments)

Notes
-----
    Overview:
    Inputs: A list of calibrated dark current ramps

    Algorithm: Create CDS images for each input ramp. Stack all of these 
    images from each input ramp together. The readnoise is the sigma-clipped 
    standard deviation through this master CDS image stack.
"""

from astropy.io import fits
from astropy.stats import sigma_clip
from multiprocessing import Pool
import numpy as np

def calculate_stddev(stack, clipping_sigma=3, max_clipping_iters=5):
    """Calculates the sigma-clipped standard deviation through a stack
    of images.

    Parameters
    ----------
    stack : numpy.ndarray
        A 3D stack of images.
    
    clipping_sigma : int
        Number of sigmas to use when sigma-clipping the input stack.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping the input 
        stack.

    Returns
    -------
    stddev : numpy.ndarray
        2D image of the sigma-clipped standard deviation through the 
        input stack.
    """

    clipped = sigma_clip(stack, sigma=clipping_sigma, 
                         maxiters=max_clipping_iters, axis=0)
    stddev = np.std(clipped, axis=0)

    return stddev

def make_cds_stack(data, group_diff_type='independent'):
    """Creates a stack of CDS images (group difference images) using the 
    input ramp data, combining multiple integrations if necessary.

    Parameters
    ----------
    data : numpy.ndarray
        The input ramp data. The data shape is assumed to be a 4D array in 
        DMS format (integration, group, y, x).
    
    group_diff_type : str
        The method for calculating group differences. Options are:
        ``independent``: Each groups is only differenced once (e.g. 6-5, 4-3, 
                         2-1)
        ``consecutive``: Each group is differenced to its neighbors (e.g. 
                         4-3, 3-2, 2-1)

    Returns
    -------
    cds_stack : numpy.ndarray
        A 3D stack of the group difference images.
    """

    n_ints, n_groups, n_y, n_x = data.shape
    for integration in range(n_ints):
        if group_diff_type == 'independent':
            cds = data[integration, 1::2, :, :] - data[integration, ::2, :, :]
        elif group_diff_type == 'consecutive':
            cds = data[integration, 1:, :, :] - data[integration, 0:-1, :, :]
        else:
            raise ValueError('Unknown group difference option: {}.'
                             .format(group_diff_type))
        if integration == 0:
            cds_stack = cds
        else:
            cds_stack = np.concatenate((cds_stack, cds), axis=0)

    return cds_stack

def make_readnoise(filenames, group_diff_type='independent', clipping_sigma=3, 
                   max_clipping_iters=5, nproc=6, slice_width=50):
    """The main function. Creates a readnoise reference file using the input 
    dark current ramps. See module docstring for more details.
    
    Parameters
    ----------
    filenames : list
        List of dark current files. These should be calibrated ramp images. 
        The data shape in these images is assumed to be a 4D array in DMS 
        format (integration, group, y, x).

    group_diff_type : str
        The method for calculating group differences. Options are:
        ``independent``: Each groups is only differenced once (e.g. 6-5, 4-3, 
                         2-1)
        ``consecutive``: Each group is differenced to its neighbors (e.g. 
                         4-3, 3-2, 2-1)

    clipping_sigma : int
        Number of sigma to use when sigma-clipping the 3D array of
        CDS images to find the readnoise.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping the 3D array 
        of CDS images to find the readnoise.

    nproc : int
        The number of processes to use during multiprocessing.

    slice_width : int
        The width (in pixels) of the image slice to use during 
        multiprocessing. The readnoise of each slice is calculatd separately 
        during multiprocessing and combined together at the end of 
        processing.
    """

    # Find the column indexes for each image slice to use during 
    # multiprocessing.
    n_x = fits.getdata(filenames[0], 'SCI').shape[3]
    columns = list(np.arange(n_x)[::slice_width])
    n_cols = len(columns)

    # Find the readnoise using the input ramp files; do this for each slice  
    # separately to allow for multiprocessing and avoiding memory issues.
    p = Pool(nproc)
    files = [filenames] * n_cols
    group_diff_types = [group_diff_type] * n_cols
    sigmas = [clipping_sigma] * n_cols
    iters = [max_clipping_iters] * n_cols
    slice_widths = [slice_width] * n_cols
    readnoise = p.map(wrapper_readnoise_by_slice, 
                      zip(files, group_diff_types, sigmas, iters, columns, 
                          slice_widths)
                      )
    readnoise = np.concatenate(readnoise, axis=1)
    p.close()

    # Convert masked array to normal numpy array and check for any missing 
    # data.
    readnoise = readnoise.filled(fill_value=np.nan)
    n_nans = len(readnoise[~np.isfinite(readnoise)])
    if n_nans > 0:
        print('Warning: Readnoise file has {} nan pixels.'.format(n_nans))
    fits.writeto('readnoise_new.fits', readnoise, overwrite=True)

def readnoise_by_slice(filenames, group_diff_type='independent', 
                       clipping_sigma=3, max_clipping_iters=5, column=0,
                       slice_width=50):
    """Calculates the readnoise for a given slice in the input dark file 
    ramps. Useful for multiprocessing and avoiding memory issues for large 
    image stacks.

    Parameters
    ----------
    filenames : list
        List of dark current files. These should be calibrated ramp images  
        with the same image shape.

    group_diff_type : str
        The method for calculating group differences. Options are:
        ``independent``: Each groups is only differenced once (e.g. 6-5, 4-3, 
                         2-1)
        ``consecutive``: Each group is differenced to its neighbors (e.g. 
                         4-3, 3-2, 2-1)

    clipping_sigma : int
        Number of sigma to use when sigma-clipping the 3D array of
        CDS images to find the readnoise.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping the 3D array 
        of CDS images to find the readnoise.

    column : int
        The index of the starting image column for each slice.

    slice_width : int
        The width of the slice in pixels.

    Returns
    -------
    readnoise : numpy.ndarray
        2D image of the calculated readnoise.
    """

    for i,filename in enumerate(filenames):

        # Get image data for the given slice; if the slice goes outside of 
        # the image border, the slice just goes up to the end of the image.
        data = fits.getdata(filename, 'SCI')
        data = data[:, :, :, column:column+slice_width]

        # Remove first 5 groups and last group for MIRI to avoid reset/rscd 
        # effects.
        instrument = fits.getheader(filename)['INSTRUME']
        if instrument == 'MIRI':
            data = data[:, 5:-1, :, :]  # need to test this

        # Create the CDS stack for the image (combining multiple 
        # integrations if necessary) and add it to the master CDS stack  
        # containing the CDS stacks for all images and integrations.
        cds_stack = make_cds_stack(data, group_diff_type=group_diff_type)
        if i == 0:
            master_cds_stack = cds_stack
        else:
            master_cds_stack = np.concatenate((master_cds_stack, cds_stack), 
                                              axis=0)

    # Find the readnoise for this column
    readnoise = calculate_stddev(master_cds_stack, 
                                 clipping_sigma=clipping_sigma, 
                                 max_clipping_iters=max_clipping_iters)

    return readnoise

def wrapper_readnoise_by_slice(args):
    """A wrapper around the readnoise_by_slice function to allow for 
    multiprocessing.
    
    Parameters
    ----------
    args : tuple
        A tuple containing the input arguments for the readnoise_by_slice 
        function. See readnoise_by_slice docstring for more details.

    Returns
    -------
    readnoise : numpy.ndarray
        2D image of the calculated readnoise (i.e. the output from the 
        readnoise_by_slice function).
    """

    return readnoise_by_slice(*args)

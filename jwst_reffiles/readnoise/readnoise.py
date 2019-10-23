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

    Algorithm:
        1. 
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

def make_cds_stack(filename, group_diff_type='independent'):
    """Creates a stack of CDS images (group difference images), combining 
    multiple integrations if necessary.

    Parameters
    ----------
    filename : str
        The calibrated dark current ramp file. The data shape in this file is 
        assumed to be a 4D array in DMS format (integration, group, y, x).
    
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

    print('Creating CDS stack for {}'.format(filename))
    data = fits.getdata(filename, 'SCI')
    instrument = fits.getheader(filename)['INSTRUME']
    
    # Remove first 5 groups from MIRI data
    if instrument == 'MIRI':
        data = data[:, 5:, :, :]

    # Create a 3D stack of CDS images that combines all integrations together
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
                   max_clipping_iters=5, nproc=6):
    """The main function. See module docstring for more details.
    
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

    nproc : int
        The number of processes to use during multiprocessing.
    """

    # Create a master CDS stack of all input files and integrations
    p = Pool(nproc)
    group_diff_types = [group_diff_type] * len(filenames)
    cds_stacks = p.map(wrapper_make_cds_stack, zip(filenames, group_diff_types))
    master_cds_stack = np.concatenate(cds_stacks, axis=0)
    p.close()

    master_cds_stack_shape = master_cds_stack.shape
    print('Master CDS stack shape: {}'.format(master_cds_stack_shape))

    # Calculate the readnoise based on the master CDS stack. Calculate this
    # column-by-column to avoid memory issues.
    print('Calculating sigma-clipped stats through master CDS stack...')
    p = Pool(nproc)
    column_stacks = np.split(master_cds_stack, master_cds_stack_shape[2], axis=2)
    sigmas = [clipping_sigma] * len(column_stacks)
    iters = [max_clipping_iters] * len(column_stacks)
    column_stddevs = p.map(wrapper_calculate_stddev, zip(column_stacks, sigmas, iters))
    readnoise = np.concatenate(column_stddevs, axis=1)
    p.close()

    # Convert masked array to normal numpy array and check for any missing data
    readnoise = readnoise.filled(fill_value=np.nan)
    n_nans = len(readnoise[~np.isfinite(readnoise)])
    if n_nans > 0:
        print('Warning: Readnoise file has {} nans.'.format(n_nans))
    fits.writeto('readnoise.fits', readnoise, overwrite=True)

def wrapper_calculate_stddev(args):
    """A wrapper around the calculate_stddev() function to allow for 
    multiprocessing.
    
    Parameters
    ----------
    args : tuple
        A tuple containing the input arguments for the calculate_stddev() 
        function (3D image stack, sigma for sigma-clipping, and maximum 
        clipping iterations for sigma-clipping). See calculate_stddev()
        docstring for more details.
    
    Returns
    -------
    stddev : numpy.ndarray
        2D image of the sigma-clipped standard deviation through the 
        input stack (i.e. the output from the calculate_stddev() function).
    """

    return calculate_stddev(*args)

def wrapper_make_cds_stack(args):
    """A wrapper around the make_cds_stack() function to allow for 
    multiprocessing.
    
    Parameters
    ----------
    args : tuple
        A tuple containing the input arguments for the make_cds_stack() 
        function (the calibrated dark ramp filename and the method for 
        calculating the group differences). See make_cds_stack() docstring 
        for more details.

    Returns
    -------
    cds_stack : numpy.ndarray
        A 3D stack of the group difference images (i.e. the output from  
        the make_cds_stack() function).
    """

    return make_cds_stack(*args)

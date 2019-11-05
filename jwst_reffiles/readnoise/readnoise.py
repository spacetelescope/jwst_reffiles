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
    Inputs: A list of calibrated (refpix-corrected) dark current ramps

    Algorithm:

    Method 1: 
    Create CDS images for each input ramp. Stack all of these 
    images from each input ramp together. The readnoise is the sigma-clipped 
    standard deviation through this master CDS image stack.
    
    Method 2:
    Calculate the readnoise in each ramp individually, and then average
    these readnoise images together to produce a final readnoise map.
"""

from astropy.io import fits
from astropy.stats import sigma_clip
from multiprocessing import Pool
import numpy as np

from jwst.datamodels import ReadnoiseModel

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
        ``independent``: Each group is only differenced once (e.g. 6-5, 4-3, 
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

def make_readnoise(filenames, method='stack', group_diff_type='independent', 
                   clipping_sigma=3.0, max_clipping_iters=3, nproc=1, 
                   slice_width=50, outfile='readnoise_jwst_reffiles.fits'):
    """The main function. Creates a readnoise reference file using the input 
    dark current ramps. See module docstring for more details.
    
    Parameters
    ----------
    filenames : list
        List of dark current files. These should be calibrated ramp images. 
        The data shape in these images is assumed to be a 4D array in DMS 
        format (integration, group, y, x).

    method : str
        The method to use when calculating the readnoise. Options are:
        ``stack``: Creates a master stack of all CDS images by combining the 
                   CDS images for each ramp and integration together. The 
                   final readnoise map is calculated by taking the stddev 
                   through this master stack.
        ``ramp``: Calculates the readnoise in each ramp individually. The 
                  final readnoise map is calculated by averaging these 
                  individual readnoise maps together.

    group_diff_type : str
        The method for calculating group differences. Options are:
        ``independent``: Each group is only differenced once (e.g. 6-5, 4-3, 
                         2-1).
        ``consecutive``: Each group is differenced to its neighbors (e.g. 
                         4-3, 3-2, 2-1).

    clipping_sigma : float
        Number of sigma to use when sigma-clipping.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping.

    nproc : int
        The number of processes to use during multiprocessing.

    slice_width : int
        The width (in pixels) of the image slice to use during 
        multiprocessing. The readnoise of each slice is calculatd separately 
        during multiprocessing and combined together at the end of 
        processing. Only relevant if method==stack.

    outfile : str
        The filename to give to the output readnoise image.
    """

    if method == 'stack':

        # Find the column indexes for each image slice to use during 
        # multiprocessing.
        n_x = fits.getdata(filenames[0], 'SCI').shape[3]
        columns = list(np.arange(n_x)[::slice_width])
        n_cols = len(columns)

        # Find the readnoise by taking the stddev through a master stack of 
        # CDS images that incorporates every input dark ramp and integration; 
        # do this in slices to allow for multiprocessing and avoiding memory 
        # issues.
        print('Calculating the readnoise in {} separate image slices...'
              .format(n_cols))
        p = Pool(nproc)
        files = [filenames] * n_cols
        group_diff_types = [group_diff_type] * n_cols
        sigmas = [clipping_sigma] * n_cols
        iters = [max_clipping_iters] * n_cols
        slice_widths = [slice_width] * n_cols
        readnoise = p.map(wrapper_readnoise_by_slice, 
                          zip(files, group_diff_types, sigmas, iters, 
                              columns, slice_widths)
                          )
        readnoise = np.concatenate(readnoise, axis=1)
        p.close()

    elif method == 'ramp':

        # Create a 3D stack of readnoise images, one for each input ramp
        print('Calculating the readnoise in each of the {} input ramps '
              'individually...'.format(len(filenames)))
        p = Pool(nproc)
        n_files = len(filenames)
        group_diff_types = [group_diff_type] * n_files
        sigmas = [clipping_sigma] * n_files
        iters = [max_clipping_iters] * n_files
        readnoise_stack = p.map(wrapper_readnoise_by_ramp,
                                zip(filenames, group_diff_types, sigmas, 
                                    iters)
                                )
        p.close()

        # Create the final readnoise map by averaging the individual  
        # readnoise images.
        print('Combining the readnoises for all {} individual ramps '
              'together...'.format(len(filenames)))
        readnoise = calculate_mean(np.array(readnoise_stack), 
                                   clipping_sigma=clipping_sigma, 
                                   max_clipping_iters=max_clipping_iters)

    else:
        raise ValueError('Unknown readnoise method: {}'.format(method))

    # Convert masked array to normal numpy array and check for any missing 
    # data.
    readnoise = readnoise.filled(fill_value=np.nan)
    n_nans = len(readnoise[~np.isfinite(readnoise)])
    if n_nans > 0:
        print('Warning: Readnoise file has {} nan pixels.'.format(n_nans))
    
    # Save the final readnoise reference file
    header = fits.getheader(filenames[0])
    instrument = header['INSTRUME']
    detector = header['DETECTOR']
    subarray = header['SUBARRAY']
    save_readnoise(readnoise, instrument, detector, subarray, outfile=outfile)

def readnoise_by_ramp(filename, group_diff_type='independent', 
                      clipping_sigma=3.0, max_clipping_iters=3):
    """Calculates the readnoise for the given input dark current ramp.
    
    Parameters
    ----------
    filename : str
        The dark current ramp file. The data shape in this image is assumed  
        to be a 4D array in DMS format (integration, group, y, x).

    group_diff_type : str
        The method for calculating group differences. Options are:
        ``independent``: Each group is only differenced once (e.g. 6-5, 4-3, 
                         2-1)
        ``consecutive``: Each group is differenced to its neighbors (e.g. 
                         4-3, 3-2, 2-1)

    clipping_sigma : float
        Number of sigma to use when sigma-clipping.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping.
    """

    print('Calculating readnoise for {}'.format(filename))

    # Get the ramp data; remove first 5 groups and last group for MIRI to 
    # avoid reset/rscd effects.
    data = fits.getdata(filename, 'SCI', uint=False)
    instrument = fits.getheader(filename)['INSTRUME']
    if instrument == 'MIRI':
        data = data[:, 5:-1, :, :]

    # Create a CDS stack for the input ramp (combining multiple integrations 
    # if necessary).
    cds_stack = make_cds_stack(data, group_diff_type=group_diff_type)

    # Calculate the readnoise
    readnoise = calculate_stddev(cds_stack, clipping_sigma=clipping_sigma, 
                                 max_clipping_iters=max_clipping_iters)

    return readnoise

def readnoise_by_slice(filenames, group_diff_type='independent', 
                       clipping_sigma=3.0, max_clipping_iters=3, column=0,
                       slice_width=50):
    """Calculates the readnoise for a given slice in the input dark file 
    ramps. Useful for multiprocessing and avoiding memory issues for large 
    image stacks.

    Parameters
    ----------
    filenames : list
        List of dark current files. The data shape in these images is assumed 
        to be a 4D array in DMS format (integration, group, y, x).

    group_diff_type : str
        The method for calculating group differences. Options are:
        ``independent``: Each group is only differenced once (e.g. 6-5, 4-3, 
                         2-1)
        ``consecutive``: Each group is differenced to its neighbors (e.g. 
                         4-3, 3-2, 2-1)

    clipping_sigma : float
        Number of sigma to use when sigma-clipping.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping.

    column : int
        The index of the starting image column for each slice.

    slice_width : int
        The width of the slice in pixels.

    Returns
    -------
    readnoise : numpy.ndarray
        2D image of the calculated readnoise.
    """

    print('Calculating readnoise in image slice [{}:{}]'.format(column, 
          column+slice_width))

    for i,filename in enumerate(filenames):

        # Get image data for the given slice; if the slice goes outside of 
        # the image border, the slice just goes up to the end of the image.
        data = fits.getdata(filename, 'SCI', uint=False)
        data = data[:, :, :, column:column+slice_width]

        # Remove first 5 groups and last group for MIRI to avoid reset/rscd 
        # effects.
        instrument = fits.getheader(filename)['INSTRUME']
        if instrument == 'MIRI':
            data = data[:, 5:-1, :, :]

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

def save_readnoise(readnoise, instrument='', detector='', subarray='GENERIC', 
                   readpatt='ANY', outfile='readnoise_jwst_reffiles.fits'):
    """Saves a readnoise image that can be used as a reference file in the 
    JWST calibration pipeline.

    Parameters
    ----------
    readnoise : numpy.ndarray
        The 2D readnoise image.

    instrument : str
        The instrument to use this reference file for.

    detector : str
        The detector to use this reference file for.

    subarray : str
        The subarray to use this reference file for.

    readpatt : str
        The readout pattern to use this reference file for.

    outfile : str
        The filename to give to the output readnoise image.
    """

    r = ReadnoiseModel()
    r.data = readnoise
    r.meta.instrument.name = instrument
    r.meta.instrument.detector = detector
    r.meta.subarray.name = subarray
    r.meta.exposure.readpatt = readpatt
    r.meta.description = 'Readnoise image'
    r.meta.useafter = '2000-01-01T00:00:00'
    r.save(outfile)

    print('Final readnoise map saved to {}'.format(outfile))

def wrapper_readnoise_by_ramp(args):
    """A wrapper around the readnoise_by_ramp function to allow for 
    multiprocessing.
    
    Parameters
    ----------
    args : tuple
        A tuple containing the input arguments for the readnoise_by_ramp 
        function. See readnoise_by_ramp docstring for more details.

    Returns
    -------
    readnoise : numpy.ndarray
        2D image of the calculated readnoise (i.e. the output from the 
        readnoise_by_ramp function).
    """

    return readnoise_by_ramp(*args)

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

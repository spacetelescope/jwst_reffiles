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
    Inputs: A list of dark current ramps

    Algorithm:

    Method 1 (method='stack'): 
    Create CDS images for each input ramp. Stack all of these 
    images from each input ramp together. The readnoise is the sigma-clipped 
    standard deviation through this master CDS image stack.
    
    Method 2 (method='ramp'):
    Calculate the readnoise in each ramp individually, and then average
    these readnoise images together to produce a final readnoise map.
"""

from multiprocessing import Pool
import os

from astropy.io import fits
from astropy.stats import sigma_clip
from jwst.datamodels import ReadnoiseModel, util
import numpy as np

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

def get_crds_info(filename, subarray, readpatt):
    """Get CRDS-required header information to put into the final readnoise 
    reference file.

    Parameters
    ----------
    filename : str
        Path to one of the files being used to generate the readnoise 
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
                   slice_width=50, single_value=False, 
                   outfile='readnoise_jwst_reffiles.fits', 
                   author='jwst_reffiles', description='CDS Noise Image', 
                   pedigree='GROUND', useafter='2015-10-01T00:00:00', 
                   history='', subarray=None, readpatt=None, save_tmp=False):
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
        multiprocessing. The readnoise of each slice is calculated separately 
        during multiprocessing and combined together at the end of 
        processing. Only relevant if method==stack.

    single_value : bool
        Option to use a single readnoise value (the average of all readnoise 
        values) for all pixels.

    outfile : str
        Name of the CRDS-formatted readnoise reference file to save the final
        readnoise map to.

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

    subarray : str
        CRDS-required subarray for which to use this reference file for.

    readpatt : str
        CRDS-required read pattern for which to use this reference file for.

    save_tmp : bool
        Option to save the final readnoise map before turning it into CRDS 
        format. This is useful if the CRDS transformation fails; in this  
        scenario, you won't lose all of the final readnoise results so all 
        of the previous processing wasn't for nothing.
    """

    # Inputs listed as None in the config file are read in as strings.
    # Change these to NoneType objects.
    subarray = none_check(subarray)
    readpatt = none_check(readpatt)

    # CRDS doesnt allow descriptions over 65 characters
    if len(description) > 65:
        raise ValueError('Description cannot exceed 65 characters.')

    # Get info needed by CRDS to put into the final readnoise reference file
    instrument, detector, subarray, readpatt, fastaxis, slowaxis, substrt1, substrt2 = \
        get_crds_info(filenames[0], subarray, readpatt)

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
        p.join()

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
        p.join()

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

    # Option to use a single readnoise value for all pixels
    if single_value:
        readnoise = use_single_value(readnoise, instrument=instrument,
                                     detector=detector, subarray=subarray, 
                                     clipping_sigma=clipping_sigma, 
                                     max_clipping_iters=max_clipping_iters)
    
    # Save final readnoise map before turning it into CRDS format in case
    # that step has issues (so all of the previous work isnt lost).
    if save_tmp:
        outfile_temp = outfile.replace('.fits','_temp.fits')
        fits.writeto(outfile_temp, readnoise, overwrite=True)
        print('Readnoise map saved to {}'.format(outfile_temp))
        print('Creating CRDS-formatted version of {}...'.format(outfile_temp))

    # Save the final readnoise reference file in CRDS format
    save_readnoise(readnoise, instrument=instrument, detector=detector, 
                   subarray=subarray, readpatt=readpatt, outfile=outfile, 
                   author=author, description=description, pedigree=pedigree, 
                   useafter=useafter, history=history, fastaxis=fastaxis, 
                   slowaxis=slowaxis, substrt1=substrt1, substrt2=substrt2, 
                   filenames=filenames)

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
                   readpatt='ANY', outfile='readnoise_jwst_reffiles.fits',
                   author='jwst_reffiles', description='CDS Noise Image', 
                   pedigree='GROUND', useafter='2015-10-01T00:00:00', 
                   history='', fastaxis=-1, slowaxis=2, substrt1=1, substrt2=1, 
                   filenames=[]):
    """Saves a CRDS-formatted readnoise reference file.

    Parameters
    ----------
    readnoise : numpy.ndarray
        The 2D readnoise image.

    instrument : str
        CRDS-required instrument for which to use this reference file for.

    detector : str
        CRDS-required detector for which to use this reference file for.

    subarray : str
        CRDS-required subarray for which to use this reference file for.

    readpatt : str
        CRDS-required read pattern for which to use this reference file for.

    outfile : str
        Name of the CRDS-formatted readnoise reference file to save the final
        readnoise map to.

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

    r = ReadnoiseModel()
    
    r.data = readnoise
    r.meta.bunit_data = 'DN'
    r.meta.instrument.name = instrument
    r.meta.instrument.detector = detector
    r.meta.subarray.name = subarray
    r.meta.exposure.readpatt = readpatt
    r.meta.author = author
    r.meta.description = description
    r.meta.pedigree = pedigree
    r.meta.useafter = useafter
    r.meta.subarray.fastaxis = fastaxis
    r.meta.subarray.slowaxis = slowaxis
    r.meta.reftype = 'READNOISE'

    yd, xd = readnoise.shape
    r.meta.subarray.xstart = substrt1
    r.meta.subarray.xsize = xd
    r.meta.subarray.ystart = substrt2
    r.meta.subarray.ysize = yd

    package_note = ('This file was created using the readnoise.py module '
                    'within the jwst_reffiles package.')
    software_dict = {'name': 'jwst_reffiles.readnoise.py', 'author': 'STScI',
                     'homepage': 'https://github.com/spacetelescope/jwst_reffiles',
                     'version': '0.0.0'}
    entry = util.create_history_entry(package_note, software=software_dict)
    r.history.append(entry)

    # Add the list of input files used to create the readnoise reference file
    r.history.append('DATA USED:')
    for f in filenames:
        f = os.path.basename(f)
        totlen = len(f)
        div = np.arange(0, totlen, 60)
        for val in div:
            if totlen > (val+60):
                r.history.append(util.create_history_entry(f[val:val+60]))
            else:
                r.history.append(util.create_history_entry(f[val:]))
    
    if history != '':
        r.history.append(util.create_history_entry(history))
    
    r.save(outfile, overwrite=True)
    print('Final CRDS-formatted readnoise map saved to {}'.format(outfile))

def use_single_value(readnoise, instrument, detector, subarray, 
                     clipping_sigma=3.0, max_clipping_iters=3):
    """Replaces all readnoise image values with the average of all readnoise
    values.

    Parameters
    ----------
    readnoise : numpy.ndarray
        The 2D readnoise image.

    instrument : str
        CRDS-required instrument for which to use this reference file for.

    detector : str
        CRDS-required detector for which to use this reference file for.

    subarray : str
        CRDS-required subarray for which to use this reference file for.

    clipping_sigma : float
        Number of sigma to use when sigma-clipping.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping.

    Returns
    -------
    readnoise_single_value : numpy.ndarray
        2D readnoise image where all of the readnoise values are the same 
        (the average of all of the original readnoise values).
    """

    readnoise_single_value = np.zeros(readnoise.shape)

    # For NIRCam full-frame darks, this option is used to create a readnoise 
    # reffile to be used for those subarrays without refpix. Therefore, 
    # only use the readnoise values within the output used to readout the 
    # subarrays. This is always output 7, i.e. x,y=[0:511,0:2047] in detector 
    # coordinates.
    if ((instrument == 'NIRCAM') & 
        ((subarray == 'FULL') | (subarray == 'GENERIC'))):
        print('Only using Output 7 values to find single readnoise value '
              '(i.e. assuming this readnoise reffile is being used for '
              'those subarrays without reference pixels).')

        # Put the data into detector coordinates
        if detector in ['NRCA1', 'NRCA3', 'NRCALONG', 'NRCB2', 'NRCB4']:
            readnoise = np.flip(readnoise, axis=1)
        elif detector in ['NRCA2', 'NRCA4', 'NRCB1', 'NRCB3', 'NRCBLONG']:
            readnoise = np.flip(readnoise, axis=0)
        else:
            print('Detector {} not recognized - keeping coordinates as is.'
                  .format(detector))
        
        # Get only the output 7 values (i.e. x,y=[0:511,0:2047] in detector 
        # coordinates).
        readnoise, *_ = np.split(readnoise, 4, axis=1)

    # Find the average of all readnoise values, and use it as the single 
    # readnoise value in the readnoise reffile.
    clipped = sigma_clip(readnoise, sigma=clipping_sigma, 
                         maxiters=max_clipping_iters)
    readnoise_single_value[:, :] = np.mean(clipped)

    return readnoise_single_value

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

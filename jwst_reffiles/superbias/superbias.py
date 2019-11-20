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
               The superbias error is the sigma-clipped stddev of all 0th 
               frames.
"""

import os

from astropy.io import fits
from astropy.stats import sigma_clip
from jwst.datamodels import dqflags, SuperBiasModel, util
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

def create_dqdef():
    """Create the DQ definition data needed to populate the final superbias 
    reference file.

    Returns
    -------
    definitions : list
        Bad pixel bit definitions
    """

    definitions = []
    standard_defs = dqflags.pixel
    for bitname in standard_defs:
        bitvalue = standard_defs[bitname]
        if bitvalue != 0:
            bitnumber = np.uint8(np.log(bitvalue)/np.log(2))
        else:
            bitnumber = 0
        newrow = (bitnumber, bitvalue, bitname, '')
        definitions.append(newrow)
    
    return definitions

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
    """

    header = fits.getheader(filename)
    instrument = header['INSTRUME']
    detector = header['DETECTOR']
    fastaxis = header['FASTAXIS']
    slowaxis = header['SLOWAXIS']

    if subarray is None:
        subarray = header['SUBARRAY']
    if readpatt is None:
        readpatt = 'ANY'

    return instrument, detector, subarray, readpatt, fastaxis, slowaxis

def make_superbias(filenames, clipping_sigma=3.0, max_clipping_iters=3, 
                   save_reffile=True, outfile='superbias_jwst_reffiles.fits', 
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
        Number of sigma to use when sigma-clipping.

    max_clipping_iters : int
        Maximum number of iterations to use when sigma-clipping.

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
    n_nans = len(superbias[~np.isfinite(superbias)])
    if n_nans > 0:
        print('Warning: Superbias file has {} nan pixels.'.format(n_nans))

    if save_reffile:

        # Inputs listed as None in the config file are read in as strings.
        # Change these to NoneType objects.
        subarray = none_check(subarray)
        readpatt = none_check(readpatt)

        # CRDS doesnt allow descriptions over 65 characters
        if len(description) > 65:
            raise ValueError('Description cannot exceed 65 characters.')

        # Get info needed by CRDS to put into the final superbias reference file
        instrument, detector, subarray, readpatt, fastaxis, slowaxis = \
            get_crds_info(filenames[0], subarray, readpatt)

        # Calculate the superbias error as the sigma-clipped stddev through 
        # the stack.
        error = calculate_stddev(stack, clipping_sigma=clipping_sigma, 
                                 max_clipping_iters=max_clipping_iters)

        # Convert masked array to normal numpy array and check for any missing 
        # data.
        error = error.filled(fill_value=np.nan)
        n_nans = len(error[~np.isfinite(error)])
        if n_nans > 0:
            print('Warning: Superbias Error array has {} nan pixels.'
                  .format(n_nans))

        # Use empty DQ array
        dq = np.zeros((n_y, n_x))
        
        # Save CRDS-formatted superbias reference file
        save_superbias(superbias, error, dq, instrument=instrument, 
                       detector=detector, subarray=subarray, readpatt=readpatt, 
                       outfile=outfile, author=author, description=description, 
                       pedigree=pedigree, useafter=useafter, history=history, 
                       fastaxis=fastaxis, slowaxis=slowaxis, 
                       filenames=filenames)

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
                   history='', fastaxis=-1, slowaxis=2, filenames=[]):
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

    filenames : list
        List of dark current files that were used to generate the reference 
        file.
    """

    s = SuperBiasModel()
    
    s.data = superbias
    s.err = error
    s.dq = dq
    s.dq_def = create_dqdef()

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
    s.meta.subarray.xstart = 1
    s.meta.subarray.xsize = xd
    s.meta.subarray.ystart = 1
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

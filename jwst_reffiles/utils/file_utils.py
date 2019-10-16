#! /usr/bin/env python

"""This module contains general functions related to working with FITS files

Author
------

Bryan Hilbert

"""


def get_file_header(file_name, extension='SCI'):
    """Return the header associated with ``extension`` of ``file_name``

    Parameters
    ----------
    file_name : str
        Name of FITS file

    extension : str
        Name of the FITS file extension containing the header

    Returns
    -------
    header : astropy.io.fits.header
        Header from the ``extension`` extension of ``file_name``
    """
    with fits.open(file_name) as hdulist:
        header = hdulist[extension].header
    return header


def read_ramp(file_name, integ_number=None, min_group=0, max_group=None, extension='SCI'):
    """Read in the science extension from a group of slope images

    Parameters
    ----------
    filename : list
        Name of FITS file containing the data

    integ_number : int
        Integration number of the exposure to be read in. If None, all
        integrations are read in

    min_group : int
        Minimum group number to be read in.

    max_group : int
        Maximum group number to be read in. If None then all groups from
        ``min_group`` to the end are read in.

    extension : str
        Extension name within the FITS file to be read in

    Returns
    -------
    data : numpy.ndarray
        Array of data
    """
    if not os.path.isfile(file_name):
        raise FileNotFoundError('ERROR: Input file {} does not exist'.format(file_name))

    with fits.open(file_name) as hdulist:
        if integ_number is None:
            if max_group is None:
                data = hdulist[extension][:, min_group:, :, :]
            else:
                data = hdulist[extension][:, min_group: max_group, :, :]
        else:
            if max_group is None:
                data = hdulist[extension][integration, min_group:, :, :]
            else:
                data = hdulist[extension][integration, min_group: max_group, :, :]
    return data

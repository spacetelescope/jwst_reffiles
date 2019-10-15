#! /usr/bin/env python

"""
This is a module containing general tools for dealing with cosmic rays
in JWST data

Author
------

Bryan Hilbert

"""

from astropy.io import fits
from jwst.datamodels import dqflags


def extract_cr_mask(filename):
    """Read in ``filename``, which is assumed to be a 4D ramp, extract
    the DQ bits associated with JUMPS, and return a map containing only
    those bits.
    """
    dq = get_groupdq(filename)
    cr_hits = flag_map(dq, 'JUMP_DET')
    return cr_hits


def flag_map(dq_array, mnemonic):
    """Return a map of cosmic ray flags given a data quality array

    Parameters
    ----------
    dq_array : numpy.ndarray
        GROUP_DQ extension of an integration

    mnemonic : str
        Name of a type of bad pixel. This must be recognized by the
        calibration pipeline (i.e. it must come from jwst.datamodels.dqflags)

    Returns
    -------
    cr_map : numpy.ndarray
        GROUP_DQ extension with all flags other than jumps removed
    """
    cr_map = (dq_array & dqflags.group[mnemonic.upper()] > 0)
    return cr_map


def get_groupdq(filename, refpix):
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
        raise ValueError("ERROR: expecting GROUPDQ to be a 4D array")

    nint, ngroup, ydim, xdim = groupdq.shape
    left, right, bottom, top = refpix
    return groupdq[:, :, bottom: ydim-top, left: xdim-right]

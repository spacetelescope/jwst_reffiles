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
import numpy as np

from jwst_reffiles.utils.file_utils import read_ramp


def collapse_cr_map(dq_map):
    """Transform a 4D array containing cosmic ray hit locations
    (1 for CR hit, 0 for no hit), into a 3D (integration, y, x)
    map that lists for each pixel the group number of the first
    CR hit. If that pixel has no CR hits in the integration, then
    it will have a value of NaN.

    Parameters
    ----------
    dq_map : numpy.ndarray
        4D array of CR hit locations

    Returns
    -------
    index_map : numpy.ndarray
        3D array of group numbers of the first CR hit in each pixel
    """
    nints, ngroups, ny, nx = dq_map.shape

    # Create an array containing all group indexes
    all_groups = np.zeros((1, ngroups, 1, 1))
    all_groups[0, :, 0, 0] = np.arange(ngroups)
    intermediate1 = np.repeat(all_groups, nints, axis=0)
    intermediate2 = np.repeat(intermediate1, ny, axis=2)
    all_indexes = np.repeat(intermediate2, nx, axis=3)

    # Array to contain only group numbers of CR hits
    hit_indexes = np.zeros_like(all_indexes) + np.nan

    # Find the CR flag locations
    hits = np.where(dq_map != 0)

    # All elements are NaN except the groups with CR hits
    hit_indexes[hits] = all_indexes[hits]

    # Find the minimum group number of the CR-hit groups for each pixel
    index_map = np.nanmin(hit_indexes, axis=1)
    return index_map


def create_refpix_map(filename):
    """Create a map of reference pixel locations in a given file.
    Reference pixels have a value of 1. Other pixels have a value
    of zero.

    Parameters
    ----------
    filename : str
        Name of FITS file containing the data quality information

    Returns
    -------
    refmap : numpy.ndarray
        2D array containing a map of reference pixels. 1 indicates a
        reference pixel, while 0 indicates a science pixel
    """
    groupdq = read_ramp(filename, integ_number=0, min_group=0, max_group=1, extension='GROUPDQ')

    # Keep only the REFERENCE_PIXEL flags
    refmap = flag_map(groupdq, 'REFERENCE_PIXEL')
    return refmap


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
    cr_map = (dq_array & dqflags.pixel[mnemonic.upper()] > 0)
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


def signals_per_group(number_of_good, ngroup):
    """Translate an array that lists the number of good groups per
    integration with good signal into an array that shows the number
    of good signal values per group (i.e. look across integrations)

    Parameters
    ----------
    number_of_good : numpy.ndarray
        3D array (integration, y, x) of the number of good groups per
        integration.

    ngroup : int
        The total number of groups in an integration

    Returns
    -------
    good_per_group : numpy.ndarray
        3D array (group, y, x) giving the number of good signal values
        within each group

    bad_per_group : numpy.ndarray
        3D array (group, y, x) giving the number of bad signal values
        within each group
    """
    nint, ny, nx = number_of_good.shape
    good_per_group = np.zeros((ngroup, ny, nx)).astype(np.int)
    bad_per_group = np.zeros((ngroup, ny, nx)).astype(np.int)

    for group in range(ngroup):
        grp_map = np.sum(group < number_of_good, axis=0)
        good_per_group[group, :, :] = grp_map
        bad_per_group[group, :, :] = nint - grp_map
    return good_per_group, bad_per_group

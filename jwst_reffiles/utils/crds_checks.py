#! /usr/bin/env python

"""This module contains general functions that can be used to make basic
CRDS-formatting related checks to various metadata

Author
------

Bryan Hilbert

"""

def validate_pedigree(value):
    """Make sure that the pedigree entry is one of the valid options

    Parameters
    ----------
    value : str
        Pedigree value
    """
    allowed_values = ['GROUND', 'INFLIGHT']

    # For INFLIGHT values, which can have dates tacked onto them,
    # strip off the dates and just check the  text
    if len(value) > 8:
        value, start_date, end_date = value.split(' ')
        if start_date[2] != '/' or start_date[5] != '/' or end_date[2] != '/' or end_date[5] != '/':
            raise ValueError(("ERROR, starting and ending dates in pedigree keyword must have format: mm/dd/yyyy"))

    if value not in allowed_values:
        raise ValueError(("ERROR: pedigree value {} is not valid. Must be one of: {}"
                          .format(value, allowed_values)))

def validate_useafter(value):
    """Make sure that the use_after entry is the correct format

    Parameters
    ----------
    value : str
        Use after value
    """
    # 1964-10-04T00:00:00
    if len(value) != 19:
        if len(value) == 10:
            print('No time given in use_after. Attaching 00:00:00 for the time and attempting to move on.')
            value = value + 'T00:00:00'
        else:
            raise ValueError(("ERROR: useafter entry is not the correct length. Should have 19 characters."))

    try:
        date, time = value.split('T')
    except ValueError:
        raise ValueError('ERROR: use_after format should be yyyy-mm-ddThh:mm:ss')

    try:
        year, month, day = date.split('-')
    except ValueError:
        raise ValueError("ERROR: date format in use_after should be yyyy-mm-dd")

    if month > 12:
        raise ValueError("ERROR: date format in use_after should be yyyy-mm-dd")

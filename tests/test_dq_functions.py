#! /usr/bin/env python

"""
This module contains test functions for the DQ functionality in utils/dq_flags.py
"""
from jwst.datamodels import dqflags
import numpy as np

from jwst_reffiles.utils import dq_flags


def test_flag_map():
    """Test the extraction of a single bit into a map
    """
    dq = np.zeros((5, 5)).astype(np.int)

    # Insert bits

    # Jump only
    dq[0, 0] = dqflags.pixel['JUMP_DET']

    # Saturated only
    dq[1, 1] = dqflags.pixel['SATURATED']

    # Jump and saturated
    dq[2, 2] = dqflags.pixel['JUMP_DET'] + dqflags.pixel['SATURATED']

    # Jump and Hot
    dq[3, 3] = dqflags.pixel['JUMP_DET'] + dqflags.pixel['HOT']

    # Hot only
    dq[4, 4] = dqflags.pixel['HOT']

    jump_map = dq_flags.flag_map(dq, 'JUMP_DET')
    hot_map = dq_flags.flag_map(dq, 'HOT')
    sat_map = dq_flags.flag_map(dq, 'SATURATED')

    jump_map_true = np.zeros((5, 5)).astype(np.int)
    jump_map_true[0, 0] = 1
    jump_map_true[2, 2] = 1
    jump_map_true[3, 3] = 1

    assert np.all(jump_map.astype(np.int) == jump_map_true)

    hot_map_true = np.zeros((5, 5)).astype(np.int)
    hot_map_true[3, 3] = 1
    hot_map_true[4, 4] = 1
    assert np.all(hot_map.astype(np.int) == hot_map_true)

    sat_map_true = np.zeros((5, 5)).astype(np.int)
    sat_map_true[1, 1] = 1
    sat_map_true[2, 2] = 1
    assert np.all(sat_map.astype(np.int) == sat_map_true)

#! /usr/bin/env python

"""Tests for bad_pixel_mask.py"""

import numpy as np

from jwst_reffiles.bad_pixel_mask import bad_pixel_mask as bpm


def test_find_open_and_low_qepixels():
    """Make sure open and low QE pixels are found correctly"""

    # Create test array with one open and one low QE pixel
    test = np.ones((10, 10))
    test[4, 4] = 0.4
    test[8, 8] = 0.4
    test[8, 9] = 1.1
    test[8, 7] = 1.1
    test[9, 8] = 1.1
    test[7, 8] = 1.1

    qe, openpix = bpm.find_open_and_low_qepixels(test)

    compare_low_qe = np.zeros((10, 10))
    compare_low_qe[4, 4] = 1
    assert qe == compare_low_qe

    compare_open = np.zeros((10, 10))
    compare_open[7:10, 7:10] = 2
    compare_open[8, 8] = 1
    assert openpix == compare_open


def test_find_miri_dead_pixels():
    """Make sure MIIR's dead pixels are found correctly"""

    test = np.ones((10, 10, 10))
    test[:, 2, 2] = 0.
    test[1:, 3, 3] = 0.
    test[2:, 4, 4] = 0
    miri_dead = bpm.test_find_miri_dead_pixels(test, min_zero_signal_fraction=0.9)

    comparison = np.zeros((10, 10))
    comparison[2, 2] = 1
    comparison[3, 3] = 1
    assert miri_dead == comparison

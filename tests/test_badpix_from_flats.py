#! /usr/bin/env python

"""Tests for bad_pixel_mask.py"""

from astropy.io import fits
import numpy as np
from jwst.datamodels import dqflags

from jwst_reffiles.bad_pixel_mask import badpix_from_flats as bpm


def test_combine_individual_maps():
    """Test the combination of various bad pixel maps"""

    dead_map = np.zeros((10, 10))
    dead_map[2, 3:6] = 1
    dead_map[8, 4:8] = 1

    low_qe_map = np.zeros((10, 10))
    low_qe_map[1, 3:6] = 1
    low_qe_map[3, 3:6] = 1

    open_map = np.zeros((10, 10))
    open_map[6, 2] = 1
    open_map[3, 8] = 1

    adj_map = np.zeros((10, 10))
    adj_map[2:5, 7:10] = 1
    adj_map[3, 8] = 0

    maps = {'DEAD': dead_map, 'LOW_QE': low_qe_map, 'OPEN': open_map, 'ADJ_OPEN': adj_map}
    do_not_use = ['DEAD', 'OPEN', 'ADJ_OPEN']
    final = bpm.combine_individual_maps(maps, do_not_use)

    manual = np.zeros((10, 10))
    manual[2, 3:6] = dqflags.pixel['DEAD'] + dqflags.pixel['DO_NOT_USE']
    manual[8, 4:8] = dqflags.pixel['DEAD'] + dqflags.pixel['DO_NOT_USE']
    manual[1, 3:6] = dqflags.pixel['LOW_QE']
    manual[3, 3:6] = dqflags.pixel['LOW_QE']
    manual[6, 2] = dqflags.pixel['OPEN'] + dqflags.pixel['DO_NOT_USE']
    manual[2:5, 7:10] = dqflags.pixel['ADJ_OPEN'] + dqflags.pixel['DO_NOT_USE']
    manual[3, 8] = dqflags.pixel['OPEN'] + dqflags.pixel['DO_NOT_USE']
    assert np.all(final == manual)


def test_find_open_and_low_qe_pixels():
    """Make sure open and low QE pixels are found correctly"""

    # Create test array with one open and one low QE pixel
    test = np.ones((10, 10))
    test[4, 4] = 0.4
    test[8, 8] = 0.4
    test[8, 9] = 1.1
    test[8, 7] = 1.1
    test[9, 8] = 1.1
    test[7, 8] = 1.1

    qe, openpix, adjopen = bpm.find_open_and_low_qe_pixels(test)

    compare_low_qe = np.zeros((10, 10))
    compare_low_qe[4, 4] = 1

    compare_open = np.zeros((10, 10))
    compare_open[8, 8] = 1

    compare_adjacent = np.zeros((10, 10))
    compare_adjacent[7:10, 7:10] = 1
    compare_adjacent[8, 8] = 0
    assert np.all(openpix == compare_open)
    assert np.all(adjopen == compare_adjacent)
    assert np.all(qe == compare_low_qe)


def test_flag_type_check():
    """Make sure that non-recognized flags are ignored"""
    flag_types = ['DEAD', 'OPEN', 'ADJ_OPEN', 'LOW_QE', 'WEIRD']
    try:
        flags = bpm.flag_type_check(flag_types)
    except ValueError as e:
        assert 'WEIRD' in e.args[0]


def test_dead_pixels_absolute_rate():
    """Make sure dead pixels defined using a normalized rate image and
    an absolute flux threshold value are found correctly"""
    test = np.ones((10, 10))
    test[2, 2] = 0.04
    test[3, 3] = 0.02
    test[4, 4] = 0.

    dead = bpm.dead_pixels_absolute_rate(test, max_dead_signal=0.05)

    comparison = np.zeros((10, 10))
    comparison[2, 2] = 1
    comparison[3, 3] = 1
    comparison[4, 4] = 1
    assert np.all(dead == comparison)


def test_dead_pixels_sigma_rate():
    """Test that dead pixels defined as having a signal N*sigma below the
    mean are found correctly"""
    test = np.random.normal(loc=1., scale=0.1, size=(10, 10))
    test[5, 5] = 0.003
    test[5, 4] = 0.001
    dev_val = np.std(test)
    mean_val = np.mean(test)
    sigma_val = 5

    dead = bpm.dead_pixels_sigma_rate(test, mean_val, dev_val, sigma=sigma_val)

    manual = (test < (mean_val - dev_val * sigma_val)).astype(int)
    assert np.all(dead == manual)


def test_extract_10th_group():
    """Test the extraction of group 10 from input integrations"""
    test = np.zeros((2, 12, 3, 3))
    test[0, 9, :, :] = 1.
    test[1, 9, :, :] = 2.
    hdu0 = fits.PrimaryHDU(test)
    hdu0.header['FILENAME'] = 'test.fits'
    hdu0.header['EXTNAME'] = 'SCI'
    hdu_list = fits.HDUList([hdu0])
    extracted = bpm.extract_10th_group(hdu_list)
    manual = np.zeros((2, 3, 3))
    manual[0, :, :] = 1.
    manual[1, :, :] = 2.
    assert np.all(extracted == manual)

    test = np.zeros((12, 3, 3))
    test[9, :, :] = 3.
    hdu0 = fits.PrimaryHDU(test)
    hdu0.header['FILENAME'] = 'test.fits'
    hdu0.header['EXTNAME'] = 'SCI'
    hdu_list = fits.HDUList([hdu0])
    extracted = bpm.extract_10th_group(hdu_list)
    manual = np.zeros((1, 3, 3))
    manual[0, :, :] = 3.
    assert np.all(extracted == manual)


def test_image_stats():
    """Be sure the sigma-clipped mean of an image is calcualted correctly"""
    sigma_val = 3
    image = np.ones((10, 10))
    image[:, 0:5] = 2
    mn, dev = bpm.image_stats(image, sigma=sigma_val)
    assert mn == 1.5
    assert dev == 0.5


def test_mean_stdev_images():
    """Check that the mean image of a stack is calculated correctly"""
    test = np.zeros(((25, 10, 10)))
    for i in range(25):
        test[i, :, :] = (i % 5)
    test[0, 0, 0] = 150.
    test[4, 5, 5] = 60.
    mean_img, dev_img = bpm.mean_stdev_images(test, sigma=3.)

    manual_mean = np.zeros((10, 10)) + 2.
    manual_mean[0, 0] = 2.08333333
    manual_mean[5, 5] = 1.91666667
    manual_dev = np.zeros((10, 10)) + np.sqrt(2.)
    manual_dev[0, 0] = 1.381927
    manual_dev[5, 5] = 1.381927
    assert np.allclose(mean_img, manual_mean, atol=1e-6, rtol=0)
    assert np.allclose(dev_img, manual_dev, atol=1e-6, rtol=0)


def test_smooth():
    """Check that 2D box smoothing is done correctly"""
    width = 3
    test = np.arange(25).reshape(5, 5)
    smoothed = bpm.smooth(test, box_width=width)
    #truth = np.array([[1.33333333,  2.33333333,  3.,  3.66666667,  2.66666667],
    #                  [3.66666667,  6.,  7.,  8.,  5.66666667],
    #                  [7., 11., 12., 13.,  9.],
    #                  [10.33333333, 16., 17., 18., 12.33333333],
    #                  [8., 12.33333333, 13., 13.66666667, 9.33333333]])
    truth = np.array([[ 8.        ,  6.33333333,  7.        ,  7.66666667,  9.33333333],
                      [ 7.66666667,  6.        ,  7.        ,  8.        ,  9.66666667],
                      [11.        , 11.        , 12.        , 13.        , 13.        ],
                      [14.33333333, 16.        , 17.        , 18.        , 16.33333333],
                      [14.66666667, 16.33333333, 17.        , 17.66666667, 16.        ]])
    assert np.allclose(smoothed, truth, atol=1e-6, rtol=0)

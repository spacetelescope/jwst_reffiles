#! /usr/bin/env python

"""Tests for badpix_from_darks.py"""

from astropy.io import fits
import numpy as np
from jwst.datamodels import dqflags
import pytest

from jwst_reffiles.dark_current import badpix_from_darks as bpd


def test_add_refpix():
    """Test adding refpix around an array
    """
    data = np.ones((10, 10))
    refpix = (2, 3, 4, 5)

    new_array = bpd.add_refpix(data, refpix)
    yd, xd = new_array.shape
    print, xd, yd

    assert yd == 19
    assert xd == 15
    assert np.all(new_array[0:5, 5] == np.array([0, 0, 0, 0, 1]))
    assert np.all(new_array[13:, 5] == np.array([1, 0, 0, 0, 0, 0]))
    assert np.all(new_array[5, 0:4] == np.array([0, 0, 1, 1]))
    assert np.all(new_array[5, 10:] == np.array([1, 1, 0, 0, 0]))


def test_apply_flags():
    """Test the translation of a bad pixel list into flags
    """
    true_value = dqflags.pixel['HOT'] + dqflags.pixel['DO_NOT_USE']

    print(true_value)

    badmap = np.zeros((10, 10), dtype=np.int)
    true_map = np.zeros((10, 10), dtype=np.uint32)
    for i in range(10):
        badmap[i, i] = 1
        true_map[i, i] =  true_value


    print(true_map)


    flag_names = ['HOT', 'DO_NOT_USE']
    pixmap = bpd.apply_flags(badmap, flag_names)


    print(pixmap)


    assert np.all(pixmap == true_map)


def test_check_metadata():
    """Test header metadata check
    """
    comp = fits.PrimaryHDU()
    comp.header['FILENAME'] = 'something.fits'
    comp.header['EFFEXPTM'] = 14.12
    comp.header['TFRAME'] = 2.35
    comp.header['TGROUP'] = 7.06
    comp.header['NFRAMES'] = 2
    comp.header['NGROUPS'] = 2
    comp.header['SUBARRAY'] = 'FULL'

    hdu = fits.PrimaryHDU()
    hdu.header['FILENAME'] = 'something.fits'
    hdu.header['EFFEXPTM'] = 14.12
    hdu.header['TFRAME'] = 2.35
    hdu.header['TGROUP'] = 7.06
    hdu.header['NFRAMES'] = 2
    hdu.header['NGROUPS'] = 2
    hdu.header['SUBARRAY'] = 'SUB640'

    # This should raise an exception in check_metadata
    with pytest.raises(Exception) as e_info:
        bpd.check_metadata(hdu.header, comp.header)

    # This should not raise an exception
    hdu.header['SUBARRAY'] = 'FULL'
    bpd.check_metadata(hdu.header, comp.header)

    # This should also raise an exception
    hdu.header['NFRAMES'] = 4
    with pytest.raises(Exception) as e_info:
        bpd.check_metadata(hdu.header, comp.header)

    # This should also raise an exception
    hdu.header['NFRAMES'] = 2
    hdu.header['TGROUP'] = 7.5
    with pytest.raises(Exception) as e_info:
        bpd.check_metadata(hdu.header, comp.header)


def test_combine_bad_pixel_types():
    """Test the combination of bad pixel maps. Simply bitwise_or
    combination of input maps
    """
    satval = dqflags.pixel['SATURATED']
    sat_map = np.zeros((5, 5), dtype=np.uint32)
    sat_map[0, 0] = satval

    rcval = dqflags.pixel['RC']
    rc_map = np.zeros((5, 5), dtype=np.uint32)
    rc_map[0, 0] = rcval
    rc_map[1, 1] = rcval

    lowqeval = dqflags.pixel['LOW_QE']
    lowqe_map = np.zeros((5, 5), dtype=np.uint32)
    for i in range(3):
        lowqe_map[i, i] = lowqeval

    hotval = dqflags.pixel['HOT']
    hot_map = np.zeros((5, 5), dtype=np.uint32)
    for i in range(4):
        hot_map[i, i] = hotval

    truthmap = np.zeros((5, 5), dtype=np.uint32)
    truthmap[0, 0] = satval + rcval + lowqeval + hotval
    truthmap[1, 1] = rcval + lowqeval + hotval
    truthmap[2, 2] = lowqeval + hotval
    truthmap[3, 3] = hotval

    combmap = bpd.combine_bad_pixel_types(sat_map, rc_map, lowqe_map, hot_map)
    assert np.all(combmap == truthmap)


def test_find_pix_with_many_jumps():
    """
    """
    map_of_jumps = np.zeros((20, 5, 5))

    # Many jumps
    map_of_jumps[:, 0, 0] = 1

    # Many early jumps -> potential RC pixel
    map_of_jumps[0:12, 1, 1] = 1

    # Many early jumps, but not enough to be RC
    map_of_jumps[0:4, 2, 2] = 1
    map_of_jumps[17: 2, 2] = 1

    map_of_jumps[0:7, 3, 3] = 1

    # Few jumps
    map_of_jumps[5, 4, 4] = 1
    map_of_jumps[12, 4, 4] = 1

    high_jump_pix, rc_pix, jump_count = bpd.find_pix_with_many_jumps(map_of_jumps, max_jump_limit=10,
                                                                     jump_ratio_threshold=5,
                                                                     early_cutoff_fraction=0.25)

    truth_jump_pix = np.zeros((5, 5), dtype=bool)
    truth_jump_pix[0, 0] = True
    truth_jump_pix[1, 1] = True
    assert np.all(high_jump_pix == truth_jump_pix)

    truth_rc_pix = np.zeros((5, 5), dtype=bool)
    truth_rc_pix[1, 1] = True
    assert np.all(rc_pix == truth_rc_pix)

    truth_jump_count = np.zeros((5, 5), dtype=bool)
    truth_jump_count = np.sum(map_of_jumps, axis=0)
    assert np.all(jump_count == truth_jump_count)


def test_slopes_not_cr():
    slope_values = np.zeros((5, 5))

    num_jumps = np.zeros((5, 5), dtype=np.int)
    num_jumps[0, 0] = 1
    num_jumps[1, 1] = 5

    truth_clean_slopes = np.zeros((5, 5), dtype=np.float)
    truth_clean_slopes[0, 0] = np.nan
    truth_clean_slopes[1, 1] = np.nan

    truth_counter = np.ones((5, 5), dtype=np.int)
    truth_counter[0, 0] = 0
    truth_counter[1, 1] = 0

    clean_slopes, counter = bpd.slopes_not_cr(slope_values, num_jumps)

    filler_value = 99.
    nans = ~np.isfinite(clean_slopes)
    clean_slopes[nans] = filler_value
    nans = ~np.isfinite(truth_clean_slopes)
    truth_clean_slopes[nans] = filler_value

    assert np.all(clean_slopes == truth_clean_slopes)
    assert np.all(counter == truth_counter)

@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_combine_clean_slopes():
    slope_img = np.ones((5, 5))

    # Pixel (0, 0) has a CR hit and has a slope of nan
    slope_img[0, 0] = np.nan

    # Create a second slope image with (0, 0) and (1, 1) with CR hits
    slope_img11 = np.ones((5, 5))
    slope_img11[0, 0] = np.nan
    slope_img11[1, 1] = np.nan

    # Create a stack of slope images. (1, 1) is CR hit in 2 of the 8 images
    stack = [slope_img11, slope_img * 1.5, slope_img * 2, slope_img * 2.5,
             slope_img11, slope_img * 1.5, slope_img * 2, slope_img * 2.5]

    # Create an image of CR hit indicators
    crs = np.ones((5, 5), dtype=np.int)

    # (0, 0) has a CR in every integration
    crs[0, 0] = 0

    # Make a CR image with a CR hit at (1, 1) as well as (0, 0)
    crs11 = np.ones((5, 5), dtype=np.int)
    crs11[0, 0] = 0
    crs11[1, 1] = 0

    # Make list of CR arrays
    cr_stack = [crs11, crs, crs, crs, crs11, crs, crs, crs]

    mean_slope, stdev_slope, num_good = bpd.combine_clean_slopes(stack, cr_stack)

    truth_mean_slope = np.zeros((5, 5)) + 1.75
    truth_mean_slope[0, 0] = np.nan
    truth_mean_slope[1, 1] = 2.
    truth_stdev_slope = np.zeros((5, 5)) + 0.5590169943749475
    truth_stdev_slope[0, 0] = np.nan
    truth_stdev_slope[1, 1] = 0.408248290463863
    truth_num_good = np.zeros((5, 5), dtype=np.int) + 8
    truth_num_good[0, 0] = 0
    truth_num_good[1, 1] = 6

    filler_value = 99.
    nan = ~np.isfinite(stdev_slope)
    stdev_slope[nan] = filler_value
    nan = ~np.isfinite(truth_stdev_slope)
    truth_stdev_slope[nan] = filler_value

    # Replace NaNs with a filler value to make comparison easier
    filler = 12345.
    mean_slope[np.isnan(mean_slope)] = filler
    truth_mean_slope[np.isnan(truth_mean_slope)] = filler

    assert np.all(np.isclose(mean_slope, truth_mean_slope, atol=1e-6, rtol=0))
    assert np.all(np.isclose(stdev_slope, truth_stdev_slope, atol=1e-6, rtol=0))
    assert np.all(num_good == truth_num_good)


def test_pedestal_stats():
    pedestal = np.reshape(np.arange(100), (10, 10))

    bad = bpd.pedestal_stats(pedestal, threshold=1)

    truth = np.zeros((10, 10), dtype=bool)
    truth[0:2, :] = True
    truth[2, 0] = True
    truth[8:, :] = True
    truth[7, 9] = True
    assert np.all(bad == truth)


def test_saturated_in_all_groups():
    pedestal = np.zeros((5, 5)) + 12.2
    pedestal[1, 1] = 0.
    pedestal[2, 2] = 0.

    group_sat = np.zeros((5, 5))
    group_sat[2, 2] = 2

    satpix = bpd.saturated_in_all_groups(pedestal, group_sat)

    truth_sat = np.zeros((5, 5), dtype=np.int)
    truth_sat[2, 2] = 1

    assert np.all(satpix == truth_sat)

.. _official_bad_pixel_mask:

Official algorithms for creating the bad pixel mask
---------------------------------------------------

The definition of the `official algorithm <https://outerspace.stsci.edu/display/JWSTCC/Algorithm+details%3A+DQ+Init>`_ for the creation of the bad pixel mask reference file currently describes methods for finding DEAD, LOW QE, OPEN, and ADJACENT TO OPEN pixels.

Using this definition, the *bad_pixel_mask.py* module has been written and added to the JWST_reffiles repository. This module is currently undergoing testing before being fully integrated into the JWST_reffiles framework.

If you encounter any problems or have any questions about the code or its use, feel free to open an issue on the `jwst_reffiles github page <https://github.com/spacetelescope/jwst_reffiles/issues>`_.

There are two separate modules that are called when running *bad_pixel_mask.py*. These are *badpix_from_flats.py* and *badpix_from_darks.py*. The former takes a collection of internal flat field exposures in order to search for DEAD, LOW QE, OPEN, and ADJACENT TO OPEN pixels. The latter takes a collection of dark current exposures in order to search for NOISY, HOT, RC, TELEGRAPH, and LOW_PEDESTAL pixels. Both of these modules expect input data in at least 2 calibration states, as detailed in the table below.

+------------------+-----------------------------------+-----------------------------+
|      Module      |              Calibration States   |   Typical filename suffixes |
+------------------+-----------------------------------+-----------------------------+
|badpix_from_flats |  slope images, uncalibrated ramps |  rate, uncal                |
+------------------+-----------------------------------+-----------------------------+
|badpix_from_darks |  slope images, "fitopt" files,    |  rate, fitopt, jump, uncal  |
|                  |  CR-flagged ramps,                |                             |
|                  |  (uncalibrated ramps: MIRI-only)  |                             |
+------------------+-----------------------------------+-----------------------------+


When running *bad_pixel_mask.py* as a standalone package, the input file lists described in the table above must be provided manually. In the future, when running via *mkrefs.py*, raw or partially-calibrated files will be allowed as inputs. *mkrefs.py* will then call the calibration pipeline and produce the needed files.

To use the bad pixel mask generator function as a standalone package, use the code shown below. Keywords in the call are linked to detailed descriptions below the code snippet. Note that all of the keywords in the call below have defaults defined in the code (and in fact the call below is using all default values), so you do not have to specify any keywords where you wish to use the default.

.. parsed-literal::

    from glob import glob
    from jwst_reffiles.bad_pixel_mask.bad_pixel_mask import bad_pixels

    flat_rate_files = sorted(glob('/path/to/flats/*_rate.fits'))
    flat_fluxcheck_files = sorted(glob('/path/to/flats/*uncal.fits'))

    dark_rate_files = sorted(glob('/path/to/darks/*_rate.fits'))
    dark_jump_files = sorted(glob('/path/to/darks/*_jump.fits'))
    dark_fitopt_files = sorted(glob('/path/to/darks/*_fitopt.fits'))
    dark_uncal_files = sorted(glob('/path/to/darks/*_uncal.fits'))

    bad_pixels(flat_slope_files=flat_rate_files,
               dead_search_ =True,
               low_qe_and_open_search_ =True,
               dead_search_type_ ='sigma_rate',
               flat_mean_sigma_threshold_ =3,
               flat_mean_normalization_method_ ='smoothed',
               smoothing_box_width_ =15,
               smoothing_type_ ='Box2D',
               dead_sigma_threshold_ =5.,
               max_dead_norm_signal_ =None,
               run_dead_flux_check_ =False,
               dead_flux_check_files_ =None,
               flux_check_ =45000,
               max_low_qe_norm_signal_ =0.5,
               max_open_adj_norm_signal_ =1.05,
               manual_flag_file_ ='default',
               flat_do_not_use_ =[],
               dark_slope_files=dark_files,
               dark_uncal_files=None,
               dark_jump_files=dark_jump_files,
               dark_fitopt_files=dark_fitopt_files,
               dark_stdev_clipping_sigma_ =5.,
               dark_max_clipping_iters_ =5,
               dark_noisy_threshold_ =5,
               max_saturated_fraction_ =0.5,
               max_jump_limit_ =10,
               jump_ratio_threshold_ =5,
               early_cutoff_fraction_ =0.25,
               pedestal_sigma_threshold_ =5,
               rc_fraction_threshold_ =0.8,
               low_pedestal_fraction_ =0.8,
               high_cr_fraction_ =0.8,
               flag_values_ ={'hot': ['HOT'], 'rc': ['RC'], 'low_pedestal': ['OTHER_BAD_PIXEL'], 'high_cr': ["TELEGRAPH"]},
               dark_do_not_use_ =['hot', 'rc', 'low_pedestal', 'high_cr'],
               plot_ =False,
               output_file_ ='./test_bpm.fits',
               author_ ='jwst_reffiles',
               description_ ='A bad pix mask',
               pedigree_ ='GROUND',
               useafter_ ='2019-04-01 00:00:00',
               history_ ='',
               quality_check_ =False)


.. _dead_search:

Dead Search
-----------

Boolean. If True, a dead pixel search (of type dead_search_type_) is performed.

.. _low_qe_and_open_search:

Low QE and Open Search
----------------------

Boolean. If True, a search for Low QE, Open, and Adjacent to Open pixels is performed.

.. _dead_search_type:

Dead Search Type
----------------

Use this string parameter to specify which type of dead pixel search to perform. Options are:

.. parsed-literal::

    'sigma_rate': Using a normalized signal rate image, dead pixels
                  are defined as those with a rate smaller than
                  dead_sigma_threshold_ standard deviations below
                  the mean.
    'absolute_rate': Using a normalized signal rate image, dead pixels
                     are defined as those with a rate less than
                     max_dead_norm_signal_.

.. _flat_mean_sigma_threshold:

Flat Mean Sigma Threshold
-------------------------

Number of standard deviations to use when sigma-clipping to calculate the mean slope image or the mean across the detector
when working with flat field images.

.. _flat_mean_normalization_method:

Flat Mean Normalization Method
------------------------------

Specify how the mean flat field image is normalized prior to searching for bad pixels. Options are:

.. parsed-literal::

    'smoothed': Mean image will be smoothed using a smoothing_box_width_ x smoothing_box_width_ box kernel. The mean
    image is then normalized by this smoothed image.

    'none': No normalization is done. Mean slope image is used as is

    'mean': Mean image is normalized by its sigma-clipped mean

.. _smoothing_box_width:

Smoothing Box Width
-------------------

Width in pixels of the 2D box kernel to use to compute the smoothed mean flat field image

.. _smoothing_type:

Smoothing Type
--------------

Type of smoothing to do when creating a smoothed mean flat field image. ``Box2D`` or ``Median``. Box2D uses an astropy Box2DKernel, while
Median uses a scipy median_filter.


.. _dead_sigma_threshold:

Dead Sigma Threshold
--------------------

Number of standard deviations below the mean at which a pixel is considered dead when using the ``sigma_rate`` :ref:`dead search type<dead_search_type>`.

.. _max_dead_norm_signal:

Maximum Dead Normalized Signal
------------------------------

Maximum normalized signal rate of a pixel that is considered dead when using the ``absolute_rate`` :ref:`dead search type<dead_search_type>`.


.. _run_dead_flux_check:

Run Dead Flux Check
-------------------

Boolean controlling whether or not to search flagged dead pixels for flux. This search potentially removes false positives, as pixels that are saturated in all groups of an integration will have a value of zero in the slope image and therefore appear dead.


.. _dead_flux_check_files:

Dead Flux Check Files
---------------------

Files to use for the dead flux check. These should be raw (i.e. unal) files.

.. _flux_check:

Flux Check
----------

Signal level threshold to use during the dead flux check test. Pixels flagged as dead and with signals less than this signal level are considered dead.

.. _max_low_qe_norm_signal:

Maximum Low QE Normalized Signal
--------------------------------

The maximum normalized signal a pixel can have and be considered low QE.

.. _max_open_adj_norm_signal:

Maximum Normalized Signal in Adjacent to Open Pixels
----------------------------------------------------

The maximum normalized signal a pixel adjacent to a low QE pixel can have in order for the low QE pixel to be reclassified as OPEN

.. _manual_flag_file:

Manual Flag File
----------------

Ascii file containing a list of pixel coordinates and bad pixel types to be added to those found in badpix_from_flats.py and placed in the output bad pixel file. If left as 'default', the bad pixel file in the jwst_reffiles repository will be used.

.. _flat_do_not_use:

Flat Do Not Use
---------------

List of bad pixel types (from the flat field files) where the DO_NOT_USE flag should also be applied (e.g. ['DEAD', 'LOW_QE', 'OPEN', 'ADJ_OPEN'])

.. _dark_stdev_clipping_sigma:

Dark Stdev Clipping Sigma
-------------------------

Number of sigma to use when sigma-clipping the 2D array of standard deviation values from the dark current slope files. The sigma-clipped mean and standard deviation are used to locate noisy pixels.

.. _dark_max_clipping_iters:

Dark Max Clipping Iterations
----------------------------

Maximum number of iterations to use when sigma clipping to find the mean and standard deviation values that are used when locating noisy pixels.


.. _dark_noisy_threshold:

Dark Noisy Threshold
--------------------

Number of sigma above the mean noise (associated with the slope) to use as a threshold for identifying noisy pixels in the dark current data.


.. _max_saturated_fraction:

Maximum Saturated Fraction
--------------------------

When identifying pixels that are fully saturated (in all groups of an integration), this is the fraction of integrations within which a pixel must be fully saturated before flagging it as HOT.


.. _max_jump_limit:

Maximum Jump Limit
------------------

The maximum number of jumps a pixel can have in an integration before it is flagged as a ``high jump`` pixel (which may be flagged as noisy later).

.. _jump_ratio_threshold:

Jump Ratio Threshold
--------------------

Cutoff for the ratio of jumps early in the ramp to jumps later in the ramp. Pixels with a ratio greater than this value (and which also have a high total number of jumps) will be flagged as potential (I)RC pixels.


.. _early_cutoff_fraction:

Early Cutoff Fraction
---------------------

Fraction of the integration to use when comparing the jump rate early in the integration to that across the entire integration. Must be <= 0.5

.. _pedestal_sigma_threshold:

Pedestal Sigma Threshold
------------------------

Used when searching for RC pixels via the pedestal image. Pixels with pedestal values more than ``pedestal_sigma_threshold`` above the mean are flagged as potential RC pixels

.. _rc_fraction_threshold:

RC Fraction Threshold
---------------------

Used when searching for RC pixels. This is the fraction of input files within which the pixel must be identified as an RC pixel before it will be flagged as a permanent RC pixel

.. _low_pedestal_fraction:

Low Pedestal Fraction
---------------------

This is the fraction of input files within which a pixel must be identified as a low pedestal pixel before it will be flagged as a permanent low pedestal pixel


.. _high_cr_fraction:

High CR Fraction
----------------

This is the fraction of input files within which a pixel must be flagged as having a high number of jumps before it will be flagged as permanently noisy


.. _flag_values:

Flag Values
-----------

This dictionary maps the types of bad pixels searched for to the flag mnemonics to use when creating the bad pixel file. Keys are the types of bad pixels searched for, and values are lists that include mnemonics recognized by the jwst calibration pipeline.

e.g. {'hot': ['HOT'], 'rc': ['RC'], 'low_pedestal': ['OTHER_BAD_PIXEL'], 'high_cr': ["TELEGRAPH"]}


.. _dark_do_not_use:

Dark Do Not Use
---------------

List of bad pixel types from the dark current data to be flagged as DO_NOT_USE.

e.g. ['hot', 'rc', 'low_pedestal', 'high_cr']


.. _plot:

Plot
----

If True, produce plots of intermediate results.


.. _output_file:

Output File
-----------

Name of the CRDS-formatted bad pixel reference file to save the final bad pixel map into

.. _author:

Author
------

CRDS-required name of the reference file author, to be placed in the referece file header

.. _description:

Description
-----------

CRDS-required description of the reference file, to be placed in the reference file header

.. _pedigree:

Pedigree
--------

CRDS-required pedigree of the data used to create the reference file

.. _useafter:

Useafter
--------

CRDS-required date of earliest data with which this referece file should be used. (e.g. '2019-04-01 00:00:00')

.. _history:

History
-------

String containing any text you wish to place in the HISTORY keyword of the output bad pixel mask reference file. Note that all input filenames will automatically be placed in the HISTORY keyword independent of the string entered here.

.. _quality_check:

Quality Check
-------------

Boolean. If True, the pipeline is run using the output reference file to be sure the pipeline doens't crash


.. _limitations:

Current Limitations
-------------------

Currently, only one type of dead pixel search can be performed for a given call.

.. _official_bad_pixel_mask:

Official algorithms for creating the bad pixel mask
---------------------------------------------------

The definition of the `official algorithm <https://outerspace.stsci.edu/display/JWSTCC/Algorithm+details%3A+DQ+Init>`_ for the creation of the bad pixel mask reference file currently describes methods for finding DEAD, LOW QE, OPEN, and ADJACENT TO OPEN pixels.

Using this definition, the *bad_pixel_mask.py* module has been written and added to the JWST_reffiles repository. This module is currently undergoing testing before being fully integrated into the JWST_reffiles framework.

If you encounter any problems or have any questions about the code or its use, feel free to open an issue on the `jwst_reffiles github page <https://github.com/spacetelescope/jwst_reffiles/issues>`_.

When running *bad_pixel_mask.py* as a stand-alone package, the input files must be count rate (slope) images. In the future, when running via *mkrefs.py*, raw or partially-processed pipeline outputs will be allowed as inputs.

To use the bad pixel mask generator function as a standalone package, use the code shown below. Keywords in the call are linked to detailed descriptions below the code snippet. Note that all of the keywords in the call below have defaults defined in the code (and in fact the call below is using all default values), so you do not have to specify any keywords where you wish to use the default.

.. parsed-literal::

    from glob import glob
    from jwst_reffiles.bad_pixel_mask import bad_pixel_mask as bpm

    list_of_input_files = glob('/my/calibration/inputs/*_rate.fits')
    bpm.find_bad_pix(list_of_input_files, dead_search_ =True, low_qe_and_open_search_ =True,
                     dead_search_type_ ='sigma_rate', sigma_threshold_ =3,
                     smoothing_box_width_ =15, dead_sigma_threshold_ =5.,
                     max_dead_norm_signal_ =0.05, max_low_qe_norm_signal_ =0.5,
                     max_open_adj_norm_signal_ =1.05,
                     dead_zero_signal_fraction_ = 0.9,
                     do_not_use_ =['DEAD', 'LOW_QE', 'OPEN', 'ADJ_OPEN'],
                     output_file_ =outfile_sig,
                     author_ ='jwst_reffiles', description_ ='A bad pix mask',
                     pedigree_ ='GROUND', useafter_ ='2019-04-01 00:00:00',
                     history_ ='doomed to be repeated', quality_check_ =False)

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
    'zero_signal': Using a stack of integrations (rather than slope images),
                   one group is extracted from each, and a pixel is flagged
                   as dead if its signal is zero in at least
                   (dead_zero_signal_fraction_ * 100) percent of the
                   extracted groups.

.. _sigma_threshold:

Sigma Threshold
---------------

Number of standard deviations to use when sigma-clipping to calculate the mean slope image or the mean across the detector


.. _smoothing_box_width:

Smoothing Box Width
-------------------

Width in pixels of the 2D box kernel to use to compute the smoothed mean image

.. _dead_sigma_threshold:

Dead Sigma Threshold
--------------------

Number of standard deviations below the mean at which a pixel is considered dead when using the ``sigma_rate`` :ref:`dead search type<dead_search_type>`.

.. _max_dead_norm_signal:

Maximum Dead Normalized Signal
------------------------------

Maximum normalized signal rate of a pixel that is considered dead when using the ``absolute_rate`` :ref:`dead search type<dead_search_type>`.

.. _dead_zero_signal_fraction:

Dead Zero Signal Fraction
----------------------------

For the ``zero_signal`` case of :ref:`dead search type<dead_search_type>`, where dead pixels are defined as having zero signal, this is the fration of input integrations in which a pixel must have zero signal for it to be flagged as dead. (i.e. 0.9 means a pixel must have no signal in 90% of the input integrations for it to be flagged as dead.)

.. _max_low_qe_norm_signal:

Maximum Low QE Normalized Signal
--------------------------------

The maximum normalized signal a pixel can have and be considered low QE.

.. _max_open_adj_norm_signal:

Maximum Normalized Signal in Adjacent to Open Pixels
----------------------------------------------------

The maximum normalized signal a pixel adjacent to a low QE pixel can have in order for the low QE pixel to be reclassified as OPEN


.. _do_not_use:

Do Not Use
----------

List of bad pixel types where the DO_NOT_USE flag should also be applied (e.g. ['DEAD', 'LOW_QE', 'OPEN', 'ADJ_OPEN'])

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

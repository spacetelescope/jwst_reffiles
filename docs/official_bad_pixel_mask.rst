.. _official_bad_pixel_mask:

Official algorithms for creating the bad pixel mask
---------------------------------------------------

The definition of the `official algorithm <https://outerspace.stsci.edu/display/JWSTCC/Algorithm+details%3A+DQ+Init>`_ for the creation of the bad pixel mask reference file currently describes methods for finding DEAD, LOW QE, OPEN, and ADJACENT TO OPEN pixels.

Using this definition, the *bad_pixel_mask.py* module has been written and added to the JWST_reffiles repository. This module is currently undergoing testing before being fully integrated into the JWST_reffiles framework.

To use the bad pixel mask generator function as a standalone package:

::

    from glob import glob
    from jwst_reffiles.bad_pixel_mask import bad_pixel_mask as bpm

    list_of_input_files = glob('file_number_*.fits')
    mark_as_do_not_use = ['DEAD', 'LOW_QE', 'OPEN', 'ADJ_OPEN']
    bpm.find_bad_pix(list_of_input_files, dead_search_type='sigma_rate', sigma_threshold=3,
                     do_not_use=mark_as_do_not_use)


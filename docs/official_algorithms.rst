.. _official_algorithms:

Official algorithms for creating refrence files
-----------------------------------------------

The `JWST Calibration Reference File Generation Tools <https://outerspace.stsci.edu/display/JWSTCC/JWST+Calibration+Reference+File+Generation+Tools>`_ group is tasked with defining a set of officially endorsed algorithms for creating JWST reference files. The goals of the group, and the defined algorithms, is to reduce duplication of effort and identify algorithms that may be used across instruments. Once defined, the Referene File Implementation subgroup will implement these algorithms in python modules which will be added to the JWST_reffiles repository so that they are available for all instrument teams to use.

The current list of complete or in-progress officially-defined algorithms include:

.. _bad_pixel_mask:

:ref:`Bad Pixel Mask (used in DQ_Init pipeline step) <official_bad_pixel_mask>`

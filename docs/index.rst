.. FALCO documentation master file, created by
   sphinx-quickstart on Thu Jun 21 20:18:02 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

JWST_REFFILES: Reference File Generator Tool
============================================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   install.rst
   organize_input_files.rst
   individual_reference_file_creation.rst
   final_reference_file_creation.rst
   plug_in_user_module.rst
   configuration_files.rst
   quickstart.ipynb
   api.rst

The **jwst_reffiles** package contains a framework to organize and automate the creation of calibration reference files for JWST instruments. These reference files are used by the `JWST calibration pipeline <https://jwst-pipeline.readthedocs.io/en/stable/>`_ for data calibration.

There are three general stages within **jwst_reffiles** for reference file creation.

1. :ref:`Identification and orgnization of input files <organize_inputs>`
2. :ref:`Creation of an individual reference file <individual_reffile>`
3. :ref:`Combining individual reference files into a final reference file <final_reffile>`

In this documentation, we use the term *individual reference file* to refer to a reference file created from a single input file or group of files. For example, if a single dark current exposure can be used to create a readnoise image, then this readnoise image, onced saved in the appropriate format, is an *individual reference file*. Conversely, we refer to the mean readnoise image created from a stack of individual readnoise images as the *final reference file*.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

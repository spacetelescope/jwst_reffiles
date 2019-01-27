.. FALCO documentation master file, created by
   sphinx-quickstart on Thu Jun 21 20:18:02 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

JWST_REFFILES: Reference File Generator Tool
============================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install.rst
   plug_in_user_module.rst
   configuration_files.rst
   quickstart.ipynb
   api.rst

The **jwst_reffiles** package contains a framework to organize and automate the creation of calibration reference files for JWST instruments. These reference files are used by the `JWST calibration pipeline <https://pipeline-link.here>`_ for data calibration.

There are three general stages within **jwst_reffiles** for reference file creation.

1. :ref:`Identification and orgnization of input files <organize_inputs>`
2. :ref:`Creation of an individual reference file <individual_reffile>`
3. :ref:`Combining individual reference files into a final reference file <final_reffile>`

In this documentation, we use the term *individual reference file* to refer to a reference file created from a single input file or group of files. For example, if a single dark current exposure can be used to create a readnoise image, then this readnoise image, onced saved in the appropriate format, is an *individual reference file*. Conversely, we refer to the mean readnoise image created from a stack of individual readnoise images as the *final reference file*.


**move the steps below to their own separate pages**


.. _organize_inputs:

Identification and organization of input files
----------------------------------------------

The primary function performed by **jwst_reffiles** is to perform the organization and bookkeeping of input files to be used in the creation of reference files. The package is designed such that users can provide lists of input files (e.g. dark current files, flat field files), along with information on which types files should be used to create a given type of reference file, as well as some simple rules on which input files can be used in combination with others. **jwst_reffiles** then performs several tasks to prepare these input files for use in creating reference files.

Group input files
+++++++++++++++++

Within a :ref:`configuration file <main_cfg>`, the user can specify which combinations of inputs are needed to create a given type of reference file. For example, the user can specify that an individual gain reference file should be created using a pair of flat field files and a pair of dark current files.

In addition, the user can define some simple rules for creating these groups of input files. Currently this includes a maximum time between files in the group. For example, to minimize any systematic variation in dark current, the user may specify that all dark current files in a single input group must have been observed within a 2 day time period. **Are there any other rules?**

Using the grouping definitions and input file rules, **jwst_reffiles** will then create groups of input files, where each group can be used to create an :ref:`individual reference file <individual_reffile>`.


Determine which pipeline calibration steps must be run
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Also within the main :ref:`configuration file <main_cfg>`, the user can specify which calibration pipeline steps are necessary to be complete on the input files prior to creating individual reference files. **jwst_reffiles** will then examine the input files, determine which pipeline steps have already been completed, and which steps remain to be done. It then generates, for each file, the appropriate call to the JWST calibration pipeline such that the output will have all necessary steps completed. **jwst_reffiles** will run the pipeline commands, save the outputs, and use these in the input file groups.



.. _individual_reffile:

Create individal reference files
--------------------------------

Once the input files have been organized and calibrated to the requested level, **jwst_reffiles** will create an individual reference file from each group of input files. The software package used to create each type of reference file is specified by the user in the :ref:`main configuration file <main_cfg>`. This package must be plugged in to **jwst_reffiles** in the manner described in the :ref:`plug-in <plugin>` page.





.. _final_reffile:

Create final reference files
----------------------------




While **jwst_reffiles** provides the framework, user-generated packages can be plugged in






Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

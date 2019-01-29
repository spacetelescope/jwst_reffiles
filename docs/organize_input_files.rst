.. _organize_inputs:

Identification and Orgnization of Input Files
=============================================

The primary function performed by **jwst_reffiles** is to perform the organization and bookkeeping of input files to be used in the creation of reference files. The package is designed such that users can provide lists of input files (e.g. dark current files, flat field files), along with information on which types of files should be used to create a given type of reference file, as well as some simple rules on which input files can be used in combination with others. **jwst_reffiles** then performs several tasks to prepare these input files for use in creating reference files.

.. _group_inputs:

Group input files
+++++++++++++++++

Within a :ref:`configuration file <main_cfg>`, the user can specify which combinations of inputs are needed to create a given type of reference file. For example, the user can specify that an individual gain reference file should be created using a pair of flat field files and a pair of dark current files.

In addition, the user can define some simple rules for creating these groups of input files. Currently this includes a maximum time between files in the group. For example, to minimize any systematic variation in dark current, the user may specify that all dark current files in a single input group must have been observed within a 2 day time period. **Are there any other rules?**

Using the grouping definitions and input file rules, **jwst_reffiles** will then create groups of input files, where each group can be used to create an :ref:`individual reference file <individual_reffile>`.

.. _define_pipeline_steps:

Determine which pipeline calibration steps must be run
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Also within the main :ref:`configuration file <main_cfg>`, the user can specify which calibration pipeline steps are necessary to be complete on the input files prior to creating individual reference files. **jwst_reffiles** will then examine the input files, determine which pipeline steps have already been completed, and which steps remain to be done. It then generates, for each file, the appropriate call to the JWST calibration pipeline such that the output will have all necessary steps completed. **jwst_reffiles** will run the pipeline commands, save the outputs, and use these in the input file groups.

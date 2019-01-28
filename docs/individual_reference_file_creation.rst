.. _individual_reffile:

Create individal reference files
--------------------------------

Once the input files have been :ref:`organized and calibrated <organize_inputs>` to the requested level, **jwst_reffiles** will create an individual reference file from each group of input files. The software package used to create each type of reference file is specified by the user in the :ref:`main configuration file <main_cfg>`, towards the bottom. An example, called :ref:`example_readnoise_module <example_readnoise_module>`, is shown in the linked configuration file. This package must be plugged in to **jwst_reffiles** in the manner described in the :ref:`plug-in <plugin>` page.

Values for package-specific variables are provided by the user in this section of the configuration file as well. There are several inputs that are required for all reference file generation scripts. These are the :ref:`reftype <reftype>`, :ref:`imtypes <imtypes>`, and :ref:`ssbsteps <ssbsteps>` entries. In addition to these, any variables that are unique to the package can be listed here. In the :ref:`example_readnoise_module <example_readnoise_module>` example, these include *verbose*, *boxsize*, and *sigma_threshold*. As described in the :ref:`Plug-in <plugin>` page, values entered here will be passed to the *example_readnoise_module* when it is called, unless they are overridden on the command line.

For each reference file generation module specified in the configuration file's :ref:`reflabels <reflabels>` section, the appropriate software will be called once for each of the :ref:`previously created <group_inputs>` group of pipeline-processed input files. Output files will be placed in the directory specified by the :ref:`reference file output directory <basename>` in the configuration file.

Note that each reference file creation module must take care of creating and saving the CRDS-`formatted <http://jwst-reffiles.stsci.edu/index.html>`_ reference file itself. This means that the plugged-in module should use the appropriate `jwst` `data model <https://jwst-pipeline.readthedocs.io/en/stable/jwst/datamodels/index.html#module-jwst.datamodels>`_, correctly populate all data and metadata, and save the file in the location provided by *jwst_reffiles*.


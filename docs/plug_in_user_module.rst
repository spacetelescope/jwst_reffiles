How to Plug In A User-Generated Module
======================================

Jwst_reffiles is simply a framework upon which reference file generation modules can be built. In order to add functionality for creating and averaging/combining reference files, users will need to plug-in their code to the jwst_reffiles framework. This page describes how to plug an existing module into jwst_reffiles. An example is provided in *jwst_reffiles/example*, where we have plugged in *example_readnoise_module.py* by creating *mkref_example_readnoise_module.py* following the steps below.

Assume you have a python script that creates a gain reference file. It takes as input several fits files containing dark current and flat field exposures, as well as optional arguments, and produces a CRDS-formatted gain reference file which is saved as a fits file. For this example, we will call this script *my_gain_module.py*.

.. _wrapping:

Plug-in Steps
-------------

The steps needed to plug *my_gain_module.py* into jwst_reffiles are:

1. Make a copy of *jwst_reffiles/templates/plugin_template.py* and rename it to the name of your module with `mkref_` prepended. (e.g. *mkref_my_gain_module.py*). This is your plug-in module.
2. Make your plug-in module executable. (e.g. ``chmod ugo+x mkref_my_gain_module.py``). This is necessary so that the plug-in module can be called via the command line without prepending ``python`` to the call. This is helpful for cases where commands are being sent to Condor.
3. Add the appropriate parameters and a call to your script (e.g. *my_gain_module.py*). Details on edits are given in the :ref:`edits list <edits_list>` below.
4. Add the appropriate parameters for your script to a new section towards the bottom (see for example the block for ``example_readnoise_module``) of the main :ref:`config file <main_cfg>`. An example of what to add is shown in the :ref:`Config File Update <config_file_update>` section below.

.. _edits_list:

Plug-in Wrapper Edits
---------------------

Within your plug-in module, make the following edits:

1. Add an import statement for your script. e.g. ``from my_package.gain import my_gain_module``
2. In ``__init__``, set the reflabel and reftype. The reflabel is simply the name of your script. (e.g. *my_gain_module*) The reftype is the jwst_reffiles shorthand for the type of reference file that your script produces. (e.g. gain) The list of allowed reftypes is in ????
3. In the ``extra_optional_arguments`` function, add any parameters for your script that you may want to call via the command line. These parameters will be available during command line calls to **mkrefs.py**. Parameters set via the command line will override the values set in the config file.
4. Within the ``callalgorithm`` function, create a call to your script (e.g. *my_gain_module.py*) When creating the call to your module, you should have access to all necessary information. The names of input fits files for a single run of your script will be in a list in ``self.inputimages['output_name']``. Similarly, the output name of the reference file to be created will be in ``self.args.outputreffilename``. Any additional parameters (as defined in the config file and/or on the command line) will be available in the ``self.parameters`` dictionary.


.. _config_file_update:

Config File Update
------------------

You may add any parameters for your module to the mkrefs config file. This is done in a block under the name of your module.

In this example, let's assume that *mkref_my_gain_module.py* takes two arguments:

1. num_boxes: The detector will be divided up into a *num_boxes*x*num_boxes* grid of squares
2. max_signal: The maximum signal value to use in the calculation of the gain

Add a section to the config file containing the arguments and corresponding values.

::

    my_gain_module:
        num_boxes: 16
        max_signal: 20000

The values for these parameters will then be read in from the config file when mkrefs.py is called. Keep in mind that if you have added these parameters to ``extra_optional_arguments`` function, then you can override the values in the config file via the command line.


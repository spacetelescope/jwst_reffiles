#! /usr/bin/env python

"""This file contains template code that can be used to build a class for
plugging a user-generated reference file creation script into jwst_reffiles.

See mkref_example_readnoise_module.py in jwst_reffiles/example for an example.

This class will wrap around the user-generated code and allow mkrefs
to call the code.

Several steps are needed in order to tailor this file for your reference
file generation script:

::

    1. Begin by making a copy of this file and renaming it to
       mkref_your_module_name.py. If you are wrapping "niriss_bpm.py"
       then rename your copy of this file to mkref_niriss_bpm.py.
    2. MAKE YOUR RENAMED COPY OF THIS FILE EXECUTABLE!!
    3. Edit your copy of this file following the guidelines below


::
    Edits to make:
    1. Add an import statement for your module
    2. In __init__, set the reflabel and reftype
    3. In extra_optional_arguments, add any parameters that you may want to call
       via the command line.
    4. In callalgorithm, call your module. Note that mkrefs.py will supply the
       input filenames, and will also specify the name of the output file that
       your module will create. All parameters will be available in the
       self.parameters dictionary.


Futher details: ``Setting parameters for your module``

All arguments/keyword inputs to the reference-file creation module
must be included in:

::

    1) the config file
    2) extra_optional_arguments below
    3) both

You may add any parameters for your module to the mkrefs config file.
This is done in a block under the name of your module.

For example, if you are wrapping a module called nircam_readnoise.py,
which takes arguments: sigma_threshold and box_size, add a section to the
config file containing the arguments and corresponding values.

:

    nircam_readnoise:
        sigma_threshold: 3
        box_size: 64

The values for these parameters will then be read in from the config
file when mkrefs.py is called. If you wish, you may also add the
options to the argument parser in extra_optional_arguments in this module.
This will allow you to override the values set in the config file when calling
mkrefs.py via the command line.


Further details: ``Calling your module within callalgorithm``

When creating the call to your module, you should have access to all necessary
information. The names of input fits files for a single run of your module
will be in a list in self.inputimages['output_name']. Similarly, the output
name of the reference file to be created will be in self.args.outputreffilename.

Any additional parameters (as defined in the config file and/or on the command
line) will be available in the self.parameters dictionary.
"""
import argparse
import os
import re
import sys
import types

from jwst_reffiles.plugin_wrapper import mkrefclass_template

# import your script!
from my_package.scripts.my_script import MyClass


class mkrefclass(mkrefclass_template):
    def __init__(self, *args, **kwargs):
        mkrefclass_template.__init__(self, *args, **kwargs)

        # self.reflabel should be set to the name of the user-generated
        # script that this class is wrapping around. For example, if your
        # code is in nircam_readnoise.py, then set self.reflabel to
        # 'nircam_readnoise'
        self.reflabel = 'example_readnoise_module_to_plug_in'

        # This must be set to the appropriate reference file label.
        # Allowed labels include:
        # e.g. 'rdnoise', 'bpm', 'gain'
        self.reftype = 'rdnoise'

    def extra_optional_arguments(self, parser):
        """Any arguments added here will give the option of overriding
        the default argument values present in the config file. To override,
        call these arguments from teh command line in the call to mkrefs.py
        """
        parser.add_argument('--parameter1',
                            help='Specifies some important feature in your module')
        return(0)

    def callalgorithm(self):
        """Call your module. The only requirement is that the output
        reference file is saved as self.args.outputreffilename

        mkrefs.py will supply the input files in self.inputimages['output_name'].
        This will be a list containing the filenames to use as input. The
        file types (e.g. dark, flat) associated with each filename are
        contained in self.inputimages['imtype']. From this, you can specify
        the appropriate file names in the call to your module.
        """
        print('*** This is the output filename of the reference file: {}'
              .format(self.args.outputreffilename))
        print('*** This is the file with the list of the input images: {}'
              .format(self.args.imageslist_filename))
        print('*** The table in this file has been loaded into self.inputimages')

        print('*** These are the filenames of the input images, and some more info:')
        for i in range(len(self.inputimages)):
            print('image {} has the type {} and taken at MJD={}'.format(self.inputimages['output_name'][i],
                                                                        self.inputimages['imtype'][i],
                                                                        self.inputimages['MJD'][i]))
        print('*** All parameters are available in self.parameters:')
        for key in self.parameters:
            print('{}: {}'.format(key, self.parameters[key]))

        # Call the wrapped module and provide the proper arguments from the
        # self.parameters dictionary.
        ron = My_Reffile_Module(self.inputimages['output_name'].data[0],
                                param1=self.parameters['parameter1_name'],
                                param2=self.parameters['parameter2_name'],
                                output_name=self.args.outputreffilename)
        return(0)


if __name__ == '__main__':
    mkref = mkrefclass()
    mkref.make_reference_file()

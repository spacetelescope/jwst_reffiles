Installing JWST_REFFILES
========================

Note that *jwst_reffiles*' main dependency, the `JWST calibration pipeline <https://github.com/spacetelescope/jwst>`_ requires Python 3.5 or above.

Environment Creation
--------------------

It may be easiest to create an environment to contain *jwst_reffiles* by first creating an environment around the JWST calibration pipeline. See the installation instructions in the `JWST repo <https://github.com/spacetelescope/jwst>`_

JWST_REFFILES installation
--------------------------

To download and install *jwst_reffiles*, clone from `GitHub <https://github.com/spacetelescope/jwst_reffiles>`_. Make sure you are in your desired environment, and then install using setup.py.

1. git clone https://github.com/spacetelescope/jwst_reffiles.git
2. cd jwst_reffiles
3. python setup.py install


Calling JWST_REFFILES
=====================

Currently, *jwst_reffiles* is designed such that it can only be called from the command line. The top-level script to be called is `mkrefs.py`. Given a :ref:`configuration file <main_cfg>`, as well as paths to collections of the necessary input files, `mkrefs.py` can be called in the following way:

::

  mkrefs.py -vvv --cfgfile main_config_file.cfg example_readnoise_module /path/to/dark/files/*DARK*_uncal.fits /path/to/flatfield/files/NRCN815A-LIN-53650[5678]*uncal.fits

To call in this way, make sure that `mkrefs.py` is in your PYTHONPATH. Otherwise, prepend `python /path/to/mkrefs/` to the call above.

The `-vvv` flag sets the verbosity. More v's leads to more details printed to the screen.

`example_readnoise_module` is the reference file creation script that you wish to run. You can create more than one reference file type by providing a comma-separated list of reference file creation scripts.

You can provide as many input file lists as you need. *jwst_reffiles* combines the input file lists into a single list of files. It then examines each of the input files and determines what type of file it is (dark or flat. *jwst_reffiles* currently only supports these two file types.).
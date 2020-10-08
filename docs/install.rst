Installing JWST_REFFILES
========================

There are two main methods for installing jwst_reffiles. Note that for both of these installation methods, the main script, mkrefs.py, will be made available from the command line (i.e. it can be called from anywhere).

Installation for Users
----------------------

If you are not planning on doing any development work for jwst_reffiles but only to use it, then installation can be done by (optionally) creating a new environment, activating it, and then using pip to install jwst_reffiles.

1 (optional): conda create -n jwst-reffiles python
2. conda activate jwst-reffiles
3. pip install git+https://github.com/spacetelescope/jwst_reffiles.git

Installation for Developers
---------------------------

For those who will be working on development of jwst_reffiles, installation should be done by `cloning the repository from GitHub <https://github.com/spacetelescope/jwst_reffiles>`_ and then installing the local copy. It's recommended that you install jwst_reffiles in a new environment. Also, when installing, it's recommended that you install using the "-e" flag, which allows you to any changes to the code without reinstalling each time.

1. conda create -n jwst-reffiles python
2. conda activate jwst-reffiles
3. git clone https://github.com/spacetelescope/jwst_reffiles.git
4. cd jwst_reffiles
5. pip install -e .


Calling JWST_REFFILES
=====================

Currently, the top-level script in *jwst_reffiles*, which is called mkrefs.py, is designed such that it can only be called from the command line. Given a :ref:`configuration file <main_cfg>`, as well as paths to collections of the necessary input files, `mkrefs.py` can be called in the following way from any directory (mkrefs.py need not be in the working directory):

::

  mkrefs.py -vvv --cfgfile main_config_file.cfg example_readnoise_module /path/to/dark/files/*DARK*_uncal.fits /path/to/flatfield/files/NRCN815A-LIN-53650[5678]*uncal.fits

The `-vvv` flag sets the verbosity. More v's leads to more details printed to the screen.

`example_readnoise_module` is the reference file creation module that you wish to run. You can create more than one reference file type by providing a comma-separated list of reference file creation modules.

You can provide as many input file lists/paths as you need. *jwst_reffiles* combines the input file lists into a single list of files. It then examines each of the input files and determines what type of file it is (dark or flat. *jwst_reffiles* currently only supports these two file types.).
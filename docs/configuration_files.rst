Configuration Files
===================

Jwst_reffiles makes use of multiple configuration files to control the behavior of the main ``mkrefs.py`` script, as well as describe the reference files to be created and define the parameters needed for each.



.. _main_cfg:

An example of the main configuration file:

.. parsed-literal::

    instrument_: NIRCam
    reflabels_: [example_readnoise_module]
    reftypes_: [bpm,rdnoise,gain,linearity]

    default_reflabels_:
       bpm: bpm
       rdnoise: example_readnoise_module
       gain: gain_armin

    inputfiles:
       dark_pattern_: Bias|Dark|BIAS|DARK|ZERO
       flat_pattern_: GAIN|LIN

       basenamepattern_: (.*)\_uncal\.fits$

       # specificy the fits keywords to determine MJD time of
       # observation. They are automatically added to the list of required
       # fits keywords. if MJD-OBS does not exist, then DATE-OBS (and
       # TIME-OBS if specified) are used to determine the MJD.
       dateobs_fitskey_: DATE-END
       timeobs_fitskey_: TIME-END
       mjdobs_fitskey_:

       requiredfitskeys_: [INSTRUME,DETECTOR,NGROUPS,NINTS,TGROUP,SUBARRAY]
       optionalfitskeys_: [FILTER,PUPIL]

    DD_:
       max_Delta_MJD_: 2.0

    FF_:
       max_Delta_MJD_: 2.0
    DDFF_:
       max_Delta_MJD_: 10.0

    # reference file output basename_: outrootdir[/outsubdir][/runID]/reflabel/reflabel[_outsubdir][_runID][.addsuffix].cmdID.reftype.fits
    # if runID=AUTO, then the latest runID in the output dir outrootdir[/outsubdir] is taken
    output:
       outrootdir_: $JWST_MKREFS_OUTDIR
       outsubdir_: test1
       runID_: AUTO
       # how many leading zeros in runID
       runIDNdigits_: 3
       addsuffix_: bla

       # directory in which ssb products are saved.
       # if not specified, default is outrootdir[/outsubdir][/runID]/ssb
       # if skip_runID_ssbdir, then the runID is skipped
       ssbdir_:
       skip_runID_ssbdir: False
       # If ssblogFlag=True, the stdout from the strun command is saved in a log
       # file, same name as reduced filename, but with suffix log.txt
       ssblogFlag_: False
       # If ssberrorlogFlag=True, the stderr from the strun command is
       # saved in a log file, same name as reduced filename, but with
       # suffix err.txt
       ssberrorlogFlag_: False

       # directories in which ssb products are looked for. comma-separated list!
       pipeline_prod_search_dir_:

    batch:
       batch_enabled_: False
       batchmode_: Condor


    test1: 4

    # values allowed: CRDS, SELF, filepattern
    reffile4ssb_:
       gain: CRDS
       bpm: CRDS
       rdnoise: CRDS

    # values allowed: CRDS, filename
    validation_:
       gainreffile: CRDS
       bpmreffile: CRDS

    example_readnoise_module_:
        reftype_: rdnoise
        imtypes_: D
        ssbsteps_: dq_init
        verbose: 3
        boxsize: 128
        sigma_threshold: 4

    bpm:
        reftype: bpm
        imtypes: D
        ssbsteps: refpix-

    rdnoise_nircam:
        reftype: rdnoise
        imtypes: DD
        ssbsteps: dq_init
        test1: 6

    gain_armin:
        reftype: gain
        imtypes: DDFF
        ssbsteps: superbias-
        reffile4ssb:
            bpmreffile: /bla/bpm.fits

    gain_bryan:
        reftype: gain
        imtypes: FF
        ssbsteps: rate-
        reffile4ssb:
            bpmreffile: /bla/bpm.fits





.. _instrument:

Instrument
----------

The name of the instrument whose reference files are being created. Must be one of the JWST instruments. Capitalization is not important.

.. _reflabels:

Reference Labels
----------------

The list of allowed reference labels. These are simply the names of the reference file creation packages that can be run. In the near future, this list will be generated automatically by searching for the wrapper scripts around plugged-in modules. At that time, reflabels will be removed from the config file.

.. _reftypes:

Reference Types
---------------

List of the allowed reference file types that can be created. At a minimum, this list must contain the reference file types corresponding to all of the entries in the reflabels_ list.

.. _default_reflabels:

Default Reference Labels
------------------------

(NOT YET IMPLEMENTED) In this section, the user can define which of the reference file creation packages (reflabels_) is the default for each reference type. This is a convenience feature for the user. If you define `my_nircam_gain_script` as the default package to use for gain files, then you can call `mkrefs.py` from the command line and simply request `gain` rather than `my_nircam_gain_script`.

.. _dark_pattern:

Dark Pattern
------------

List of strings that will be used to identify dark current files within the list of input files. mkrefs will look for these strings in the filenames. Note that the values shown in the example configuration file above are designed around ground testing files. JWST data obtained in flight will not have unique features in the filenames that identify them as darks, flats, etc. Future improvements to *jwst_reffiles* will address this.

.. _flat_pattern:

Flat Pattern
------------

List of strings that will be used to identify flat field files within the list of input files. *jwst_reffiles* will look for these strings in the filenames. Note that the values shown in the example configuration file above are designed around ground testing files. JWST data obtained in flight will not have unique features in the filenames that identify them as darks, flats, etc. Future improvements to *jwst_reffiles* will address this.

.. _basenamepattern:

Basename Pattern
----------------

String to use when looking for input files. It may be that there are different versions (e.g. outputs from different points within the pipeline) of a file specified by the user in the input directory. If the basename pattern is set, *jwst_reffiles* will look only for the files matching that pattern. In the example above, we limit inputs to uncalibrated files.

.. _dateobs_fitskey:

Date-Obs Fits Header Keyword
----------------------------

The header keyword in the input files that contains the date of the observation. Dates are needed to enforce rules when pairing darks and flats. For JWST, set this to **DATE-END**.

.. _timeobs_fitskey:

Time-Obs Fits Header Keyword
----------------------------

The header keyword in the input files that contains the time of the observation. Times are needed to enforce rules when pairing darks and flats. For JWST, set this to **TIME-END**.

.. _mjdobs_fitskey:

Obervation MJD Fits Header Keyword
----------------------------------

The header keyword in the input files that contains the time of the observation in MJD. If specified, the value associated with this keyword is used to determinte the modified julian date (MJD) of the observation.

.. _requiredfitskeys:

Required Fits File Header Keywords
----------------------------------

List of header keywords which must be present in the input files. Values from these keywords will be copied into the master table created by *jwst_reffiles*. If any input files do not contain all of these keywords, an error will be raised.

.. _optionalfitskeys:

Optional FitsFile Header Keywords
---------------------------------

List of optional header keywords in the input files. Values from these keywords will be copied into the master table created by *jwst_reffiles*. If any of the keywords are missing in any of the input files, they are simply not copied to the master table, and *jwst_reffiles* proceeds without raising an error.

.. _DD:

DD
---

Stands for dark dark, and indicates a situation where a pair of dark current ramps are needed to produce an individual reference file.

.. _FF:

FF
---

Stands for flat flat, and indicates a situation where a pair of flat field ramps are needed to produce an individual reference file.

.. _DDFF:

DDFF
----

Stands for dark dark flat flat, and indicates a situation where a pair of dark current ramps and a pair of flat field ramps are needed to produce an individual reference file.

.. _max_Delta_MJD:

Max_Delta_MJD
-------------
The maximum time allowed, in days, between input observations when creating pairs/groups of files. For example, to minimize the chances of dark current varying enough to impact reference file creation, you can set `max_Delta_MJD` to (e.g.) 2. In this case, when pairing dark current files (i.e. DD_) *jwst_reffiles* will not pair observations taken more than 2 days apart.

.. _basename:

Reference File Output Basename
------------------------------

Format of the output names for individual reference files. Output names will be automatically generated by *jwst_reffiles* to ensure accurate bookkeeping. The overall format of the reference file output names follows the convention shown below. **cmdID** is an ID assigned by *jwst_reffiles* for a particular run of the package. This helps to ensure unique file and directory names for the outputs.

outrootdir_[/outsubdir_][/runID_]/reflabel/reflabel[_outsubdir_][_runID_][.addsuffix_].cmdID.reftype_.fits


.. _outrootdir:

Output Root Directory
---------------------

Path to the top level output directory for *jwst_reffiles*. Default is to define this within the **JWST_MKREFS_OUTDIR** environment variable, but any valid path is acceptable.

.. _outsubdir:

Output Subdirectory
-------------------

Subdirectory name to add to the Output Root Directory when creating outputs.

.. _runID:

Run ID
------

An integer that will be used to create a unique subdirectory for *jwst_reffiles* outputs for a given "run" of the software. With this parameter, it is easy to organize the outputs and prevent conflicts for multiple runs using the same input files. Leaving this set to the default value of **AUTO** will cause *jwst_reffiles* to search for the most recent/highest existing run ID, and add 1 for the next run.

.. _runIDNdigits:

Number of Digits in the Run ID
------------------------------

The total number of digits in the run ID. Leading zeros are added as necessary. The default is 3.

.. _addsuffix:

Add Suffix
----------

Optional suffix that can be added to the output files from *jwst_reffiles*.

.. _ssbdir:

SSB Directory
-------------

Directory where JWST calibration pipeline output files, which are often created in the process of running *jwst_reffiles*, will be saved. If this entry is left blank, the default value of outrootdir_[/outsubdir_][/runID_]/ssb/ will be used.

.. _skip_runID_ssbdir:

Skip RunID for SSB Directory
----------------------------

Boolean. If True, the output directory (ssbdir_) for the pipeline output files will follow the convention outrootdir_[/outsubdir_]/ssb/.

.. _ssblogFlag:

SSB Log Flag
------------

If ssblogFlag=True, the stdout from the JWST calibration pipeline call is saved in a log file, same name as pipeline output filename, but with suffix log.txt.

.. _ssberrorlogFlag:

SSB Error Log Flag
------------------

If ssberrorlogFlag=True, the stderr from the JWST calibration pipeline call is saved in a log file, same name as reduced filename, but with suffix err.txt.

.. _pipeline_prod_search_dir:

Pipeline Products Search Directories
------------------------------------

This is comma-separated list of directories. Prior to calling the calibration pipeline for a given input fits file, *jwst_reffiles* will search the directories listed here to see if the proper pipeline-processed file already exists. If blank, no searching will be done.

.. _batch_enabled:

Batch Enabled
-------------

Boolean entry. If True, *jwst_reffiles* will be run in batch mode. (BATCH MODE NOT YET IMPLEMENTED)

.. _batchmode:

Batch Mode
----------

The batch system to use when running in batch mode. Default is Condor.

.. _reffile4ssb:

Reference Files for SSB
-----------------------

NOT YET IMPLEMENTED. In this section, list the types of reference files to use in the calls to the calibration pipeline. Options are **CRDS**, **SELF**, or a file pattern. If **CRDS** is used, then the appropriate reference files will be selected from the Calibration Reference Data System (`CRDS <https://jwst-pipeline.readthedocs.io/en/stable/jwst/introduction.html#crds>`_). This system contains only officially delivered reference files.

If **SELF** is used, then calls to the calibration pipeline will use reference files generated from the current run of *jwst_reffiles*. Note that in this case, running *jwst_reffiles* becomes an iterative process. For example, run once to produce a superbias reference file. Then run again to use this superbias reference file to calibrate inputs when creating dark current reference files.

Finally, you can enter a file pattern (e.g. /my/files/reffiles/gain/*gain.fits), in which case *jwst_reffiles* will use that file.

.. _validation:

Validation
----------

NOT YET IMPLEMENTED. Can be either **CRDS** or a filename. When comparing an output reference file to a previous version, this controls where the comparison file comes from. If set to **CRDS**, the most recent matching file in CRDS is used for comparison. If set to a filename, that file is used.

.. _example_readnoise_module:

Example Readnoise Module (example of a reference file creation module)
----------------------------------------------------------------------

This point in the config file is where you define the options that are specific to each reference file creation module. Any modules that are going to be run must be listed in the reflabels_ list. Any modules not in the reflabels_ list will be ignored.

.. _reftype:

Reference Type
--------------

Define the type of reference file created by this module. (e.g. gain, rdnoise)

.. _imtypes:

Image Types
-----------

The type of inputs required by this package. (e.g. DD for a pair of dark ramps. FF for a pair of flat field ramps. DDFF for a pair of each.)

.. _ssbsteps:

SSB Steps
---------

Comma-separated list of calibration pipeline steps which must be complete on the input files prior to creating the reference file. Note that if the input files have been partially processed by the pipeline, the full list of completed steps must still be given here. For convenience, there is also a "-" shorthand that can be used. If the input files require all pipeline steps up to and including dark current subtraction, then you can enter "dark-". The pipeline steps currently recognized by *jwst_reffiles* includes all of those in the calwebb_detetor1 pipeline, and are called using the values from the list below. Other values will not be recognized.

::
  group_scale, dq_init, saturation, ipc, firstframe, lastframe, superbias, refpix, linearity, persistence, rscd, dark_current, jump, rate.

Example entries:

dq_init, saturation, refpix
refpix-

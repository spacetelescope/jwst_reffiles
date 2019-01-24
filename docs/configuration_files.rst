Configuration Files
===================

Jwst_reffiles makes use of multiple configuration files to control the behavior of the main ``mkrefs.py`` script, as well as describe the reference files to be created and define the parameters needed for each.



.. _main_cfg:

An example of the main configuraiton file:

::

    instrument: NIRCam
    reflabels: [example_readnoise_module]
    reftypes: [bpm,rdnoise,gain,linearity]

    default_reflabels:
       bpm: bpm
       rdnoise: example_readnoise_module
       gain: gain_armin

    inputfiles:
       dark_pattern: Bias|Dark|BIAS|DARK|ZERO
       flat_pattern: GAIN|LIN

       basenamepattern: (.*)\_uncal\.fits$

       # specificy the fits keywords to determine MJD time of
       # observation. They are automatically added to the list of required
       # fits keywords. if MJD-OBS does not exist, then DATE-OBS (and
       # TIME-OBS if specified) are used to determine the MJD.
       dateobs_fitskey: DATE-END
       timeobs_fitskey: TIME-END
       mjdobs_fitskey:

       requiredfitskeys: [INSTRUME,DETECTOR,NGROUPS,NINTS,TGROUP,SUBARRAY]
       optionalfitskeys: [FILTER,PUPIL]

    DD:
       max_Delta_MJD: 2.0

    FF:
       max_Delta_MJD: 2.0
    DDFF:
       max_Delta_MJD: 10.0

    # reference file output basename: outrootdir[/outsubdir][/runID]/reflabel/reflabel[_outsubdir][_runID][.addsuffix].cmdID.reftype.fits
    # if runID=AUTO, then the latest runID in the output dir outrootdir[/outsubdir] is taken
    output:
       outrootdir: $JWST_MKREFS_OUTDIR
       outsubdir: test1
       runID: AUTO
       # how many leading zeros in runID
       runIDNdigits: 3
       addsuffix: bla

       # directory in which ssb products are saved.
       # if not specified, default is outrootdir[/outsubdir][/runID]/ssb
       # if skip_runID_ssbdir, then the runID is skipped
       ssbdir:
       skip_runID_ssbdir: False
       # If ssblogFlag=True, the stdout from the strun command is saved in a log
       # file, same name as reduced filename, but with suffix log.txt
       ssblogFlag: False
       # If ssberrorlogFlag=True, the stderr from the strun command is
       # saved in a log file, same name as reduced filename, but with
       # suffix err.txt
       ssberrorlogFlag: False

       # directories in which ssb products are looked for. comma-separated list!
       pipeline_prod_search_dir:

    batch:
       batch_enabled: False
       batchmode: Condor


    test1: 4

    # values allowed: CRDS, SELF, filepattern
    reffile4ssb:
       gain: CRDS
       bpm: CRDS
       rdnoise: CRDS

    # values allowed: CRDS, filename
    validation:
       gainreffile: CRDS
       bpmreffile: CRDS

    example_readnoise_module:
        reftype: rdnoise
        imtypes: D
        ssbsteps: dq_init
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


#!/usr/bin/env python
'''
create reference files for JWST instruments
A. Rest
'''
import argparse
import glob
import os
import re
import sys
import types

import astropy.io.fits as fits
from astropy.io import ascii
import astropy
import numpy as np
import scipy

from jwst_reffiles.mkref import mkrefclass
from jwst_reffiles.pipeline.calib_prep import CalibPrep
from jwst_reffiles.pipeline import pipeline_steps
from jwst_reffiles.utils.tools import astrotableclass, yamlcfgclass
from jwst_reffiles.utils.tools import makepath

# get the root dir of the code. This is needed only if the scripts are not installed as a module!
#if 'JWST_MKREFS_SRCDIR' in os.environ:
#    rootdir = os.environ['JWST_MKREFS_SRCDIR']
#    sys.path.extend([rootdir,
#                     '%s/gain' % rootdir,
#                     '%s/badpix_map' % rootdir])
#else:
rootdir = os.path.dirname(os.path.realpath(__file__))


class cmdsclass(astrotableclass):
    def __init__(self):
        astrotableclass.__init__(self)


class mkrefsclass(astrotableclass):
    def __init__(self):
        astrotableclass.__init__(self)

        # config file
        self.cfg = None

        self.verbose = 0
        self.debug = False
        self.onlyshow = False

        self.basedir = None
        self.basename = None
        self.outsubdir = None
        self.ssbdir = None
        self.runID = None
        self.runlabel = None

        self.imtable = astrotableclass()
        self.images4ssb = astrotableclass()
        #self.darks = None
        #self.flats = None

        #
        self.DDtable = astrotableclass()
        self.FFtable = astrotableclass()
        self.DDFFtable = astrotableclass()

        self.cmdtable = cmdsclass()

        self.allowed_reflabels = ['bpm', 'rdnoise_nircam', 'gain_armin']

    def define_options(self, parser=None, usage=None, conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
        parser.add_argument("reflabels_and_imagelist", nargs='+', help=("list of ref types to be done"
                                                                        "and image (or file patterns) lists"))
        parser.add_argument('-b', '--batchmode', help="run the commands in batch mode (default=%(default)s)",
                            action="store_true", default=False)
        parser.add_argument('--onlyshow', help="only show what would be done, but don't do it."
                                               "(default=%(default)s)", action="store_true", default=False)

        parser.add_argument('--outrootdir', default=None,
                            help=('output root directory. outrootdir[/outsubdir][/runID]/reflabel/'
                                  'reflabel[_outsubdir][_runID][.addsuffix].cmdID.reftype.fits '
                                  '(default=%(default)s)'))
        parser.add_argument('--outsubdir', default=None,
                            help=('subdir added to the output root directory (and filename) '
                                  '(default=%(default)s)'))
        parser.add_argument('--addsuffix', default=None,
                            help='suffix added to the output basename (default=%(default)s)')

        parser.add_argument('--runID', default=None,
                            help=('runID subdir added to output root directory (and filename) '
                                  '(default=%(default)s)'))
        parser.add_argument('--runIDNdigits', default=None,
                            help='how many leading zeros in runID (default=%(default)s)')
        parser.add_argument('-n', '--newrunID', help="get the next available runID",
                            action="store_true", default=False)
        parser.add_argument('--skip_runID_ssbdir', help="Don't use runID in the default ssb dir",
                            action="store_true", default=False)

        parser.add_argument('--ssbdir', default=None,
                            help='Directory into which ssb files are stored (default=%(default)s)')
        #parser.add_argument('--imdates2subdir', help=("add firstdate_lastdate subdir of the input images "
        #                                              "used (default=%(default)s)"), action="store_true",default=False)

        # Loop through all allowed reflabels, and add the options
        for reflabel in self.allowed_reflabels:
            print(reflabel)
            mkrefpackage = __import__("mkref_%s" % reflabel)
            mkref = mkrefpackage.mkrefclassX()
            # Only load defaultoptions once...
            if reflabel == self.allowed_reflabels[0]:
                mkref.defaultoptions(parser)
            # get the reflabel specific options
            mkref.extraoptions(parser)
        return(parser)

    def loadcfgfiles(self, maincfgfile, extracfgfiles=None, params=None, params4all=None,
                     params4sections=None, requireParamExists=True):
        if self.cfg is None:
            self.cfg = yamlcfgclass()
        if self.cfg.loadcfgfiles(maincfgfile, extracfgfiles=extracfgfiles,
                                 params=params, params4all=params4all, params4sections=params4sections,
                                 requireParamExists=requireParamExists, verbose=self.verbose):
            raise RuntimeError("Something went wrong when loading config files!")
        return(0)

    def getbasedir(self, outrootdir=None, outsubdir=None, runlabel=None):
        """ basedir based on options and cfg file: outrootdir[/outsubdir][/runID]"""
        if outrootdir is None:
            if self.cfg.params['output']['outrootdir'] == '':
                raise RuntimeError("No output rootdir!")
            else:
                outrootdir = self.cfg.params['output']['outrootdir']

        basedir = os.path.abspath(os.path.expandvars(outrootdir))
        basename = ''

        # outsubdir: test if outsubdir is passed or if it is specified in the config file ...
        if outsubdir is None and self.cfg.params['output']['outsubdir'] != '':
            outsubdir = self.cfg.params['output']['outsubdir']
        # ... and if yes, add it
        if outsubdir is not None:
            basedir += '/%s' % outsubdir
            if basename != '':
                basename += '_'
            basename += outsubdir

        if runlabel is not None:
            basedir += '/%s' % runlabel
            if basename != '':
                basename += '_'
            basename += runlabel

        if basename == '':
            basename = 'B'

        return(basedir, '%s/%s' % (basedir, basename), outsubdir)

    def getbasename(self):
        s = self.basedir+'/'
        if self.outsubdir is not None:
            s += '%s/' % self.outsubdir
        if self.runID is not None:
            s += '%s/' % self.outsubdir

    def getssbdir(self, ssbdir, basedir, skip_runID_ssbdir=False):
        """ get the ssbdir. It can be hardcoded to a directory with
        --ssbdir or in the cfg file. If it is not specified, then
        basedir/ssb is used"""
        if ssbdir is not None and ssbdir != '':
            return(ssbdir)
        elif self.cfg.params['output']['ssbdir'] is not None and self.cfg.params['output']['ssbdir'] != '':
            return(self.cfg.params['output']['ssbdir'])
        else:
            return('%s/ssb' % self.basedir)
        return(None)

    def get_highest_runID(self, basedir):
        """ find all subdirs witn run\d+ within the passed basedir, and extract the highest runID"""
        files = glob.glob('%s/run*' % basedir)
        if len(files) == 0:
            print('No runIDs yet!')
            return(None)
        runIDs = []
        runIDpattern = re.compile('run(\d+)$')
        for fname in files:
            print(os.path.basename(fname))
            m = runIDpattern.search(fname)
            if m is not None:
                runID = int(m.groups()[0])
                runIDs.append(runID)
        if len(runIDs) > 0:
            runID = max(runIDs)
        else:
            runID = None
        return(runID)

    def set_runID(self, outrootdir=None, outsubdir=None, runID=None, runIDNdigits=None, newrunID=False):
        """ set the runID depending on the options and config file """

        # if no runID is specified with the options, then use the one from the cfg file
        if runID is None:
            runID = mkrefs.cfg.params['output']['runID']

        if runID == '' or runID is None:
            # If no runID specified, make sure it is set to None!
            self.runID = None
        elif isinstance(runID, int) or runID.isdigit():
            # if runID is an integer, set it to it!
            self.runID = int(runID)
        else:
            # if runID is AUTO, get the current highest runID, and increment if --new_runID (-n)
            if runID.upper() == 'AUTO':
                (basedir, basename, subdir) = self.getbasedir(outrootdir=outrootdir, outsubdir=outsubdir,
                                                              runlabel=None)
                self.runID = self.get_highest_runID(basedir)
                if self.runID is None:
                    self.runID = 0
                else:
                    if newrunID:
                        self.runID += 1
            else:
                raise RuntimeError("Incorrect runID=%s requested!" % runID)

        if self.runID is not None:
            if runIDNdigits is None:
                runIDNdigits = self.cfg.params['output']['runIDNdigits']
            if runIDNdigits is None:
                runIDNdigits = 1

            b = 'run%0'+'%d' % runIDNdigits+'d'
            self.runlabel = b % (self.runID)
        else:
            self.runlabel = ''

        return(self.runID, self.runlabel)

    #def getssbdir(self, ssbdir=None, **kwarg):
    #    if ssbdir != None and ssbdir != '':
    #        return(ssbdir)
    #    ssbdir = '%s/ssb' % self.getbasedir(**kwarg)
    #    return(ssbdir)

    def set_dirs(self, outrootdir=None, outsubdir=None, runID=None, runIDNdigits=1, newrunID=False,
                 ssbdir=None, skip_runID_ssbdir=False):
        # get the runID
        mkrefs.set_runID(outrootdir=outrootdir, outsubdir=outsubdir,
                         runID=runID, runIDNdigits=runIDNdigits, newrunID=newrunID)

        (self.basedir, self.basename, self.outsubdir) = self.getbasedir(outrootdir=outrootdir,
                                                                        outsubdir=outsubdir,
                                                                        runlabel=self.runlabel)

        self.ssbdir = self.getssbdir(ssbdir, self.basedir, skip_runID_ssbdir=skip_runID_ssbdir)

        # Create the directories if they don't already exist
        if not self.onlyshow:
            state_basedir = makepath(self.basedir)
            state_ssbdir = makepath(self.ssbdir)

        if self.verbose >= 3:
            print('#### Directories:')
            print('basedir:', self.basedir)
            print('basename:', self.basename)
            print('ssbdir:', self.ssbdir)
            print('outsubdir:', self.outsubdir)
            print('runID:', self.runID)
            print('runlabel:', self.runlabel)

        return(0)

    def setoutbasename(self, reftype, imagelist, outbasename=None, outrootdir=None, outsubdir=None,
                       addsuffix=None):
        #if outbasename is not None:
        #    self.outbasename = outbasename
        #    return(outbasename)

        if len(imagelist) > 1:
            print('### more than 1 input file in setoutbasename not yet implemented!!!! exiting .... ###')
            sys.exit(0)

        # if outdir is passed or specified, use it, if not use the directiroy of the input filename
        (imdir, imbasename) = os.path.split(os.path.abspath(imagelist[0]))
        if outrootdir is None:
            if self.cfg.params['output']['outrootdir'] == '':
                raise RuntimeError("No output rootdir!")
            else:
                outrootdir = self.cfg.params['output']['outrootdir']

        outbasename = '%s' % outrootdir

        # outsubdir: test if outsubdir is passed or if it is specified in the config file ...
        if outsubdir is None and self.cfg.params['output']['outsubdir'] != '':
            outsubdir = self.cfg.params['output']['outsubdir']
        # ... and if yes, add it
        if outsubdir is not None:
            outbasename += '/%s' % outsubdir

        outbasename += '/'
        if reftype is not None:
            outbasename += '%s_' % (reftype)
        outbasename += re.sub('\.fits$', '', imbasename)

        # addsuffix: test if addsuffix is passed or if it is specified in the config file ...
        if addsuffix is None and self.cfg.params['output']['addsuffix'] != '':
            addsuffix = self.cfg.params['output']['addsuffix']
        # ... and if yes, add it
        if addsuffix is not None:
            if not re.search('\.', addsuffix):
                outbasename += '.'
            outbasename += '%s' % addsuffix

        if self.verbose > 1:
            print('Output basenbame:', outbasename)

        self.outbasename = outbasename

        return(outbasename)

    def trim_imagelist(self, imagelist, basenamepattern=None):
        '''
        only keep fits file that match basenamepattern
        '''
        if self.verbose > 1:
            print("#### Trimming images...")
        # print imagelist
        if imagelist is None or len(imagelist) == 0:
            print('Nothing to do, no images!!!')
            return(imagelist)

        if basenamepattern is not None:
            newimagelist = []
            basenamepattern_compiled = re.compile(basenamepattern)
            # make sure the input files are all ok!
            for i in range(len(imagelist)):
                m = basenamepattern_compiled.search(imagelist[i])
                if m is None:
                    print('SKIPPING', imagelist[i])
                else:
                    #if len(m.groups())==0:
                    #    raise RuntimeError,"%s is matching %s, but does not return basename" % (basenamepattern_compiled.pattern,imagelist[i])
                    #basename = m.groups()[0]
                    #print 'VVV',imagelist[i]
                    newimagelist.append(imagelist[i])
            if len(imagelist) != len(newimagelist):
                if self.verbose > 1:
                    print('skipping {} out of {} images, {} left'.format(len(imagelist)-len(newimagelist),
                                                                         len(imagelist), len(newimagelist)))
            imagelist = newimagelist

        return(imagelist)

    def parse_reflabels_images(self, reflabels_and_imagelist, basenamepattern=None):
        reflabellist = []
        imagelist = []
        for s in reflabels_and_imagelist:
            if s in self.cfg.params['reflabels']:
                reflabellist.append(s)
            else:
                if not os.path.isfile(s):
                    raise RuntimeError("ERROR: file %s does not exist, thus not a viable input file" % s)
                imagelist.append(os.path.abspath(s))

        imagelist = self.trim_imagelist(imagelist, basenamepattern)
        return(reflabellist, imagelist)

    def getimtypes(self):
        if not ('imtype' in self.imtable.t.colnames):
            self.imtable.t['imtype'] = None

        darkpattern = re.compile(self.cfg.params['inputfiles']['dark_pattern'])
        flatpattern = re.compile(self.cfg.params['inputfiles']['flat_pattern'])

        for i in range(len(self.imtable.t)):
            shortfilename = os.path.basename(self.imtable.t['fitsfile'][i])
            if darkpattern.search(shortfilename):
                self.imtable.t['imtype'][i] = 'dark'
            elif flatpattern.search(shortfilename):
                self.imtable.t['imtype'][i] = 'flat'
            else:
                raise RuntimeError('ERROR: image type of image %s is unknown!' % shortfilename)

    def getimageinfo(self, imagelist, dateobsfitskey=None, timeobsfitskey=None, mjdobsfitskey=None):

        if self.verbose > 1:
            print("Getting image info")

        # self.imtable['fitsfile'].format('%s')
        self.imtable.t['fitsfile'] = imagelist
        self.imtable.t['imID'] = list(range(len(imagelist)))
        self.imtable.t['imtype'] = None
        self.imtable.t['skip'] = False

        if self.verbose > 1:
            print("Table created")

        requiredfitskeys = self.cfg.params['inputfiles']['requiredfitskeys']
        if requiredfitskeys is None:
            requiredfitskeys = []
        if type(requiredfitskeys) == bytes:
            requiredfitskeys = [requiredfitskeys]
        if dateobsfitskey is not None:
            requiredfitskeys.append(dateobsfitskey)
        if timeobsfitskey is not None:
            requiredfitskeys.append(timeobsfitskey)
        if mjdobsfitskey is not None:
            requiredfitskeys.append(mjdobsfitskey)
            mjdcol = mjdobsfitskey
        else:
            mjdcol = 'MJD'

        self.imtable.fitsheader2table('fitsfile',
                                      requiredfitskeys=requiredfitskeys,
                                      optionalfitskey=self.cfg.params['inputfiles']['optionalfitskeys'],
                                      raiseError=True, skipcolname='skip')
        self.imtable.dateobs2mjd(dateobsfitskey, mjdcol, mjdobscol=mjdobsfitskey, timeobscol=timeobsfitskey)

        self.getimtypes()
        # sort by MJD
        self.imtable.t.sort('MJD')

        #self.darks = self.imtable.t[np.where(self.imtable.t['imtype']=='dark')]
        #self.flats = self.imtable.t[np.where(self.imtable.t['imtype']=='flat')]
        return(0)

    def organize_inputfiles(self, reflabels_and_imagelist):
        # parse teh command line arguments for reflabels and images
        (reflabellist, imagelist) = self.parse_reflabels_images(reflabels_and_imagelist,
                                                                basenamepattern=self.cfg.params['inputfiles']['basenamepattern'])
        self.reflabellist = reflabellist

        # make the image table and populate it with info. Also get teh darks and flats table
        self.getimageinfo(imagelist,
                          dateobsfitskey=self.cfg.params['inputfiles']['dateobs_fitskey'],
                          timeobsfitskey=self.cfg.params['inputfiles']['timeobs_fitskey'],
                          mjdobsfitskey=self.cfg.params['inputfiles']['mjdobs_fitskey'])

        self.detectors = set(self.imtable.t['DETECTOR'])

        if self.verbose:
            print('#################\n### %d images found!' % len(self.imtable.t))
            print('### %d darks, %d flats' % (len(np.where(self.imtable.t['imtype'] == 'dark')[0]),
                                              len(np.where(self.imtable.t['imtype'] == 'flat')[0])))
            print('### %d detectors:' % (len(self.detectors)), ", ".join(self.detectors))
            if self.verbose > 1:
                print(self.imtable.t)

    def check_inputfiles(self):
        print('basedir:', self.basedir)
        print('basename:', self.basename)

        # test if output dir already exists.
        #if not os.path.isdir(self.basedir):
        #    # if it doesn't exist: write the image list into it if !onlyshow!
        #    if not self.onlyshow:
        #        if self.verbose > 1:
        #            print('Creating directory %s' % self.basedir)
        #        os.makedirs(self.basedir)
        #        if not os.path.isdir(self.basedir):
        #            raise RuntimeError('Cannot create directory %s' % self.basedir)
        #    else:
        #        print('*** only showing: directory %s would be created here!' % self.basedir)

        imtbl_filename = '%s.im.txt' % self.basename

        # check: does the saved inputim file has the same input images
        # than the current list? If yes, then all is good, if not
        # throw an error!
        if os.path.isfile(imtbl_filename):
            imtable_before = astrotableclass()
            imtable_before.load(imtbl_filename)

            if (len(self.imtable.t['fitsfile']) != len(imtable_before.t['fitsfile'])):
                raise RuntimeError(('length %d of new list of input images is different than the one from '
                                    'the saved table (%d), cannot proceed! either get a new runID with -n, '
                                    'or remove the files in %s' % (len(self.imtable.t['fitsfile']),
                                                                   len(imtable_before.t['fitsfile']),
                                                                   self.basedir)))

            for i in range(len(self.imtable.t['fitsfile'])):
                if self.imtable.t['fitsfile'][i] != imtable_before.t['fitsfile'][i]:
                    raise RuntimeError(('input file %s in new list is different than previous input '
                                        'file %s' % (self.imtable.t['fitsfile'][i],
                                                     imtable_before.t['fitsfile'][i])))

            if self.verbose:
                print('new input list matches previous input list! Continuing...')

        self.imtable.write(imtbl_filename, verbose=True, clobber=True)

        # just some paranoia...
        if not os.path.isfile(imtbl_filename):
            raise RuntimeError('Could not write %s! permission issues?' % (imtbl_filename))

        return(0)

    def get_optional_arguments(self, args, sysargv):
        fitspattern = re.compile('\.fits$')

        opt_arg_list = []
        for i in range(1, len(sysargv)):
            if sysargv[i] in args.reflabels_and_imagelist:
                if sysargv[i] in self.reflabellist:
                    print('this is a reflabel', sysargv[i])
                else:
                    # test if it is an image
                    if not fitspattern.search(sysargv[i]):
                        print('not a fits file! filepattern?', sysargv[i])
            else:
                opt_arg_list.append(sysargv[i])
        print('optional arguments:', opt_arg_list)
        return(opt_arg_list)

    def getDlist(self, detector):
        '''
        returns list of Dark indeces, where the indices refer to the self.darks table
        '''

        #self.imtable.t['skip'][7]=True
        # indices for dark frames
        dindex, = np.where(self.imtable.t['imtype'] == 'dark')
        # indices for detector and not skipped
        dindex4detector = dindex[np.where(np.logical_and(self.imtable.t['DETECTOR'][dindex] == detector,
                                                         np.logical_not(self.imtable.t['skip'][dindex])))]

        if self.verbose > 2:
            print('Possible %d Darks for detector %s' % (len(dindex4detector), detector))
            print(self.imtable.t[dindex4detector])

        #D = astrotableclass(names=('D1index','D1imID'),dtype=('i4', 'i4'))
        D = astrotableclass()
        D.t['D1index'] = dindex4detector
        D.t['D1imID'] = self.imtable.t['imID'][dindex4detector]
        D.t['D1index', 'D1imID'].format = '%d'
        #print D.t
        #sys.exit(0)
        return(D, ('D1',))

    def getDDlist(self, detector, max_Delta_MJD=None):
        '''
        returns list of Dark-Dark pair indeces,  where the indices refer to the self.imtable table
        '''
        if self.verbose > 1:
            print('# Getting DD list')
        #self.imtable.t['skip'][7]=True
        #print  self.imtable.t[6:11]

        # indices for dark frames
        dindex, = np.where(self.imtable.t['imtype'] == 'dark')
        # indices for detector and not skipped
        dindex4detector = dindex[np.where(np.logical_and(self.imtable.t['DETECTOR'][dindex] == detector,
                                                         np.logical_not(self.imtable.t['skip'][dindex])))]
        if self.verbose > 2:
            print('Possible %d Darks for detector %s' % (len(dindex4detector), detector))
            print(self.imtable.t[dindex4detector])

        DD = astrotableclass(names=('D1index', 'D2index', 'D1imID', 'D2imID'), dtype=('i4', 'i4', 'i4', 'i4'))
        i = 0
        while i < len(dindex4detector)-1:
            if max_Delta_MJD is not None:
                if self.verbose > 2:
                    print('Checking if imID=%d and %d can be DD' % (self.imtable.t['imID'][dindex4detector[i]],
                                                                    self.imtable.t['imID'][dindex4detector[i+1]]))
                dMJD = self.imtable.t['MJD'][dindex4detector[i+1]]-self.imtable.t['MJD'][dindex4detector[i]]
                print('dMJD:', dMJD)
                if dMJD > max_Delta_MJD:
                    if self.verbose > 2:
                        print(('Skipping imID=%d (MJD=%f) since imID=%d is not within timelimit (Delta '
                               'MJD = %f>%f)!' % (self.imtable.t['imID'][dindex4detector[i]],
                                                  self.imtable.t['MJD'][i],
                                                  self.imtable.t['imID'][dindex4detector[i+1]],
                                                  dMJD, max_Delta_MJD)))
                    i += 1
                    continue
            if self.verbose > 1:
                print('Adding DD pair with imID=%d and %d' % (self.imtable.t['imID'][dindex4detector[i]],
                                                              self.imtable.t['imID'][dindex4detector[i+1]]))
            DD.t.add_row({'D1index': dindex4detector[i], 'D2index': dindex4detector[i+1],
                          'D1imID': self.imtable.t['imID'][dindex4detector[i]],
                          'D2imID': self.imtable.t['imID'][dindex4detector[i+1]]})
            i += 2

        if self.verbose > 2:
            print(DD.t)
        return(DD, ('D1', 'D2'))

    def getFFlist(self, detector, max_Delta_MJD=None):
        '''
        returns list of Flat-Flat pair indeces, where the indices refer to the self.imtable table
        '''
        if self.verbose > 1:
            print('# Getting FF list')
        #self.imtable.t['skip'][7]=True
        #print  self.imtable.t[6:11]

        # indices for dark frames
        findex, = np.where(self.imtable.t['imtype'] == 'flat')
        # indices for detector and not skipped
        findex4detector = findex[np.where(np.logical_and(self.imtable.t['DETECTOR'][findex] == detector,
                                                         np.logical_not(self.imtable.t['skip'][findex])))]
        if self.verbose > 2:
            print('Possible %d Flats for detector %s' % (len(findex4detector), detector))
            print(self.imtable.t[findex4detector])

        FF = astrotableclass(names=('F1index', 'F2index', 'F1imID', 'F2imID'), dtype=('i4', 'i4', 'i4', 'i4'))
        i = 0
        while i < len(findex4detector)-1:
            if max_Delta_MJD is not None:
                if self.verbose > 2:
                    print('Checking if imID=%d and %d can be FF' % (self.imtable.t['imID'][findex4detector[i]],
                                                                    self.imtable.t['imID'][findex4detector[i+1]]))
                dMJD = self.imtable.t['MJD'][findex4detector[i+1]]-self.imtable.t['MJD'][findex4detector[i]]
                print('dMJD:', dMJD)
                if dMJD > max_Delta_MJD:
                    if self.verbose > 2:
                        print(('Skipping imID=%d (MJD=%f) since imID=%d is not within '
                               'timelimit (Delta MJD = %f>%f)!') % (self.imtable.t['imID'][findex4detector[i]],
                                                                    self.imtable.t['MJD'][i],
                                                                    self.imtable.t['imID'][findex4detector[i+1]],
                                                                    dMJD, max_Delta_MJD))
                    i += 1
                    continue
            if self.verbose > 1:
                print('Adding FF pair with imID=%d and %d' % (self.imtable.t['imID'][findex4detector[i]],
                                                              self.imtable.t['imID'][findex4detector[i+1]]))
            FF.t.add_row({'F1index': findex4detector[i], 'F2index': findex4detector[i+1],
                          'F1imID': self.imtable.t['imID'][findex4detector[i]],
                          'F2imID': self.imtable.t['imID'][findex4detector[i+1]]})
            i += 2

        if self.verbose > 2:
            print(FF.t)
        return(FF, ('F1', 'F2'))

    def getDDFFlist(self, detector, DD_max_Delta_MJD=None, FF_max_Delta_MJD=None, DDFF_max_Delta_MJD=None):
        '''
        returns list of Flat-Flat pair indeces, where the indices refer to the self.imtable table
        '''
        if self.verbose > 0:
            print('\n### Getting DDFF list')
        DD, DDimtypes = self.getDDlist(detector, max_Delta_MJD=DD_max_Delta_MJD)
        FF, FFimtypes = self.getFFlist(detector, max_Delta_MJD=FF_max_Delta_MJD)
        if self.verbose > 0:
            print('### DD and FF lists created, no matching them!!!')

        DDFF = astrotableclass(names=('F1index', 'F2index', 'F1imID', 'F2imID', 'D1index', 'D2index', 'D1imID',
                                      'D2imID'), dtype=('i4', 'i4', 'i4', 'i4', 'i4', 'i4', 'i4', 'i4'))
        ddcount = np.zeros(len(DD.t))

        for f in range(len(FF.t)):
            if self.verbose > 2:
                print('# Finding DD pair for FF pair with imID=%d and %d' % (FF.t['F1imID'][f], FF.t['F2imID'][f]))
            if DDFF_max_Delta_MJD is not None:
                FF_MJDmin = np.amin(np.array((self.imtable.t['MJD'][FF.t['F1index'][f]],
                                              self.imtable.t['MJD'][FF.t['F2index'][f]])))
                FF_MJDmax = np.amax(np.array((self.imtable.t['MJD'][FF.t['F1index'][f]],
                                              self.imtable.t['MJD'][FF.t['F2index'][f]])))
                if self.verbose > 3:
                    print('FF MJDs:', FF_MJDmin, FF_MJDmax)

            d_best = None
            ddcountmin = None
            for d in range(len(DD.t)):
                if ddcountmin is None or ddcount[d] < ddcountmin:
                    if DDFF_max_Delta_MJD is not None:
                        if self.verbose > 2:
                            print('# Testing DD pair with imID=%d and %d' % (DD.t['D1imID'][d], DD.t['D2imID'][d]))
                        DD_MJDmin = np.amin(np.array((self.imtable.t['MJD'][DD.t['D1index'][d]],
                                                      self.imtable.t['MJD'][DD.t['D2index'][d]])))
                        DD_MJDmax = np.amax(np.array((self.imtable.t['MJD'][DD.t['D1index'][d]],
                                                      self.imtable.t['MJD'][DD.t['D2index'][d]])))
                        print('DD MJD range', DD_MJDmin, DD_MJDmax)
                        dMJD = np.fabs([DD_MJDmax-FF_MJDmin, DD_MJDmin-FF_MJDmax, DD_MJDmax-FF_MJDmax,
                                        DD_MJDmin-FF_MJDmin])
                        max_dMJD = np.amax(dMJD)
                        if max_dMJD > DDFF_max_Delta_MJD:
                            print('DD pair with imID=%d and %d cannot be used, dMJD=%f>%f' % (DD.t['D1imID'][d],
                                                                                              DD.t['D2imID'][d],
                                                                                              max_dMJD, DDFF_max_Delta_MJD))
                            continue
                    ddcountmin = ddcount[d]
                    d_best = d

            if d_best is None:
                print('SKIPPING following FF pair since there is no matching DD pair!')
                print(FF.t[f])
            else:
                if self.verbose > 2:
                    print('SUCCESS! DD pair with imID %d and %d found for FF pair with imID %d and %d' % (DD.t['D1imID'][d_best],
                                                                                                          DD.t['D2imID'][d_best],
                                                                                                          FF.t['F1imID'][f],
                                                                                                          FF.t['F2imID'][f]))
                ddcount[d_best] += 1
                DDFF.t.add_row({'D1index': DD.t['D1index'][d_best],
                                'D2index': DD.t['D2index'][d_best],
                                'D1imID': DD.t['D1imID'][d_best],
                                'D2imID': DD.t['D2imID'][d_best],
                                'F1index': FF.t['F1index'][f],
                                'F2index': FF.t['F2index'][f],
                                'F1imID': FF.t['F1imID'][f],
                                'F2imID': FF.t['F2imID'][f]})

        print(DDFF.t)
        return(DDFF, ('D1', 'D2', 'F1', 'F2'))

    def get_inputimage_sets(self, reflabel, detector, DD_max_Delta_MJD=None, FF_max_Delta_MJD=None,
                            DDFF_max_Delta_MJD=None):
        imtypes = self.cfg.params[reflabel]['imtypes']
        print("imtypes is {}".format(imtypes))
        imagesets = []
        if imtypes == 'D':
            imagesets, imagelabels = self.getDlist(detector)
        elif imtypes == 'DD':
            imagesets, imagelabels = self.getDDlist(detector, max_Delta_MJD=DD_max_Delta_MJD)
        elif imtypes == 'FF':
            imagesets, imagelabels = self.getFFlist(detector, max_Delta_MJD=FF_max_Delta_MJD)
        elif imtypes == 'DDFF' or imtypes == 'FFDD':
            imagesets, imagelabels = self.getDDFFlist(detector, DD_max_Delta_MJD=DD_max_Delta_MJD,
                                                      FF_max_Delta_MJD=FF_max_Delta_MJD,
                                                      DDFF_max_Delta_MJD=DDFF_max_Delta_MJD)
        else:
            raise RuntimeError("ERROR: imtypes=%s not yet implemented!" % imtypes)
        return(imagesets, imagelabels)

    def mk_ssb_cmds(self,force_redo_ssb=False):
        '''
        Construct the reflabel commands, get the input files
        '''
        if self.verbose:
            print('\n##################################\n### Constructing commands\n##################################')

        # Expand ssbstep list if a shorthand is used
        for reflabel in self.reflabellist:
            step_name = self.cfg.params[reflabel]['ssbsteps'].strip()
            if step_name[-1] == '-':
                self.cfg.params[reflabel]['ssbsteps'] = pipeline_steps.step_minus(step_name, self.cfg.params['instrument'])
            elif step_name[-1] == '+':
                self.cfg.params[reflabel]['ssbsteps'] = pipeline_steps.step_plus(step_name, self.cfg.params['instrument'])

        # this is just to get the correct dtype for the reflabel  columns
        dummyreflabel = astrotableclass()
        dummyreflabel.t['reflabel'] = self.reflabellist
        dummyreflabel.t['ssbsteps'] = [self.cfg.params[reflabel]['ssbsteps'] for reflabel in self.reflabellist]

        self.cmdtable = astrotableclass(names=('reflabel', 'detector', 'cmdID', 'Nim'),
                                        dtype=(dummyreflabel.t['reflabel'].dtype,
                                               self.imtable.t['DETECTOR'].dtype, 'i8', 'i8'))
        #self.cmdtable = astrotableclass(names=('reflabel','detector','cmdID','Nim'))
        self.cmdtable.t['cmdID', 'Nim'].format = '%5d'
        self.cmdtable.t['reflabel', 'detector'].format = '%s'
        self.inputimagestable = astrotableclass(names=('cmdID', 'reflabel', 'imlabel', 'imtype', 'detector',
                                                       'ssbsteps', 'imindex', 'imID', 'MJD', 'fitsfile'),
                                                dtype=('i8', dummyreflabel.t['reflabel'].dtype, 'S40',
                                                       self.imtable.t['imtype'].dtype,
                                                       self.imtable.t['DETECTOR'].dtype,
                                                       dummyreflabel.t['ssbsteps'].dtype,
                                                       'i8', 'i8', 'f8', self.imtable.t['fitsfile'].dtype))
        self.inputimagestable.t['cmdID', 'imindex', 'imID'].format = '%5d'
        self.inputimagestable.t['MJD'].format = '%.8f'

        cmdID = 0
        for reflabel in self.reflabellist:

            print('SSB steps for reflabel %s: %s' % (reflabel, self.cfg.params[reflabel]['ssbsteps']))
            #continue

            counter = 0
            for detector in self.detectors:
                if self.verbose:
                    print('\n#############################\n### Constructing %s commands for detector %s' % (reflabel, detector))
                inputimagesets, inputimagelabels = self.get_inputimage_sets(reflabel, detector,
                                                                            DD_max_Delta_MJD=self.cfg.params['DD']['max_Delta_MJD'],
                                                                            FF_max_Delta_MJD=self.cfg.params['FF']['max_Delta_MJD'],
                                                                            DDFF_max_Delta_MJD=self.cfg.params['DDFF']['max_Delta_MJD'])

                if self.verbose:
                    print('### %d image sets for %s for detector %s' % (len(inputimagesets.t), reflabel, detector))
                if self.verbose > 1:
                    print(inputimagesets.t)
                for i in range(len(inputimagesets.t)):
                    if self.verbose > 1:
                        print('image set index %d' % i)
                    for inputimagelabel in inputimagelabels:

                        # This is the index to the image in imtable
                        imindex = inputimagesets.t['%sindex' % inputimagelabel][i]

                        # sanity test: imindex and imID need to agree!!
                        if inputimagesets.t['%simID' % inputimagelabel][i] != self.imtable.t['imID'][imindex]:
                            raise RuntimeError(('BUG!!! This should not happen! the imIDs %d and %d need to '
                                                'match!') % (inputimagesets.t['%simID' % inputimagelabel][i],
                                                             self.imtable.t['imID'][imindex]))

                        if self.verbose > 3:
                            print(self.imtable.t[inputimagesets.t['%sindex' % inputimagelabel][i]])

                        dict2add = {'cmdID': cmdID, 'imindex': imindex, 'reflabel': reflabel,
                                    'detector': detector, 'imlabel': inputimagelabel,
                                    'ssbsteps': self.cfg.params[reflabel]['ssbsteps']}
                        for key in ['fitsfile', 'imID', 'imtype', 'MJD']:
                            dict2add[key] = self.imtable.t[key][imindex]
                        self.inputimagestable.t.add_row(dict2add)

                    self.cmdtable.t.add_row({'reflabel': reflabel, 'detector': detector, 'cmdID': cmdID,
                                             'Nim': len(inputimagelabels)})
                    cmdID += 1

        if self.verbose > 1:
            print('\n### COMMANDS:')
            print(self.cmdtable.t)
            print('\n### INPUT FILES:')
            print(self.inputimagestable.t)

        print('**** ADD: additional check if input images are the same!! ****')
        self.inputimagestable.write('%s.inputim.txt' % self.basename, verbose=True, clobber=True)

        mmm = CalibPrep(self.cfg.params['instrument'])
        mmm.inputs = self.inputimagestable.t

        print('Bryan: we need to look for ssb files in the ssb output dir, and then also in the optional pipeline_prod_search_dir. Let me know how I should pass this inof!')
        mmm.search_dir = self.ssbdir
        # add additional search dirs!
        if self.cfg.params['output']['pipeline_prod_search_dir'] is not None:
            mmm.search_dir += ',%s' % (self.cfg.params['output']['pipeline_prod_search_dir'])

        mmm.output_dir = self.ssbdir
        mmm.prepare()
        print('Table column names:')
        print(mmm.proc_table.colnames)

        #print('BACK IN MKREFS:')
        #print(mmm.proc_table['index', 'cmdID', 'reflabel', 'output_name', 'steps_to_run', 'repeat_of_index_number', 'index_contained_within'])
        #print(mmm.proc_table['strun_command'][-1])
        #sys.exit()

        self.ssbcmdtable = astrotableclass()
        self.ssbcmdtable.t = mmm.proc_table

        self.check_if_inputfiles_exists()

        self.ssbcmdtable.write('%s.ssbcmds.txt' % self.basename,verbose=True,clobber=True)

        self.cmdtable.write('%s.refcmds.txt' % self.basename,verbose=True,clobber=True)

        sys.exit(0)
        
        #print('strun commands:')
        #print(mmm.strun)

        return(0)

    def check_if_inputfiles_exists(self,outcol='file_exists'):
        print('XXXX %d' % len(self.ssbcmdtable.t))
        file_exists = np.full((len(self.ssbcmdtable.t)),False)
        print(file_exists)
        sys.exit(0)
        
    def run_ssb_cmds(self,batchmode=False):
        print("### run ssb commands: NOT YET IMPLEMENTED!!!")
        return(0)

    def run_ref_cmds(self,batchmode=False):
        print("### run ref commands: NOT YET IMPLEMENTED!!!")
        return(0)

    def submitbatch(self):
        print("### submitbatch: NOT YET IMPLEMENTED!!!")
        sys.exit(0)

    def mk_ref_cmds(self):
        for i in range(len(self.cmdtable.t)):

            print('### cmd ID %d: building cmd for %s' % (self.cmdtable.t['cmdID'][i],
                                                          self.cmdtable.t['reflabel'][i]))
            # (1) get teh input images from self.inputimagestable
            # (2) parse through the options

    def combinerefs(self):
        print("### combinerefs: NOT YET IMPLEMENTED!!!")
        sys.exit(0)

    def overview(self):
        print("### overview: NOT YET IMPLEMENTED!!!")
        sys.exit(0)


if __name__ == '__main__':

    mkrefs = mkrefsclass()
    parser = mkrefs.define_options()
    args = parser.parse_args()

    print("Input files:")
    print(args.reflabels_and_imagelist)

    #print(args)
    #for b in vars(args):
    #    print(b,vars(args)[b])
    #sys.exit()

    # set verbose, debug, and onlyshow level
    mkrefs.verbose = args.verbose
    mkrefs.debug = args.debug
    mkrefs.onlyshow = args.onlyshow

    # Load config files
    mkrefs.loadcfgfiles(args.cfgfile,
                        extracfgfiles=args.extracfgfile,
                        params=args.params,
                        params4all=args.pall,
                        params4sections=args.pp)

    # set the basedir, basename, ssbdir, outsubdir, runID, runlabel
    mkrefs.set_dirs(outrootdir=args.outrootdir, outsubdir=args.outsubdir,
                    runID=args.runID, runIDNdigits=args.runIDNdigits, newrunID=args.newrunID,
                    ssbdir=args.ssbdir, skip_runID_ssbdir=args.skip_runID_ssbdir)

    # get the inputfile list and reflabel list. For input files, get into!
    mkrefs.organize_inputfiles(args.reflabels_and_imagelist)

    # check if current input files are consistent with previous ionput files of same output directory!
    mkrefs.check_inputfiles()
    #sys.exit(0)

    # create the ssb commands
    mkrefs.mk_ssb_cmds()

    # run the ssb commands
    mkrefs.run_ssb_cmds(batchmode=args.batchmode)

    # create the reference file commands
    mkrefs.mk_ref_cmds()

    # run the reference file  commands
    mkrefs.run_ref_cmds(batchmode=args.batchmode)

    mkrefs.combinerefs()

    mkrefs.overview()


#!/usr/bin/env python

'''
Given a pair of input flat field integrations,
calculate the gain map and associated error array

This module uses Armin's delta PTC method
and can be run on raw (uncalibrated) integrations
'''


import sys, os
import re
import glob
import copy
import argparse
import math
import numpy as np
import scipy
from scipy import ndimage
from astropy.io import fits
from .refpixcorr import refpixcorrclass,frameresetclass
from ....tools.math.sigmacut import calcaverageclass
from ....tools.misc.texttable import txttableclass
from ....tools.misc.tools import rmfile


lowedir = '/grp/jwst/wit/nircam/hilbert/epoxy_thinned_area_definition/'
lowepoxy = {}
lowepoxy['NRCA1'] = os.path.join(lowedir,'A1_CV3_low_epoxy_region_DMSorientation.fits')
lowepoxy['NRCA3'] = os.path.join(lowedir,'A3_CV3_low_epoxy_region_DMSorientation.fits')
lowepoxy['NRCA4'] = os.path.join(lowedir,'A4_CV3_low_epoxy_region_DMSorientation.fits')
lowepoxy['NRCALONG'] = os.path.join(lowedir,'ALONG_CV3_low_epoxy_region_DMSorientation.fits')
lowepoxy['NRCB1'] = os.path.join(lowedir,'B1_CV3_low_epoxy_region_DMSorientation.fits')
lowepoxy['NRCB3'] = os.path.join(lowedir,'B3_CV3_low_epoxy_region_DMSorientation.fits')
lowepoxy['NRCB4'] = os.path.join(lowedir,'B4_CV3_low_epoxy_region_DMSorientation.fits')


class ExposureGain:
    def __init__(self, flatfile1, flatfile2, darkfile1=None, darkfile2=None, readnoise_reffile=None,
                 agressive_masking=True, gain_per_amp=True, boxsize=64, output_dir='./', verbos=True):
        self.flatfile1 = flatfile1
        self.flatfile2 = flatfile2
        #self.Nx = None set this after reading in the flats
        #self.Ny = None set this after reading in the flats

        if readnoise_reffile is None and (darkfile1 is None or darkfile2 is None):
            raise ValueError('Either readnoise_reffile or darkfile1 and darkfile2 must be specified.')

        if readnoise_reffile is not None and (darkfile1 is not None or darkfile2 is not None):
            print("Both readnoise reference file and input dark files supplied. The reference file takes precendence.")

        if readnoise_reffile.lower() == 'crds':
            Retrieve the reference file from crds


        self.verbose = verbose
        self.aggressivemasking = agressive_masking
        #self.test = False
        self.gain_per_amp = gain_per_amp  # force a single gain value per amplifier. overrides boxsize.
        self.boxsize = boxsize
        self.darks = (darkfile1, darkfile2)
        self.readnoise_reffile = readnoise_reffile # If None, use the code to calculate readnoise. If a filename, read in and use
        self.gmax4dark = None
        self.xmin = 0
        self.xmax = 2040
        self.ymin = 0
        self.ymax = 2040
        self.Ngroupskeyword = 'NGROUPS'
        self.gmin = 0
        self.gmax = None
        self.Smeanmax4fit = 20000.
        self.skipcorrect4framereset = False
        self.refpixboxsize = 12
        self.skiprefpixcorr = False

        # array list with readnoise info
        self.rdnoise={}
        # set the extensions that are required to None
        for e in ['im','err','flag']:
            self.rdnoise[e]=None

        self.readnoisemethod=None


    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('flatfile1', help='Name of first flat integration to use',default=None)
        parser.add_option('flatfile2', help='Name of second flat integration to use',default=None)
        parser.add_option('-v', '--verbose', action='count', dest='verbose',default=0)
        parser.add_option('--aggressivemasking', default=False, action="store_true",
                          help='make aggressive bpm masks, up to 20-30% of masking')
        parser.add_option('-t','--test', default=False, action="store_true",
                          help='test')
        parser.add_option('--outsubdir'  , default=None , type="string",
                          help='Add this subdir to path (default=%default)')
        parser.add_option('-b','--boxsize'  , default=64 , type="int",
                          help='boxsize for stats (default=%default)')
        parser.add_option('--darks', default=(None,None), nargs=2, type="string",
                          help='specify the two darks to be used for the readnoise determination')
        parser.add_option('-a','--autodark4readnoise', default=False, action="store_true",
                          help='search for two images in workspace subdir with *DARK* fits and determine the readnoise from these images.')
        parser.add_option('--gmax4dark'  , default=None , type="int",
                          help='maximum group for dark frame (default=%default)')
        parser.add_option('--xmin'  , default=0 , type="int",
                          help='xmin for stats (default=%default)')
        parser.add_option('--xmax'  , default=2048 , type="int",
                          help='xmax for stats (default=%default)')
        parser.add_option('--ymin'  , default=0 , type="int",
                          help='ymin for stats (default=%default)')
        parser.add_option('--ymax'  , default=2048 , type="int",
                          help='ymax for stats (default=%default)')
        parser.add_option('--Ngroupskeyword'  , default='NGROUPS' , type="string",
                          help='keyword for # of groups (default=%default)')
#        parser.add_option('-r','--readnoise'  , default=7.084/math.sqrt(2) , type="float",
#                          help='readnoise (default=%default)')
        parser.add_option('--gmin'  , default=0 , type="int",
                          help='minimum group (default=%default)')
        parser.add_option('--gmax'  , default=None , type="int",
                          help='maximum group (default=%default)')
        parser.add_option('--Smeanmax4fit'  , default=20000.0 , type="float",
                          help='all groups with a mean signal>Smeanmax4fit are not used for stdev and gain linear fit (default=%default)')
        parser.add_option('--skipcorrect4framereset' , default=False, action="store_true",
                          help='do not correct for frame reset (default=%default)')
        parser.add_option('--refpixboxsize'  , default=12, type="int",
                          help='boxsize for boxcar smoothing of biasdrift vector (default=%default)')
        parser.add_option('--skiprefpixcorr' , default=False, action="store_true",
                          help='do not correct 1/f with side refpix (default=%default)')
        return(parser)


    def calcgain(self,boxmap,box_location,tsub,imtype,Smean_max=None):

        gainim = scipy.zeros((self.Ny,self.Nx),dtype=float)
        gainim_err = scipy.zeros((self.Ny,self.Nx),dtype=float)

        keys_all  =  tsub.CUT_inrange('S_mean',None,Smean_max)

        for boxval in box_location:
            x, y = box_location[boxval]
            #x=xinfo[0]
            #while x<xinfo[1] and x+xinfo[2]<=self.Nx:
            keys_xcut  = tsub.CUT_inrange('xmin',x,x,keys=keys_all)
            #y=yinfo[0]
            #while y<yinfo[1] and y+yinfo[2]<=self.Ny:
            keys  = tsub.CUT_inrange('ymin',y,y,keys=keys_xcut)
            keys  = tsub.CUT_none_vals('gain',keys=keys)

            # fits a slope through stdev2 versus g for dsub2
            (slope_b,slope_b_err,offset_a,offset_a_err) = tsub.straightline(keys,'geff','gain',dyval=None)
            if self.verbose:
                print('%4d,%4d: %s fit to gain versus geff: slope=%.4f(%.4f) offset=%.4f(%.4f)' % (x,y,imtype,slope_b,slope_b_err,offset_a,offset_a_err))

            # the true gain at t=0 is at g=-1, since at g=0, already delta t time has passed
            gain = -1.0*slope_b+offset_a
            gain_err = math.sqrt((-1)*(-1)*slope_b_err*slope_b_err + offset_a_err*offset_a_err)

            goodpix = boxmap == boxval
            gainim[goodpix] = gain
            gainim_err[goodpix] = gain_err

            #y+=yinfo[2]
            #x+=xinfo[2]

        gainfilename = self.outfilebasename+'.gain.%s.fits' % imtype
        phdu = fits.PrimaryHDU(data=None)
        hdu_gain = fits.ImageHDU(gainim,name='gain')
        phdu = self.set_header_info(phdu)
        hdu_gain_err = fits.ImageHDU(gainim_err,name='gain_err')
        hdulist = fits.HDUList([phdu,hdu_gain,hdu_gain_err])
        print('Saving {}'.format(gainfilename))
        rmfile(gainfilename)
        hdulist.writeto(gainfilename)

        del gainim, gainim_err, hdulist, hdu_gain, hdu_gain_err


    def set_header_info(self,hdu):
        '''
        Populate necessary header information in the output file
        '''
        hdu.header['FILETYPE'] = 'GAIN'
        hdu.header['DETECTOR'] = self.hdr1['DETECTOR']
        hdu.header['INSTRUME'] = 'NIRCAM'
        hdu.header['TELESCOP'] = 'JWST'
        hdu.header['FLATFIL1'] = self.flatfile1
        hdu.header['FLATFIL2'] = self.flatfile2
        hdu.header['DARKFIL1'] = self.darks[0]
        hdu.header['DARKFIL2'] = self.darks[1]
        return hdu


    def correct_stdev2_with_rdnoiseterms(self,boxmap,box_location,tsubik,tsubg0,tdsub1,
                                         tdsub2,tdsubij,Smean_max=None,Pmax=20.0,Nmin=3):

        print('### Correcting for readnoise')

        if self.readnoisemethod=='dsub12':
            if self.verbose:
                print(('Determining readnoise from dsub1 and dsub2 for correction'
                       ', readnoisemethod="dsub12"!!'))
        else:
            if self.verbose:
                print('readnoisemethod="{}"'.format(self.readnoisemethod))

        # arrays for the readnoise determined from dsub1 and dsub2.
        # This is only to save the readnoise image
        rdnoiseim = scipy.zeros((self.Ny,self.Nx),dtype=float)
        rdnoiseim_err = scipy.zeros((self.Ny,self.Nx),dtype=float)
        rdnoiseim_Pskip = scipy.zeros((self.Ny,self.Nx),dtype=float)
        rdnoiseim_Nused = scipy.zeros((self.Ny,self.Nx),dtype=np.int32)
        rdnoiseim_flag = scipy.zeros((self.Ny,self.Nx),dtype=np.int16)

        sigmacut = calcaverageclass()
        keyssubik_all   =  tsubik.CUT_inrange('S_mean',None,Smean_max)
        keyssubg0_all =  tsubg0.CUT_inrange('S_mean',None,Smean_max)
        keyssub1_all  =  tdsub1.CUT_inrange('S_mean',None,Smean_max)
        keyssub2_all  =  tdsub2.CUT_inrange('S_mean',None,Smean_max)
        keyssubij_all  =  tdsubij.CUT_inrange('S_mean',None,Smean_max)

        for boxval in box_location:
            x, y = box_location[boxval]
            goodpix = boxmap == boxval

        #x=xinfo[0]
        #while x<xinfo[1] and x+xinfo[2]<=self.Nx:
            keyssub1_xcut  = tdsub1.CUT_inrange('xmin',x,x,keys=keyssub1_all)
            keyssub2_xcut  = tdsub2.CUT_inrange('xmin',x,x,keys=keyssub2_all)
            keyssubij_xcut  = tdsubij.CUT_inrange('xmin',x,x,keys=keyssubij_all)
            keyssubg0_xcut = tsubg0.CUT_inrange('xmin',x,x,keys=keyssubg0_all)
            keyssubik_xcut   = tsubik.CUT_inrange('xmin',x,x,keys=keyssubik_all)
            #y=yinfo[0]
            #while y<yinfo[1] and y+yinfo[2]<=self.Ny:
            keyssub1  = tdsub1.CUT_inrange('ymin',y,y,keys=keyssub1_xcut)
            keyssub2  = tdsub2.CUT_inrange('ymin',y,y,keys=keyssub2_xcut)
            keyssubij  = tdsubij.CUT_inrange('ymin',y,y,keys=keyssubij_xcut)
            keyssubg0 = tsubg0.CUT_inrange('ymin',y,y,keys=keyssubg0_xcut)
            keyssubik   = tsubik.CUT_inrange('ymin',y,y,keys=keyssubik_xcut)

            tdsub2.printtxttable(keys=keyssub2,cols=['geff','stdev2'])

            #print(x,y)
            #print(keyssub2_xcut)
            #print(keyssub2)

            # fits a slope through stdev2 versus g for dsub2
            if len(keyssub2) < 2:
                print("keyssub2 is too short for boxval {}".format(boxval))
                print(keyssub2)
                sys.exit()
            else:
                (slope_b,slope_b_err,offset_a,offset_a_err) = tdsub2.straightline(keyssub2,'geff','stdev2',dyval=None)
            if self.verbose:
                print('%4d,%4d: dsub1 fit to stdev2 versus geff: slope=%.4f(%.4f) offset=%.4f(%.4f)' % (x,y,slope_b,slope_b_err,offset_a,offset_a_err))

            # subtract stdev2 from dsub2 slope fit from stdev2 of dsub1.
            # This difference should be 2*readnoise^2
            stdev2array = scipy.zeros((len(keyssub1)),dtype=float)
            counter=0
            for key1 in keyssub1:
                stdev2diff = tdsub1.getentry(key1,'stdev2')-(tdsub1.getentry(key1,'geff')*slope_b+offset_a)
                stdev2array[counter]=stdev2diff
                counter+=1
            #calculate the readnoise^2 with dsub12 method
            sigmacut.calcaverage_sigmacutloop(stdev2array,Nsigma=3.0,verbose=0)
            if sigmacut.mean!=None and sigmacut.mean>0.0:
                rdnoise2_dsub12 = sigmacut.mean*0.5
                rdnoise_dsub12 = math.sqrt(rdnoise2_dsub12)
                rdnoise_dsub12_err = math.sqrt(0.5/(4.0*sigmacut.mean)*sigmacut.mean_err*sigmacut.mean_err)
                #rdnoiseim[y:y+yinfo[2],x:x+xinfo[2]]=rdnoise_dsub12
                #rdnoiseim_err[y:y+yinfo[2],x:x+xinfo[2]]=rdnoise_dsub12_err
                rdnoiseim[goodpix]=rdnoise_dsub12
                rdnoiseim_err[goodpix]=rdnoise_dsub12_err
                Pskip=100.0*sigmacut.Nskipped/(sigmacut.Nused+sigmacut.Nskipped)
                #rdnoiseim_Nused[y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.Nused
                #rdnoiseim_Pskip[y:y+yinfo[2],x:x+xinfo[2]]=Pskip
                rdnoiseim_Nused[goodpix]=sigmacut.Nused
                rdnoiseim_Pskip[goodpix]=Pskip
                if  Pskip<=Pmax and sigmacut.Nused>=Nmin:
                    #rdnoiseim_flag[y:y+yinfo[2],x:x+xinfo[2]]=0
                    rdnoiseim_flag[goodpix]=0
                else:
                    #rdnoiseim_flag[y:y+yinfo[2],x:x+xinfo[2]]=1
                    rdnoiseim_flag[goodpix]=1
            else:
                rdnoise_dsub12 = np.nan

            # Which readnoise to take for correction?
            if self.rdnoise['im'] is None:
                rdnoise = rdnoise_dsub12
            else:
                #rdnoise = self.rdnoise['im'][int(y+0.5*yinfo[2]),int(x+0.5*xinfo[2])]
                #rdnoise = self.rdnoise['im'][int(y+0.5*self.boxsize),int(x+0.5*self.boxsize)]
                rdnoise = self.rdnoise['im'][y,x]
            rdnoise2=rdnoise*rdnoise

            if self.verbose:
                print('rdnoise=%.3f' % (rdnoise))
                if rdnoise == 0:
                    if self.rdnoise['im'] is None:
                        print('rdnoise_dsub12',rdnoise_dsub12)
                    else:
                        print(y,self.boxsize,y+0.5*self.boxsize,x,self.boxsize,x+0.5*self.boxsize,self.rdnoise['im'][int(y+0.5*self.boxsize),int(x+0.5*self.boxsize)])
                    print("WARNING: readnoise is zero here!!")
                    print("There is a problem. Quitting.")
                    sys.exit()


            # get the correction for the default diffim stdev2
            keyssubik_g0 = tsubik.CUT_inrange('g',0,0,keys=keyssubik)
            if len(keyssubik_g0)!=1:
                print('ERROR keyssubik_g0!',keyssubik_g0)
                tsubik.printtxttable(keys=keyssubik)
                print('g0')
                tsubik.printtxttable(keys=keyssubik_g0)
                sys.exit(0)
            correction4subik = tsubik.getentry(keyssubik_g0[0],'stdev2')

            # correct all the tables with the readnoise^2
            for key in keyssubik:
                tsubik.setentry(key,'stdev2_rdnoisecorr',tsubik.getentry(key,'stdev2')-correction4subik)
            for key in keyssubg0:
                tsubg0.setentry(key,'stdev2_rdnoisecorr',tsubg0.getentry(key,'stdev2')-4.0*rdnoise2)
            for key in keyssub1:
                tdsub1.setentry(key,'stdev2_rdnoisecorr',tdsub1.getentry(key,'stdev2')-6.0*rdnoise2)
            for key in keyssub2:
                tdsub2.setentry(key,'stdev2_rdnoisecorr',tdsub2.getentry(key,'stdev2')-4.0*rdnoise2)
            for key in keyssubij:
                tdsubij.setentry(key,'stdev2_rdnoisecorr',tdsubij.getentry(key,'stdev2')-4.0*rdnoise2)

            #y+=yinfo[2]
            #x+=xinfo[2]

        rdnoisefilename = self.outfilebasename+'.rdnoise_dsub.fits'
        readnoisehdulist = fits.HDUList([fits.PrimaryHDU(data=None),
                                         fits.ImageHDU(rdnoiseim,name='rdnoise'),
                                         fits.ImageHDU(rdnoiseim_err,name='rdnoiseerror'),
                                         fits.ImageHDU(rdnoiseim_flag,name='flag'),
                                         fits.ImageHDU(rdnoiseim_Nused,name='Nused'),
                                         fits.ImageHDU(rdnoiseim_Pskip,name='Pskip')])
        print('Saving {}'.format(rdnoisefilename))
        rmfile(rdnoisefilename)
        readnoisehdulist.writeto(rdnoisefilename)


    def nanunique(self,x):
        '''Return the unique non-NaN elements of an array'''
        a = np.unique(x)
        r = []
        for i in a:
            if np.isnan(i):
            #if i in r or (np.isnan(i) and np.any(np.isnan(r))):
                continue
            else:
                r.append(i)
        return np.array(r)


    def calcgain_grid(self,diffim,t,g,boxmap,box_locations,mask=None,S=None,
                      dS=None,outfile=None,fixmean=None,imtype=None):
        keys=[]
        sigmacut = calcaverageclass()
        #x=xinfo[0]

        print('### calculating gain for imtype {}'.format(imtype))
        #while x<xinfo[1] and x+xinfo[2]<=self.Nx:
        #    y=yinfo[0]
        #    while y<yinfo[1] and y+yinfo[2]<=self.Ny:
        for boxnum in self.nanunique(boxmap):
            #print('boxnum in calcgain_grid: {}'.format(boxnum))
            pixels = boxmap == boxnum

            # Get the central location of the box. These
            # will go into the text table
            x, y = box_locations[boxnum]

            mask2use = None
            if not (mask is None):
                #mask2use = mask[y:y+yinfo[2],x:x+xinfo[2]]
                mask2use = mask[pixels]

            #print 'x:%4d %4d %4d y:%4d %4d %4d ' % (x,x+xinfo[2],xinfo[2],y,y+yinfo[2],yinfo[2])
            #sigmacut.calcaverage_sigmacutloop(diffim[y:y+yinfo[2],x:x+xinfo[2]],mask=mask2use,Nsigma=3.0,saveused=True,verbose=0,fixmean=fixmean)
            sigmacut.calcaverage_sigmacutloop(diffim[pixels],mask=mask2use,Nsigma=3.0,saveused=True,verbose=0,fixmean=fixmean)
            key = sigmacut.results2texttable(t)
            keys.append(key)

            # For S and dS calculation: use exactly same pixels used for the stdev calculations!
            if not (mask is None):
                mask4S = np.logical_or(mask2use,sigmacut.clipped)
            else:
                mask4S = sigmacut.clipped

            if imtype!=None:
                stdev2 = t.getentry(key,'stdev')*t.getentry(key,'stdev')
                t.setentry(key,'stdev2',stdev2)

            if imtype in ['dsub2','dsubij']:
                t.add2row(key,{'g':g,'geff':g-0.5,'xmin':x,'xdelta':self.boxsize,'ymin':y,'ydelta':self.boxsize})
            else:
                t.add2row(key,{'g':g,'geff':g,'xmin':x,'xdelta':self.boxsize,'ymin':y,'ydelta':self.boxsize})
            if not (S is None):
                #sigmacut.calcaverage_sigmacutloop(S[y:y+yinfo[2],x:x+xinfo[2]],mask=mask2use,Nsigma=3.0,verbose=0)
                #sigmacut.calcaverage_sigmacutloop(S[y:y+yinfo[2],x:x+xinfo[2]],mask=mask4S,Nsigma=0.0,verbose=0)
                sigmacut.calcaverage_sigmacutloop(S[pixels],mask=mask4S,Nsigma=0.0,verbose=0)
                sigmacut.results2texttable(t,key=key,meancol='S_mean',meanerrcol='S_mean_err',stdevcol='S_stdev',Nusedcol='S_Nused',Nskippedcol='S_Nskipped',Pskippedcol='S_Pskipped',convergedcol='S_converged',iterationcol='S_i')
            if not (dS is None):
                #sigmacut.calcaverage_sigmacutloop(dS[y:y+yinfo[2],x:x+xinfo[2]],mask=mask2use,Nsigma=3.0,verbose=0)
                #sigmacut.calcaverage_sigmacutloop(dS[y:y+yinfo[2],x:x+xinfo[2]],mask=mask4S,Nsigma=0.0,verbose=0)
                sigmacut.calcaverage_sigmacutloop(dS[pixels],mask=mask4S,Nsigma=0.0,verbose=0)
                sigmacut.results2texttable(t,key=key,meancol='dS_mean',meanerrcol='dS_mean_err',stdevcol='dS_stdev')

                #y+=yinfo[2]
            #x+=xinfo[2]

        if outfile != None:
            t.save2file(outfile,keys=keys,verbose=True)
        return(keys)


    def make_box_map(self,boxsize,detector):
        '''Create and return a map of boxes across the
           detector. A separate gain value will be
           calculated within each box'''

        #numboxesx = np.int(510 / boxsize)
        #numboxesy = np.int(2040 / boxsize)
        numboxesx = np.int(512 / boxsize)
        numboxesy = np.int(2048 / boxsize)
        boxesperquad = numboxesx*numboxesy
        xstart = np.arange(numboxesx+1) * boxsize
        ystart = np.arange(numboxesy+1) * boxsize
        #quadstart = np.array([0,508,1020,1532,2044])
        quadstart = np.array([0,512,1024,1536,2044])

        # create a map of the detector, and label the pixels
        # that go with each box
        #map = np.zeros((2040,2040))
        map = np.zeros((2048,2048))
        for q,quad in enumerate(quadstart[:-1]):
            for i,xs in enumerate(xstart[:-1]):
                xend = xstart[i+1]
                if i == (numboxesx-1):
                    if quadstart[q+1] > xend:        # <? was > before!!!!
                        xend = quadstart[q+1]

                for j,ys in enumerate(ystart[:-1]):
                    yend = ystart[j+1]
                    if j == (numboxesy-1):
                        yend = 2044
                    map[ys:yend,quad+xs:quad+xend] = q*10000 + (j*numboxesx+i)

        # now we need to deal with the low-epoxy region, if any
        if detector in ['NRCA1','NRCA3','NRCA4','NRCALONG','NRCB1','NRCB3','NRCB4']:
            if self.verbose:
                print("Reading in low epoxy map")
            efile = lowepoxy[detector]
            ehdu = fits.open(efile)
            epoxymask = ehdu[1].data
            #epoxymask = epoxymask[4:2044,4:2044]
            ehdu.close()

            previousmax = np.max(map)
            #if os.path.isfile('post_epoxymask_map.fits'):
            #    os.remove('post_epoxymask_map.fits')
            map = map + (epoxymask*50000)
            #pphdu=fits.PrimaryHDU(map)
            #pphdu.writeto('post_epoxymask_map.fits')

            # Make a copy of map because we need to refer to
            # the original map even as we're making changes
            newmap = copy.deepcopy(map)

            # Clean up. Since the epoxy masks have very irregular boundaries,
            # we need to deal with areas where a given box contains only a few pixels
            # Let's identify these boxes and then move their pixels into an appropriate
            # neighboring box.

            # Loop over boxes and check how many pixels are in each one
            # if you find a box that has less than 50% of the nominally expected
            # number of pixels, then shift those pixels into a neighboring box.
            expected = boxsize * boxsize
            for boxval in self.nanunique(map):
                numpix = len(np.where(newmap == boxval)[0])
                if self.verbose:
                    print("BOX NUMBER: {}, PIXFRAC: {}".format(boxval,1.*numpix/(1.*expected)))
                #if (1.*numpix)/(1.*expected) >= 0.5:
                    #print("BOX OK: Box value: {}, number of pixels in box: {}, expected: {}".format(boxval,numpix,expected))
                if (1.*numpix)/(1.*expected) < 0.5:
                    if self.verbose:
                        print("Box value: {}, number of pixels in box: {}, expected: {}".format(boxval,numpix,expected))
                    orig = newmap == boxval
                    plusone = newmap == (boxval+1)
                    minusone = newmap == (boxval-1)
                    plus4 = newmap == (boxval+numboxesx)
                    minus4 = newmap == (boxval-numboxesx)
                    if numboxesx != 1:
                        plus3 = newmap == (boxval+numboxesx-1)
                        minus3 = newmap == (boxval-(numboxesx-1))
                        plus5 = newmap == (boxval+numboxesx+1)
                        minus5 = newmap == (boxval-(numboxesx+1))
                    else:
                        plus3 = plusone
                        minus3 = minusone
                        plus5 = plusone
                        minus5 = minusone

                    # Center of mass of box in question
                    ycenter,xcenter = ndimage.measurements.center_of_mass(orig.astype(int))

                    # Center of mass of surrounding boxes
                    yp1center,xp1center = ndimage.measurements.center_of_mass(plusone.astype(int))
                    ym1center,xm1center = ndimage.measurements.center_of_mass(minusone.astype(int))
                    yp4center,xp4center = ndimage.measurements.center_of_mass(plus4.astype(int))
                    ym4center,xm4center = ndimage.measurements.center_of_mass(minus4.astype(int))
                    yp3center,xp3center = ndimage.measurements.center_of_mass(plus3.astype(int))
                    ym3center,xm3center = ndimage.measurements.center_of_mass(minus3.astype(int))
                    yp5center,xp5center = ndimage.measurements.center_of_mass(plus5.astype(int))
                    ym5center,xm5center = ndimage.measurements.center_of_mass(minus5.astype(int))

                    # Distances
                    dp1 = np.sqrt((ycenter-yp1center)*(ycenter-yp1center) +
                                  (xcenter-xp1center)*(xcenter-xp1center))
                    dm1 = np.sqrt((ycenter-ym1center)*(ycenter-ym1center) +
                                  (xcenter-xm1center)*(xcenter-xm1center))
                    dp4 = np.sqrt((ycenter-yp4center)*(ycenter-yp4center) +
                                  (xcenter-xp4center)*(xcenter-xp4center))
                    dm4 = np.sqrt((ycenter-ym4center)*(ycenter-ym4center) +
                                  (xcenter-xm4center)*(xcenter-xm4center))
                    dp3 = np.sqrt((ycenter-yp3center)*(ycenter-yp3center) +
                                  (xcenter-xp3center)*(xcenter-xp3center))
                    dm3 = np.sqrt((ycenter-ym3center)*(ycenter-ym3center) +
                                  (xcenter-xm3center)*(xcenter-xm3center))
                    dp5 = np.sqrt((ycenter-yp5center)*(ycenter-yp5center) +
                                  (xcenter-xp5center)*(xcenter-xp5center))
                    dm5 = np.sqrt((ycenter-ym5center)*(ycenter-ym5center) +
                                  (xcenter-xm5center)*(xcenter-xm5center))
                    distances = np.array([dp1,dm1,dp4,dm4,dp3,dm3,dp5,dm5])
                    if numboxesx != 1:
                        compboxvals = np.array([boxval+1,boxval-1,boxval+numboxesx,boxval-numboxesx,
                                                boxval+numboxesx-1,boxval-(numboxesx-1),boxval+numboxesx+1,
                                                boxval-(numboxesx+1)])
                    else:
                        compboxvals = np.array([boxval+1,boxval-1,boxval+numboxesx,boxval-numboxesx,boxval+1,
                                                boxval-1,boxval+1,boxval-1])
                    #if there are two or fewer boxes across in the x-dimension, then don't
                    #use boxes that are diagonal from the box in question.
                    mindistidx = np.where(distances == np.nanmin(distances))[0]

                    if len(mindistidx) > 0:
                        mindistidx = mindistidx[0]
                        mindist = distances[mindistidx]
                        updatedval = compboxvals[mindistidx]
                        if self.verbose:
                            print('new box value is:',updatedval)
                        if ~np.isnan(mindist):
                            newmap[orig] = updatedval

                            #case where there are no surrounding boxes within +/- 5
                            #Try searching outward for a different box value
                    else:
                        width = boxsize*2
                        xx = int(xcenter)
                        yy = int(ycenter)
                        surrounding = np.absolute(newmap - boxval)
                        #nx,ny = 2040,2040
                        nx,ny = 2048,2048
                        x, y = np.meshgrid(np.arange(nx)-xcenter, np.arange(ny)-ycenter)
                        surrounding_dist = np.sqrt(x**2 + y**2)
                        # Now find the closest pixel with a value that is not:
                        #boxval in a different quadrant
                        #across the low-epoxy region boundary
                        otherboxes = (surrounding != 0) & (surrounding < 9000)

                        #if os.path.isfile('surround.fits'):
                        #    os.remove('surround.fits')
                        #if os.path.isfile('surround_dist.fits'):
                        #    os.remove('surround_dist.fits')
                        #if os.path.isfile('otherboxes.fits'):
                        #    os.remove('otherboxes.fits')
                        if self.test:
                            thdu=fits.PrimaryHDU(surrounding)
                            thdu.writeto('surround.fits',overwrite=True)
                            t2hdu=fits.PrimaryHDU(surrounding_dist)
                            t2hdu.writeto('surround_dist.fits',overwrite=True)
                            t3hdu=fits.PrimaryHDU(otherboxes.astype(int))
                            t3hdu.writeto('otherboxes.fits',overwrite=True)

                        if len(np.where(otherboxes == True)[0]) > 0:
                            minidist = np.min(surrounding_dist[otherboxes])
                            nearestbox = np.where((otherboxes == True) & (surrounding_dist == minidist))
                            bval = newmap[nearestbox[0][0],nearestbox[1][0]] #+ boxval
                            if self.verbose:
                                print("{}, nearest box value is {}.".format(boxval,bval))
                            newmap[orig] = bval
                        else:
                            if self.verbose:
                                print('For {}, no nearby boxes found!'.format(boxval))
                            # To be here, all boxes in the rest of the quadrant must be bad.
                            # this can be true if we're talking about the edge of the low epoxy
                            # region that is just barely sticking into quad 3, so there aren't
                            # any other low epoxy boxes in quad 3. In this case, we can
                            # keep the pixels in their own small box, or push them across
                            # the low-epoxy boundary...
                            # for B1, the "box" has only 29 pixels...
                            if self.verbose:
                                print('Must be at the edge of a low-epoxy area that just')
                                print('crosses quad boundaries.')
                            if numpix < 5000:
                                if self.verbose:
                                    print('Pushing these {} pixels to adjacent box across low-epoxy region'.format(numpix))
                                newmap[orig] = newmap[orig] - 50000
                            else:
                                if self.verbose:
                                    print('Keeping these {} pixels as their own, small box.'.format(numpix))

        else: #this is for the detectors that have no low-epoxy area
            newmap = map

        if self.verbose:
            print("box map complete")

        #if os.path.isfile('post_epoxymask_fix_map.fits'):
        #    os.remove('post_epoxymask_fix_map.fits')
        #pphdu=fits.PrimaryHDU(newmap)
        #pphdu.writeto('post_epoxymask_fix_map.fits',clobber=True)

        # Now find the center of mass for each box in the map.
        # Return a list of this position for all boxes in
        # addition to the boxmap itself. These will be used
        # for interactions with text tables
        if self.verbose:
            print("creating box_centers")
        # Insert mask into 2048x2048 array, with the refpix set to NaN
        ffmap = np.zeros((2048,2048))
        ffmap[:,:] = np.nan
        #ffmap[4:2044,4:2044] = newmap
        ffmap[4:2044,4:2044] = newmap[4:2044,4:2044]

        box_centers = {}
        for bnum in self.nanunique(ffmap):
            box = ffmap == bnum
            ycent,xcent = ndimage.measurements.center_of_mass(box.astype(int))
            xcent = int(xcent)
            ycent = int(ycent)
            if ffmap[ycent,xcent] == bnum:
                box_centers[bnum] = (int(xcent),int(ycent))
            else:
                goodpix = np.where(ffmap == bnum)
                dist = np.sqrt((xcent-goodpix[1])**2 + (ycent-goodpix[0])**2)
                mindist = np.where(dist == np.min(dist))
                box_centers[bnum] = (goodpix[0][mindist[0][0]],goodpix[1][mindist[0][0]])

        if np.nan in box_centers.keys():
            print('nan in box_centers!!')
            sys.exit()
        return ffmap, box_centers


    def calc_flattenfactor(self,S1,S2,mask=None):
        mask2use = scipy.where(S2<=0.0,True,False)
        if not (mask is None):
            mask2use |= mask
        ratio = S1/S2
        sigmacut = calcaverageclass()
        sigmacut.calcaverage_sigmacutloop(ratio,mask=mask2use,Nsigma=3.0,verbose=0)
        if self.verbose:
            print(sigmacut.__str__())
        return(sigmacut.mean)


    def calc_cum_stdev2(self,t,box_location,gmin,gmax):
        t.configcols(['cum_stdev2_rdnoisecorr','ave_Smean','cum_avestdev2'],'f','%.3f',visible=1)

        for boxval in box_location:
            x, y = box_location[boxval]

        #x=xinfo[0]
        #print 'CCC',xinfo,yinfo
        #while x<xinfo[1] and x+xinfo[2]<=self.Nx:
            #y=yinfo[0]
            #while y<yinfo[1] and y+yinfo[2]<=self.Ny:
            keys = t.CUT_inrange('xmin',x,x)
            keys = t.CUT_inrange('ymin',y,y,keys=keys)
            keys = t.sortkeysbycols(keys,'g',asstring=0)
            sum=0.0
            for key in keys:
                if t.getentry(key,'stdev2_rdnoisecorr')==None:
                    break
                sum+=t.getentry(key,'stdev2_rdnoisecorr')
                t.setentry(key,'cum_stdev2_rdnoisecorr',sum)
            #y+=yinfo[2]
        #x+=xinfo[2]

        for g in range(gmin,gmax):
            keys = t.CUT_inrange('g',g,g)
            keys = t.CUT_none_vals('cum_stdev2_rdnoisecorr',keys=keys)
            if len(keys)>0:
                cum_avestdev2 = t.calcaverage(keys,'cum_stdev2_rdnoisecorr')
                ave_Smean = t.calcaverage(keys,'S_mean')
                for key in keys:
                    t.setentry(key,'ave_Smean',ave_Smean)
                    t.setentry(key,'cum_avestdev2',cum_avestdev2)

        #t.printtxttable(keys=keys,cols=['g','xmin','xdelta','ymin','ydelta','stdev2_rdnoisecorr','cum_stdev2_rdnoisecorr'])

    #def readnoise_sq_4imtype(self,imtype):
    #    if imtype!=None:
    #        rdnoise2=self.options.readnoise*self.options.readnoise
    #        if imtype=='default':rdnoise2_tot = 2.0*rdnoise2
    #        elif imtype=='subg0':rdnoise2_tot = 4.0*rdnoise2
    #        #elif imtype=='subgm1':rdnoise2_tot = 2.0*rdnoise2
    #        elif imtype=='dsub1':rdnoise2_tot = 6.0*rdnoise2
    #        elif imtype=='dsub2':rdnoise2_tot = 4.0*rdnoise2
    #        else: raise RuntimeError,"Wrong imtype %s" % imtype
    #        return(rdnoise2_tot)
    #    else:
    #        return(None)

    #def calc_varying_readnoise(self,t,imtype,g0=2):
    #    t.configcols(['rdnoise2_var'],'f','%.3f',visible=1)
    #    keys = t.CUT_inrange('g',g0,g0)
    #    dS0  = t.calcaverage(keys,'dS_mean')
    #    print 'dS0:',dS0
    #    for key in t.allrowkeys:
    #        if t.getentry(key,'dS_mean')>=10.0:
    #            rdnoise2_tot = self.readnoise_sq_4imtype(imtype)
    #            rdnoise2_var = rdnoise2_tot * (t.getentry(key,'dS_mean')/dS0*t.getentry(key,'dS_mean')/dS0)
    #            stdev2 = t.getentry(key,'stdev2')
    #            t.setentry(key,'rdnoise2_var',rdnoise2_var)
    #            t.setentry(key,'stdev2_rdnoisecorrX',stdev2-rdnoise2_var)

    def calcstdev2_model(self,t,gmin,gmax,g0=2):
        t.configcols(['stdev2_rdnc_model','ave_Smean','ave_stdev2_rdnc_model'],'f','%.3f',visible=1)
        keys  = t.CUT_inrange('g',g0,g0)
        stdev2_g0= t.calcaverage(keys,'stdev2_rdnoisecorr')
        dS0   = t.calcaverage(keys,'dS_mean')
        if self.verbose:
            print('stdev2_g0:',stdev2_g0)
            print('dS0:',dS0)
        for key in t.allrowkeys:
            if t.getentry(key,'dS_mean')>=10.0:
                dSratio = t.getentry(key,'dS_mean')/dS0
                t.setentry(key,'stdev2_rdnc_model',stdev2_g0*dSratio*dSratio)

        for g in range(gmin,gmax):
            keys = t.CUT_inrange('g',g,g)
            keys = t.CUT_none_vals('stdev2_rdnc_model',keys=keys)
            if len(keys)>0:
                ave_stdev2_rdnc_model = t.calcaverage(keys,'stdev2_rdnc_model')
                ave_Smean = t.calcaverage(keys,'S_mean')
                for key in keys:
                    t.setentry(key,'ave_Smean',ave_Smean)
                    t.setentry(key,'ave_stdev2_rdnc_model',ave_stdev2_rdnc_model)


    def calcgainmodel(self,t,g0=2):
        t.configcols(['gain_theo','gain0'],'f','%.3f',visible=1)
        keys = t.CUT_inrange('g',g0,g0)
        gain0= t.calcaverage(keys,'gain')
        dS0  = t.calcaverage(keys,'dS_mean')
        if self.verbose:
            print('gain0:',gain0)
            print('dS0:',dS0)
        for key in t.allrowkeys:
            if t.getentry(key,'dS_mean')>=10.0:
                t.setentry(key,'gain_theo',gain0*dS0/t.getentry(key,'dS_mean'))
            if t.getentry(key,'gain')!=None:
                t.setentry(key,'gain0',t.getentry(key,'gain')*t.getentry(key,'dS_mean')/dS0)


    def calceffgain_dS(self,t):
        t.configcols(['gain'],'f','%.3f',visible=1)
        for key in t.allrowkeys:
            #if t.getentry(key,'S_mean')<=35000.0:
            if t.getentry(key,'stdev2_rdnoisecorr')>=3.0:
                gain = 2.0*t.getentry(key,'dS_mean')/t.getentry(key,'stdev2_rdnoisecorr')
                t.setentry(key,'gain',gain)


    def subtract_stdev2_g0(self,t):
        keys0 = t.CUT_inrange('g',0,0)
        stdev2_g0 = t.calcaverage(keys0,'stdev2_rdnoisecorr')
        keys = t.CUT_none_vals('stdev2_rdnoisecorr')
        for key in keys:
            t.setentry(key,'stdev2_rdnoisecorr',t.getentry(key,'stdev2_rdnoisecorr')-stdev2_g0)


    def getfitsfile(self,filename):
        if self.verbose:
            print('# Loading {}'.format(filename))
        #self.filename=filename
        data, hdr = fits.getdata(filename, 0, header=True)

        if not ('NGROUPS' in hdr):
            hdr = fits.getheader(filename,0)
        if not ('NGROUPS' in hdr):
            print('ERROR, could not find NGROUPS in header!!')
            sys.exit(0)

        if len(data.shape)<3:
            print('ERROR: data cube needs to have at least 3 dimensions')
            print('for {}, only {}!!!'.format(filename,len(data.shape)))
            sys.exit(0)

        Ngroups = hdr['NGROUPS']
        if self.verbose>1:
             print('Ngroups ',Ngroups)

        Nint = hdr['NINTS']
        if self.verbose>1:
             print('Nint ',Nint)

        if len(data.shape)==3:
            if self.verbose>1:
                print('CV data format')
                # make sure data is in float!
                data=data*1.0
        elif len(data.shape)==4:
            if self.verbose>1:
                print('SSB data format')
            scinew=scipy.zeros((Nint*Ngroups,data.shape[2],data.shape[3]), dtype=float)
            for i in range(Nint):
                scinew[i:i+Ngroups,:,:]=data[i,:,:,:]
            data = scinew
        else:
            sys.exit(0)

        # Make sure the NAXIS1,2 are consistent for all images
        if self.Nx==None:
            self.Nx = data.shape[2]
        else:
            if data.shape[2] !=self.Nx:
                print('ERROR! Nx {}!={}'.format(data.shape[2],self.Nx))
                sys.exit(0)

        if self.Ny==None:
            self.Ny = data.shape[1]
        else:
            if data.shape[1]!=self.Ny:
                print('ERROR! Nx {}!={}'.format(data.shape[1],self.Ny))
                sys.exit(0)

        return(data,hdr)


    def darks4readnoise(self,darks4readnoiselist,boxmap,box_location,mask=None,gmax=None,
                        Pmax=20.0,Nmin=40,saveramps=False,
                        skiprefpixcorr=False,refpixboxsize=10,skipcorrect4framereset=False):

        print('\n###### Calculating the readnoise from the double difference of two dark frames')
        if self.readnoisemethod!='doublediff':
            raise RuntimeError('ERROR: attempting to calculate the readnoise from the dark frames, but readnoise methos = {}!'.format(self.readnoisemethod))

        if len(darks4readnoiselist)<2:
            print('Not enough darks, at least 2 needed!',darks4readnoiselist)
            sys.exit(0)

        data1,hdr1 = self.getfitsfile(darks4readnoiselist[0])
        data2,hdr2 = self.getfitsfile(darks4readnoiselist[1])

        if gmax == None:
            gmax = data1.shape[0]
        else:
            if data1.shape[0]<gmax:
                gmax = data1.shape[0]

        if self.verbose>1:
            print('gmax for dark frames:',gmax)

        #print 'VVVVVVVVV11'
        #fits.writeto('test.diffim.fits',data1[1:,:,:]-data1[0,:,:],clobber=True)
        # make the frame reset correcton and side ref pix correction if wanted
        self.make_refpixcorr(data=data1,gmax=gmax,skiprefpixcorr=skiprefpixcorr,
                             refpixboxsize=refpixboxsize,
                             skipcorrect4framereset=skipcorrect4framereset)
        #print 'VVVVVVVVV'
        #fits.writeto('test.finaldiffim.fits',data1[1:,:,:]-data1[0,:,:],clobber=True)
        #sys.exit(0)
        self.make_refpixcorr(data=data2,gmax=gmax,skiprefpixcorr=skiprefpixcorr,
                             refpixboxsize=refpixboxsize,
                             skipcorrect4framereset=skipcorrect4framereset)

        # double difference. stepsize is two groups so that the differences are completely independent!
        diffimdoublediff = (data1[1:gmax:2,:,:]-data1[0:gmax-1:2,:,:]) - (data2[1:gmax:2,:,:]-data2[0:gmax-1:2,:,:])
        #fits.writeto('diffimdoublediff.fits',data=diffimdoublediff)

        # define image array ramps for readnoise determined from the double difference of dark groups
        rdnoiseramp={}
        rdnoiseramp['im']= scipy.zeros(diffimdoublediff.shape,dtype=float)
        rdnoiseramp['err']= scipy.zeros(diffimdoublediff.shape,dtype=float)
        rdnoiseramp['flag']= scipy.ones(diffimdoublediff.shape,dtype=np.int16)
        rdnoiseramp['meanval']=scipy.zeros(diffimdoublediff.shape,dtype=float)
        rdnoiseramp['Pskip']= scipy.zeros(diffimdoublediff.shape,dtype=float)
        rdnoiseramp['Nused']=scipy.zeros(diffimdoublediff.shape,dtype=np.int32)

        # define image arrays for average readnoise determined from the different groups of the double difference of dark groups
        self.rdnoise['im']= scipy.zeros((self.Ny,self.Nx),dtype=float)
        self.rdnoise['err']= scipy.zeros((self.Ny,self.Nx),dtype=float)
        self.rdnoise['flag']= scipy.ones((self.Ny,self.Nx),dtype=np.int16)
        self.rdnoise['stdev']=scipy.zeros((self.Ny,self.Nx),dtype=float)
        self.rdnoise['Pskip']= scipy.zeros((self.Ny,self.Nx),dtype=float)
        self.rdnoise['Nused']=scipy.zeros((self.Ny,self.Nx),dtype=np.int32)

        # for consistency check: readnoise determined from the single difference of dark groups
        diffimsinglediff = (data1[1:gmax:2,:,:]-data1[0:gmax-1:2,:,:])
        rdnoiseramp['singlediff']= scipy.zeros(diffimsinglediff.shape,dtype=float)
        self.rdnoise['singlediff']= scipy.zeros((self.Ny,self.Nx),dtype=float)

        sigmacut = calcaverageclass()

        for boxval in box_location:
            x, y = box_location[boxval]
            goodpix = boxmap == boxval
            #x=xinfo[0]
            # loop through all boxes
            #while x<xinfo[1] and x+xinfo[2]<=self.Nx:
            #y=yinfo[0]
            #print('x={}'.format(x))
            #while y<yinfo[1] and y+yinfo[2]<=self.Ny:
            mask2use = None
            if not (mask is None):
                #mask2use = mask[y:y+yinfo[2],x:x+xinfo[2]]
                mask2use = mask[goodpix]

            # loop through all boxes and groups
            for g in range(0,diffimdoublediff.shape[0]):
                # get statistics for double diff
                #sigmacut.calcaverage_sigmacutloop(diffimdoublediff[g,y:y+yinfo[2],x:x+xinfo[2]],
                #                                  mask=mask2use,Nsigma=3.0,verbose=3)
                temp = diffimdoublediff[g,:,:]
                sigmacut.calcaverage_sigmacutloop(temp[goodpix],
                                                  mask=mask2use,Nsigma=3.0,verbose=3)

                #print sigmacut.__str__()
                if sigmacut.mean!=None:
                    readim = rdnoiseramp['im'][g,:,:]
                    readim[goodpix] = sigmacut.stdev*0.5
                    rdnoiseramp['im'][g,:,:] = readim

                    readim = rdnoiseramp['Nused'][g,:,:]
                    readim[goodpix] = sigmacut.Nused
                    rdnoiseramp['Nused'][g,:,:] = readim

                    readim = rdnoiseramp['err'][g,:,:]
                    readim[goodpix] = sigmacut.stdev_err*0.5
                    rdnoiseramp['err'][g,:,:] = readim

                    Pskip=100.0*sigmacut.Nskipped/(sigmacut.Nused+sigmacut.Nskipped)
                    readim = rdnoiseramp['Pskip'][g,:,:]
                    readim[goodpix] = Pskip
                    rdnoiseramp['Pskip'][g,:,:] = readim

                    readim = rdnoiseramp['meanval'][g,:,:]
                    readim[goodpix] = sigmacut.mean
                    rdnoiseramp['meanval'][g,:,:] = readim

                    #rdnoiseramp['im'][g,y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.stdev*0.5
                    #rdnoiseramp['Nused'][g,y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.Nused
                    #rdnoiseramp['err'][g,y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.stdev_err*0.5
                    #rdnoiseramp['Pskip'][g,y:y+yinfo[2],x:x+xinfo[2]]=Pskip
                    #rdnoiseramp['meanval'][g,y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.mean

                    readim = rdnoiseramp['flag'][g,:,:]
                    if Pskip<=Pmax and sigmacut.Nused>=Nmin:
                        #rdnoiseramp['flag'][g,y:y+yinfo[2],x:x+xinfo[2]]=0
                        readim[goodpix] = 0
                    else:
                        #rdnoiseramp['flag'][g,y:y+yinfo[2],x:x+xinfo[2]]=1
                        readim[goodpix] = 1
                    rdnoiseramp['flag'][g,:,:] = readim

                # get statistics for single diff
                temp = diffimsinglediff[g,:,:]
                sigmacut.calcaverage_sigmacutloop(temp[goodpix],
                                                  mask=mask2use,Nsigma=3.0,verbose=0)
                temprd = rdnoiseramp['singlediff'][g,:,:]
                temprd[goodpix] = sigmacut.stdev*0.70710678118655
                rdnoiseramp['singlediff'][g,:,:] = temprd
                #rdnoiseramp['singlediff'][g,y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.stdev*0.70710678118655

            #(ymid,xmid) = (int(y+0.5*yinfo[2]),int(x+0.5*xinfo[2]))
            sigmacut.calcaverage_sigmacutloop(rdnoiseramp['im'][:,y,x],
                                              mask=rdnoiseramp['flag'][:,y,x],
                                              Nsigma=3.0,verbose=0)
            if sigmacut.mean is not None:
                #self.rdnoise['im'][y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.mean
                #self.rdnoise['Nused'][y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.Nused
                #self.rdnoise['err'][y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.mean_err
                #Pskip=100.0*sigmacut.Nskipped/(sigmacut.Nused+sigmacut.Nskipped)
                #self.rdnoise['Pskip'][y:y+yinfo[2],x:x+xinfo[2]]=Pskip
                #self.rdnoise['stdev'][y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.stdev
                self.rdnoise['im'][goodpix]=sigmacut.mean
                self.rdnoise['Nused'][goodpix]=sigmacut.Nused
                self.rdnoise['err'][goodpix]=sigmacut.mean_err
                Pskip=100.0*sigmacut.Nskipped/(sigmacut.Nused+sigmacut.Nskipped)
                self.rdnoise['Pskip'][goodpix]=Pskip
                self.rdnoise['stdev'][goodpix]=sigmacut.stdev
                if sigmacut.Nused>=Nmin:
                    #self.rdnoise['flag'][y:y+yinfo[2],x:x+xinfo[2]]=0
                    self.rdnoise['flag'][goodpix]=0
                else:
                    #self.rdnoise['flag'][y:y+yinfo[2],x:x+xinfo[2]]=1
                    self.rdnoise['flag'][goodpix]=1
            sigmacut.calcaverage_sigmacutloop(rdnoiseramp['singlediff'][:,y,x],
                                              mask=rdnoiseramp['flag'][:,y,x],
                                              Nsigma=3.0,verbose=0)
            #self.rdnoise['singlediff'][y:y+yinfo[2],x:x+xinfo[2]]=sigmacut.mean
            self.rdnoise['singlediff'][goodpix]=sigmacut.mean

            #y+=yinfo[2]
            #x+=xinfo[2]

        # save the results in an hdulist
        rdnoisefilename = self.outfilebasename+'.doublediffrdnoise.fits'
        readnoisehdulist = fits.HDUList([fits.PrimaryHDU(data=None),
                                         fits.ImageHDU(self.rdnoise['im'],name='rdnoise'),
                                         fits.ImageHDU(self.rdnoise['err'],name='rdnoiseerror'),
                                         fits.ImageHDU(self.rdnoise['stdev'],name='stdev'),
                                         fits.ImageHDU(self.rdnoise['flag'],name='flag'),
                                         fits.ImageHDU(self.rdnoise['Nused'],name='Nused'),
                                         fits.ImageHDU(self.rdnoise['Pskip'],name='Pskip'),
                                         fits.ImageHDU(self.rdnoise['singlediff'],name='rdnoiseimsinglediff')])
        print('Saving rdnoise',rdnoisefilename)
        rmfile(rdnoisefilename)
        readnoisehdulist.writeto(rdnoisefilename)

        if saveramps:
            # save the results in an hdulist
            rdnoiserampfilename = self.outfilebasename+'.doublediffrdnoiseramp.fits'
            readnoiseramphdulist = fits.HDUList([fits.PrimaryHDU(data=None),
                                                 fits.ImageHDU(rdnoiseramp['im'],name='rdnoise'),
                                                 fits.ImageHDU(rdnoiseramp['err'],name='rdnoiseerror'),
                                                 fits.ImageHDU(rdnoiseramp['flag'],name='flag'),
                                                 fits.ImageHDU(rdnoiseramp['Nused'],name='Nused'),
                                                 fits.ImageHDU(rdnoiseramp['Pskip'],name='Pskip'),
                                                 fits.ImageHDU(rdnoiseramp['meanval'],name='meanval'),
                                                 fits.ImageHDU(rdnoiseramp['singlediff'],name='rdnoiseimsinglediff')])
            print('Saving rdnoise ramp',rdnoiserampfilename)
            rmfile(rdnoiserampfilename)
            readnoiseramphdulist.writeto(rdnoiserampfilename)
            del readnoiseramphdulist

        del data1,data2,hdr1,hdr2,rdnoiseramp,diffimdoublediff,diffimsinglediff,readnoisehdulist


    def mkbpmmask(self,flatdata,hdr,boxmap,box_location,fluxcutfraction=0.7,
                  minmedian=35000.0,maxmedian=40000.0):
        #print('\n##### Making bpm mask',xinfo,yinfo)
        bpm = scipy.zeros((self.Ny,self.Nx),dtype=np.int16)
        flatdata_corr = flatdata[1:,:,:]-flatdata[0,:,:]

        g0 = None
        for g in range(flatdata_corr.shape[0]-1,0,-1):
            if self.verbose>1:
                print('g={}'.format(g))
            median = np.median(flatdata_corr[g,:,:])
            if (median>minmedian) and (median<maxmedian):
                g0=g
            if (median<=minmedian):
                if g0==None:
                    g0=g
                break
        if g0==None:
            raise RuntimeError('Could not determine a good g for which the median flux is between %.0f and %.0f' % (minmedian,maxmedian))

        flatdata_corr0=flatdata_corr[g0,:,:]

        flatdata4bpmfilename= self.outfilebasename+'.flatdata4bpm.fits'
        print('Saving flatdata4bpm:',flatdata4bpmfilename)
        rmfile(flatdata4bpmfilename)
        fits.writeto(flatdata4bpmfilename,data=flatdata_corr0)

        # crude cut: get rid of overscan and the most blatant bad pixels
        median0 = np.median(flatdata_corr0)
        bpm[scipy.where(flatdata_corr0<fluxcutfraction*median0)]=1
        bpm[scipy.where(np.isnan(flatdata_corr0))]=1
        bpm[:4,:]=1
        bpm[-4:,:]=1
        bpm[:,:4]=1
        bpm[:,-4:]=1

        #bpmfilename= self.outfilebasename+'.bpmmaskcrude.fits'
        #print 'Saving bpm:',bpmfilename
        #rmfile(bpmfilename)
        #fits.writeto(bpmfilename,data=bpm)

        sigmacut = calcaverageclass()
        for boxval in box_location:
            if np.isfinite(boxval):
                x, y = box_location[boxval]
                goodpix = boxmap == boxval
                #x=xinfo[0]
                # loop through all boxes
                #while x<xinfo[1] and x+xinfo[2]<=self.Nx:
                #y=yinfo[0]
                #while y<yinfo[1] and y+yinfo[2]<=self.Ny:

                #flatbox = flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]]
                flatbox = flatdata_corr0[goodpix]

                #medianflux = np.median(flatdata_corr[g0,y:y+yinfo[2],x:x+xinfo[2]])
                #medianflux = np.median(flatbox[scipy.where(bpm[y:y+yinfo[2],x:x+xinfo[2]]==0)])
                #if self.options.verbose>1:
                #    print 'x=%d, y=%d, median flux=%f' % (x,y,medianflux)
                #print 'TEST',np.median(flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]])

                #flatvals = sorted(flatbox[scipy.where(bpm[y:y+yinfo[2],x:x+xinfo[2]]==0)])
                flatvals = sorted(flatbox[scipy.where(bpm[goodpix] == 0)])
                N = len(flatvals)
                if N == 0:
                    print("No good points in flat: box {}".format(boxval))
                    sys.exit()
                medianflux = flatvals[int(N*0.5)]
                fluxplus = flatvals[int(N*(0.5+0.5*0.6827))]
                fluxminus = flatvals[int(N*(0.5-0.5*0.6827))]
                if self.verbose>1:
                    print('median flux: %.1f, sigma+=%.1f, sigma-=%.1f' % (medianflux,np.abs(medianflux-fluxplus),np.abs(medianflux-fluxminus)))
                sigma=min([np.abs(medianflux-fluxplus),np.abs(medianflux-fluxminus)])
                Nsigma=3.0
                if self.aggressivemasking:
                    #bpm[y:y+yinfo[2],x:x+xinfo[2]] |= scipy.where(flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]]<medianflux-Nsigma*sigma,2,0)
                    #bpm[y:y+yinfo[2],x:x+xinfo[2]] |= scipy.where(flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]]>medianflux+Nsigma*sigma,4,0)
                    bpm[goodpix] |= scipy.where(flatdata_corr0[goodpix]<medianflux-Nsigma*sigma,2,0)
                    bpm[goodpix] |= scipy.where(flatdata_corr0[goodpix]>medianflux+Nsigma*sigma,4,0)
                else:
                    #bpm[y:y+yinfo[2],x:x+xinfo[2]] |= scipy.where(flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]]<medianflux-Nsigma*np.abs(medianflux-fluxminus),2,0)
                    #bpm[y:y+yinfo[2],x:x+xinfo[2]] |= scipy.where(flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]]>medianflux+Nsigma*np.abs(medianflux-fluxplus),4,0)
                    bpm[goodpix] |= scipy.where(flatdata_corr0[goodpix]<medianflux-Nsigma*np.abs(medianflux-fluxminus),2,0)
                    bpm[goodpix] |= scipy.where(flatdata_corr0[goodpix]>medianflux+Nsigma*np.abs(medianflux-fluxplus),4,0)

                if self.verbose>1:
                    print('{}% of pixels masked in this box'.format(100.0*np.count_nonzero(bpm[goodpix])/np.sum(goodpix)))

                # old way
                # get statistics for double diff
                #sigmacut.calcaverage_sigmacutloop(flatdata_corr0[y:y+yinfo[2],x:x+xinfo[2]],mask=bpm[y:y+yinfo[2],x:x+xinfo[2]],fixmean=medianflux,Nsigma=3.0,saveused=True,verbose=0)
                #print sigmacut.__str__()
                #if sigmacut.mean!=None:
                #    bpm[y:y+yinfo[2],x:x+xinfo[2]] |= sigmacut.clipped*2

                #y+=yinfo[2]
                #x+=xinfo[2]

        bpmfilename= self.outfilebasename+'.bpmmask.fits'
        print('Saving bpm:',bpmfilename)
        rmfile(bpmfilename)
        fits.writeto(bpmfilename,data=bpm)
        return(bpm)


    def make_refpixcorr(self,data,gmax=None,skiprefpixcorr=False,refpixboxsize=10,
                        skipcorrect4framereset=False):
        test=False
        if gmax == None:
            gmax = data.shape[0]

        framereset = frameresetclass()
        refpixcorr = refpixcorrclass()

        if not skipcorrect4framereset:
            framereset.correct4framereset(data,method='first')

        if not skiprefpixcorr:
            refpixcorr.outfilebasename='test'
            refpixcorr.verbose=0
            refpixcorr.mkframeresetcorrection(data,boxsize=refpixboxsize,testflag=False)

        if not skipcorrect4framereset:
            framereset.correct4framereset(data,method='next')


    def make_refpixcorr_delme(self,data,gmax=None,skiprefpixcorr=False,refpixboxsize=10,
                              skipcorrect4framereset=False):
        test=False
        if gmax == None:
            gmax = data.shape[0]

        if not skiprefpixcorr:
            # refpix correction
            refpixcorr = refpixcorrclass()
            refpixcorr.outfilebasename='test'
            refpixcorr.verbose=3
            refpixcorr.data = data
            refpixcorr.mkdiffim(savediffimflag=True)
            refpixcorr.calcbiasdriftvec_allgroups(boxsize=refpixboxsize)
            refpixcorr.savebiasdriftvecs_as_txt(boxsize=refpixboxsize)
            refpixcorr.savebiasdriftvecs(boxsize=refpixboxsize)
            refpixcorr.correctimage4leftrightrefpix(boxsize=refpixboxsize)
            refpixcorr.mkdiffimcorr(savediffimflag=True)

        if not skipcorrect4framereset:
            sigmacutrefpixbottom = calcaverageclass()
            sigmacutrefpixtop = calcaverageclass()
            print('Frame reset correction!')
            for g in range(self.gmin+1,gmax):
                if self.verbose:
                    print('g',g)

                for amp in range(1,5):
                    if self.verbose:
                        print('amp',amp)
                    xmin4amp = (amp-1)*512
                    xmax4amp = amp*512

                    # image 1
                    refpixbottom = data[g,:4,xmin4amp:xmax4amp]-data[g-1,:4,xmin4amp:xmax4amp]
                    sigmacutrefpixbottom.calcaverage_sigmacutloop(refpixbottom,Nsigma=3,verbose=3)
                    if sigmacutrefpixbottom.mean==None:
                        raise RuntimeError("Could not correct for frame reset!")

                    refpixtop = data[g,-4:,xmin4amp:xmax4amp]-data[g-1,-4:,xmin4amp:xmax4amp]
                    if test: refpixtop = data[g,:4,xmin4amp:xmax4amp]-data[g-1,:4,xmin4amp:xmax4amp]
                    sigmacutrefpixtop.calcaverage_sigmacutloop(refpixtop,Nsigma=3,verbose=3)
                    if sigmacutrefpixtop.mean==None:
                        raise RuntimeError("Could not correct for frame reset!")

                    frameresetcorrection = 0.5*(sigmacutrefpixbottom.mean+sigmacutrefpixtop.mean)
                    data[g:,:,xmin4amp:xmax4amp] -= frameresetcorrection

    def initcols(self,ts):
        for t in ts:
            t.configcols(['g'],'d','%d',visible=1)
            t.configcols(['geff'],'f','%.1f',visible=1)
            t.configcols(['xmin','xdelta','ymin','ydelta'],'d','%d',visible=1)
            t.configcols(['stdev2','stdev2_rdnoisecorr'],'f','%.3f',visible=1)
        return ts


    #def gainim(self,infile1,infile2,darks4readnoiselist=None,pattern='.fits',format='g%03d',imtype=None,
    #           skiprefpixcorr=False,
    #           refpixboxsize=12,
    #           skipcorrect4framereset=False):

    def gainim(self,infile1,infile2,pattern='.fits',format='g%03d',imtype=None,refpixboxsize=12):
        # Lazy updates to allow for more use of class variables rather than
        # needing so many in the function definition
        skiprefpixcorr = self.skiprefpixcorr
        skipcorrect4framereset = self.skipcorrect4framereset
        darks4readnoiselist = self.darks
        self.flatfile1 = infile1
        self.flatfile2 = infile2

        # Get the proper darks
        if self.darks[0] is not None:
            darkfilelist = [os.path.abspath(self.darks[0]),os.path.abspath(self.darks[1])]
            print('Using the following dark frames:',darkfilelist)
        elif self.autodark4readnoise:
            m = re.search('(\d+)_SE',os.path.basename(infile1))
            if m is None:
                raise RuntimeError('file {}'.format(filename))
            SCA = int(m.groups()[0])
            path = os.path.dirname(os.path.abspath(infile1))
            darkfilelist = glob.glob(path+'/*DARK*_%d_SE*.fits'.format(SCA))
            if len(darkfilelist) == 0:
                raise RuntimeError('Could not find dark frames *DARK*.fits in %s'.format(path))
            print('Found the following dark frames: {}'.format(darkfilelist))
        else:
            darkfilelist = None

        if darks4readnoiselist[0] is not None:
            self.readnoisemethod='doublediff'
        else:
            self.readnoisemethod='dsub12'

        self.outfilebasename = re.sub('\.fits$','',infile1)+'.%s' % self.readnoisemethod
        if self.outsubdir!=None:
            (head,tail)=os.path.split(os.path.abspath(self.outfilebasename))
            #newoutdir = '%s/%s' % (head,self.outsubdir)
            newoutdir = os.path.join(head,self.outsubdir)
            if not os.path.isdir(newoutdir):
                os.makedirs(newoutdir)
            if not os.path.isdir(newoutdir):
                raise RuntimeError('ERROR: Cannot create directory {}'.format(newoutdir))

            #self.outfilebasename = '%s/%s' % (newoutdir,tail)
            self.outfilebasename = os.path.join(newoutdir,tail)

        if self.verbose:
            print('outfilebasename: {}'.format(self.outfilebasename))

        if self.test:
            boxsize = 32
            xmin = 1024
            ymin = 0
            xmax = int(xmin+2*boxsize)
            ymax = int(ymin+2*boxsize)
            if self.gmax==None:
                self.gmax = 5
            self.gmax4dark=10
        else:
            xmin = self.xmin
            xmax = self.xmax
            ymin = self.ymin
            ymax = self.ymax
            boxsize = self.boxsize

        print('x:%d-%d  y:%d-%d' % (xmin,xmax,ymin,ymax))

        # load the flats
        print("Reading in flats")
        self.data1,self.hdr1 = self.getfitsfile(infile1)
        self.data2,self.hdr2 = self.getfitsfile(infile2)
        detector = self.hdr1['DETECTOR']

        # Generate a grid of boxes to use for the calculations
        # This includes the effects of the low epoxy region
        # if applicable
        print("Generating box map")
        self.boxmap, self.box_centers = self.make_box_map(boxsize,detector)

        #hhh = fits.PrimaryHDU()
        #hh1 = fits.ImageHDU(self.boxmap)
        #hhl = fits.HDUList([hhh,hh1])
        #hhl.writeto('test_boxmap.fits',overwrite=True)


        # make the mask files
        if self.verbose:
            print("Making bad pixel map")
        self.bpmmask = self.mkbpmmask(self.data1,self.hdr1,self.boxmap,self.box_centers)

        # get the readnoise from the dark frames if wanted...
        if darks4readnoiselist[0] is not None:
            #self.darks4readnoise(darks4readnoiselist,(xmin,xmax,boxsize),
            #                     (ymin,ymax,boxsize),mask=self.bpmmask,
            #                     gmax=self.gmax4dark,skiprefpixcorr=skiprefpixcorr,
            #                     refpixboxsize=refpixboxsize )
            print("Calculating readnoise from darks")
            self.darks4readnoise(darks4readnoiselist,self.boxmap,
                                 self.box_centers,mask=self.bpmmask,
                                 gmax=self.gmax4dark,skiprefpixcorr=skiprefpixcorr,
                                 refpixboxsize=refpixboxsize )

        Ngroups = self.hdr2[self.Ngroupskeyword]
        if self.verbose:
             print('Ngroups ',Ngroups)
        tsubik = txttableclass()
        #tsubik_subgm1 = txttableclass()
        tsubg0 = txttableclass()
        tdsub1 = txttableclass()
        tdsub2 = txttableclass()
        tdsubij = txttableclass()
        tsubik,tsubg0,tdsub1,tdsub2,tdsubij = self.initcols([tsubik,tsubg0,tdsub1,tdsub2,tdsubij])

        if self.gmax == None:
            gmax = Ngroups
        else:
            gmax = min(self.gmax,Ngroups)

        if not skiprefpixcorr:
            print("Performing refpix correction")
        self.make_refpixcorr(data=self.data1,gmax=gmax,skiprefpixcorr=skiprefpixcorr,
                             refpixboxsize=refpixboxsize,
                             skipcorrect4framereset=skipcorrect4framereset)
        self.make_refpixcorr(data=self.data2,gmax=gmax,skiprefpixcorr=skiprefpixcorr,
                             refpixboxsize=refpixboxsize,
                             skipcorrect4framereset=skipcorrect4framereset)
        #print 'VVVVVVVVV'
        #fits.writeto('test.finaldiffim.fits',self.data1[1:,:,:]-self.data1[0,:,:],clobber=True)
        #sys.exit(0)

        #if not skipcorrect4framereset:
        #    sigmacutrefpixbottom = calcaverageclass()
        #    sigmacutrefpixtop = calcaverageclass()
        #    print 'Frame reset correction!'
        #    for g in range(self.options.gmin,gmax):
        #        if self.options.verbose:
        #            print 'g',g

        #        for amp in range(1,5):
        #            xmin4amp = (amp-1)*512
        #            xmax4amp = amp*512
        #
        #            # image 1
        #            refpixbottom = self.data1[g,:4,xmin4amp:xmax4amp]-self.data1[g-1,:4,xmin4amp:xmax4amp]
        #            sigmacutrefpixbottom.calcaverage_sigmacutloop(refpixbottom,Nsigma=3,verbose=3)
        #            if sigmacutrefpixbottom.mean==None:
        #                raise RuntimeError,"Could not correct for frame reset!"

        #            #refpixtop = self.data1[g,:4,xmin4amp:xmax4amp]-self.data1[g-1,:4,xmin4amp:xmax4amp]
        #            refpixtop = self.data1[g,-4:,xmin4amp:xmax4amp]-self.data1[g-1,-4:,xmin4amp:xmax4amp]
        #            sigmacutrefpixtop.calcaverage_sigmacutloop(refpixtop,Nsigma=3,verbose=3)
        #            if sigmacutrefpixtop.mean==None:
        #                raise RuntimeError,"Could not correct for frame reset!"

        #            frameresetcorrection = 0.5*(sigmacutrefpixbottom.mean+sigmacutrefpixtop.mean)
        #            self.data1[g:,:,xmin4amp:xmax4amp] -= frameresetcorrection

        #            # image 2
        #            refpixbottom = self.data2[g,:4,xmin4amp:xmax4amp]-self.data2[g-1,:4,xmin4amp:xmax4amp]
        #            sigmacutrefpixbottom.calcaverage_sigmacutloop(refpixbottom,Nsigma=3,verbose=3)
        #            if sigmacutrefpixbottom.mean==None:
        #                raise RuntimeError,"Could not correct for frame reset!"

        #            #refpixtop = self.data2[g,:4,xmin4amp:xmax4amp]-self.data2[g-1,:4,xmin4amp:xmax4amp]
        #            refpixtop = self.data2[g,-4:,xmin4amp:xmax4amp]-self.data2[g-1,-4:,xmin4amp:xmax4amp]
        #            sigmacutrefpixtop.calcaverage_sigmacutloop(refpixtop,Nsigma=3,verbose=3)
        #            if sigmacutrefpixtop.mean==None:
        #                raise RuntimeError,"Could not correct for frame reset!"

        #            frameresetcorrection = 0.5*(sigmacutrefpixbottom.mean+sigmacutrefpixtop.mean)
        #            self.data2[g:,:,xmin4amp:xmax4amp] -= frameresetcorrection


        print('\n##### Looping through boxes, calculating stdevs....')
        for g in range(self.gmin,gmax):
            if self.verbose:
                print('g',g)

            #print('calculating diffims')
            #diffim =  self.data2[g,:,:]-self.data1[g,:,:]
            diffim =  self.data1[g,:,:]-self.data2[g,:,:]

            if g>0:
                diffim_subg0 =  (self.data1[g,:,:]-self.data1[0,:,:])-(self.data2[g,:,:]-self.data2[0,:,:])
                # calculate the average count between g and g-1
                diffim_dsubij =  (self.data1[g,:,:]-self.data1[g-1,:,:])-(self.data2[g,:,:]-self.data2[g-1,:,:])
                dS_dsubij = 0.5*(self.data1[g,:,:]-self.data1[g-1,:,:]+self.data2[g,:,:]-self.data2[g-1,:,:])
                S_dsubij = 0.25*(self.data1[g,:,:]+self.data1[g-1,:,:]+self.data2[g,:,:]+self.data2[g-1,:,:]-2.0*self.data1[0,:,:]-2.0*self.data2[0,:,:])

            if g>0 and g<Ngroups-1:
                diffim_dsub1 =  (self.data1[g+1,:,:]-self.data1[g,:,:])-(self.data1[g,:,:]-self.data1[g-1,:,:])
                dS_dsub1 = 0.5*(self.data1[g+1,:,:]-self.data1[g-1,:,:])

            if g>1 and g<Ngroups-1:
                diffim_dsub2 =  (self.data1[g+1,:,:]-self.data1[g,:,:])-(self.data1[g-1,:,:]-self.data1[g-2,:,:])
                dS_dsub2 = 1.0/3.0*(self.data1[g+1,:,:]-self.data1[g-2,:,:])
                S_dsub2 = 0.5*(self.data1[g,:,:]-self.data1[0,:,:]+self.data1[g-1,:,:]-self.data1[0,:,:])

            #im = 0.5*(self.data1[g,:,:]-self.data1[0,:,:]+self.data2[g,:,:]-self.data2[0,:,:])
            #im_subgm1 = 0.5*(self.data1[g,:,:]-self.data1[0,:,:]+self.data1[g-1,:,:]-self.data1[0,:,:])
            S = self.data1[g,:,:]-self.data1[0,:,:]
            #dS = self.data1[g,:,:]-self.data1[g-1,:,:]

            #self.calcgain_grid(diffim,(xmin,xmax,boxsize),(ymin,ymax,boxsize),
            #                   tsubik,g,mask=self.bpmmask,S=S,imtype="default")
            self.calcgain_grid(diffim,tsubik,g,self.boxmap,self.box_centers,mask=self.bpmmask,S=S,imtype="default")

            if g>0:
                #self.calcgain_grid(diffim_subg0,(xmin,xmax,boxsize),(ymin,ymax,boxsize),
                #                   tsubg0,g,mask=self.bpmmask,S=S,imtype="subg0")
                #self.calcgain_grid(diffim_dsubij,(xmin,xmax,boxsize),(ymin,ymax,boxsize),
                #                   tdsubij,g,mask=self.bpmmask,S=S_dsubij,dS=dS_dsubij,imtype="dsubij")
                self.calcgain_grid(diffim_subg0,tsubg0,g,self.boxmap,self.box_centers,mask=self.bpmmask,S=S,imtype="subg0")
                self.calcgain_grid(diffim_dsubij,tdsubij,g,self.boxmap,self.box_centers,mask=self.bpmmask,S=S_dsubij,dS=dS_dsubij,imtype="dsubij")

            if g>0 and g<Ngroups-1:
                #self.calcgain_grid(diffim_dsub1,(xmin,xmax,boxsize),(ymin,ymax,boxsize),
                #                   tdsub1,g,mask=self.bpmmask,S=S,dS=dS_dsub1,imtype="dsub1")
                self.calcgain_grid(diffim_dsub1,tdsub1,g,self.boxmap,self.box_centers,mask=self.bpmmask,S=S,dS=dS_dsub1,imtype="dsub1")
            if g>1 and g<Ngroups-1:
                #self.calcgain_grid(diffim_dsub2,(xmin,xmax,boxsize),(ymin,ymax,boxsize),
                #                   tdsub2,g,mask=self.bpmmask,S=S_dsub2,dS=dS_dsub2,imtype="dsub2")
                self.calcgain_grid(diffim_dsub2,tdsub2,g,self.boxmap,self.box_centers,mask=self.bpmmask,S=S_dsub2,dS=dS_dsub2,imtype="dsub2")

            if 0==1:
                outname = self.outfilebasename+'.%02d.diff.fits' % g
                print('Saving ',outname)
                rmfile(outname)
                fits.writeto(outname,diffim,self.hdr1)

                outname = self.outfilebasename+'.%02d.diff.subg0.fits' % g
                print('Saving ',outname)
                rmfile(outname)
                fits.writeto(outname,diffim_subg0,self.hdr1)

                #outname = self.outfilebasename+'.%02d.diff.subgm1.fits' % g
                #print 'Saving ',outname
                #rmfile(outname)
                #fits.writeto(outname,diffim_subgm1,self.hdr1)

            print('cleaning up')
            del(diffim,S)
            if g>0:
                del(diffim_subg0,diffim_dsubij,dS_dsubij,S_dsubij)
            if g>0 and g<Ngroups-1:
                del(diffim_dsub1,dS_dsub1)
            if g>1 and g<Ngroups-1:
                del(diffim_dsub2,dS_dsub2,S_dsub2)

        #tdsub1.save2file('tdsub1.txt')
        #tdsub2.save2file('tdsub2.txt')
        #tdsubij.save2file('tdsubij.txt')
        #print('tdsub printed for testing')

        self.correct_stdev2_with_rdnoiseterms(self.boxmap,self.box_centers,
                                              tsubik,tsubg0,tdsub1,tdsub2,tdsubij,
                                              Smean_max=self.Smeanmax4fit)

        # calculate the effective gain for each group
        self.calceffgain_dS(tdsub1)
        self.calceffgain_dS(tdsub2)
        self.calceffgain_dS(tdsubij)

        #self.calcgain((xmin,xmax,boxsize),(ymin,ymax,boxsize),tdsub1,'dsub1',
        #              Smean_max=self.options.Smeanmax4fit)
        #self.calcgain((xmin,xmax,boxsize),(ymin,ymax,boxsize),tdsub2,'dsub2',
        #              Smean_max=self.options.Smeanmax4fit)
        #self.calcgain((xmin,xmax,boxsize),(ymin,ymax,boxsize),tdsubij,'dsubij',
        #              Smean_max=self.options.Smeanmax4fit)
        self.calcgain(self.boxmap,self.box_centers,tdsub1,'dsub1',
                      Smean_max=self.Smeanmax4fit)
        self.calcgain(self.boxmap,self.box_centers,tdsub2,'dsub2',
                      Smean_max=self.Smeanmax4fit)
        self.calcgain(self.boxmap,self.box_centers,tdsubij,'dsubij',
                      Smean_max=self.Smeanmax4fit)

        self.calcstdev2_model(tdsub1,self.gmin,gmax,g0=1)
        self.calcstdev2_model(tdsub2,self.gmin,gmax,g0=2)
        self.calcstdev2_model(tdsubij,self.gmin,gmax,g0=1)

        self.calcgainmodel(tdsub1)
        self.calcgainmodel(tdsub2,g0=2)
        self.calcgainmodel(tdsubij)

        #self.calc_cum_stdev2(tdsub1,(xmin,xmax,boxsize),(ymin,ymax,boxsize),self.options.gmin,gmax)
        #self.calc_cum_stdev2(tdsub2,(xmin,xmax,boxsize),(ymin,ymax,boxsize),self.options.gmin,gmax)
        #self.calc_cum_stdev2(tdsubij,(xmin,xmax,boxsize),(ymin,ymax,boxsize),self.options.gmin,gmax)
        self.calc_cum_stdev2(tdsub1,self.box_centers,self.gmin,gmax)
        self.calc_cum_stdev2(tdsub2,self.box_centers,self.gmin,gmax)
        self.calc_cum_stdev2(tdsubij,self.box_centers,self.gmin,gmax)

        tsubik.save2file('%s.all.subij' % (self.outfilebasename),verbose=True)
        tsubg0.save2file('%s.all.subg0' % (self.outfilebasename),verbose=True)
        tdsub1.save2file('%s.all.dsub1' % (self.outfilebasename),verbose=True)
        tdsub2.save2file('%s.all.dsub2' % (self.outfilebasename),verbose=True)
        tdsubij.save2file('%s.all.dsubij' % (self.outfilebasename),verbose=True)


if __name__=='__main__':

    usagestring='USAGE: gain.py flatfilename1 flatfilename2'

    g = gainimclass()
    parser = g.add_options(usage = usagestring)
    args = parser.parse_args(namespace = g)
    #gainim.options,  args = parser.parse_args()

    #if len(args)!=2:
    #    sys.argv.append('--help')
    #    options,  args = parser.parse_args()
    #    sys.exit(0)

    #(infile1,infile2) = args

    #if re.search('DARK',infile1) or re.search('DARK',infile2):
    #    print 'DARK frame in infiles, exiting!!'
    #    print "SUCCESS calcgainim" # Say success here, makes the pipeline not barf...
    #    sys.exit(0)

    g.gainim(g.flatfile1,g.flatfile2,g.darks,
                  skiprefpixcorr=g.skiprefpixcorr,
                  refpixboxsize=g.refpixboxsize,
                  skipcorrect4framereset=g.skipcorrect4framereset)

    print("SUCCESS calcgainim")

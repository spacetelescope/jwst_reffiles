#! /usr/bin/env python

"""
This is an example module for creating individual readnoise refrence files.
This will serve as an example case for how to plug-in a user-created
module into jwst_reffiles.

In this example, we say that a readnoise imge is calcualted by creating
a series of CDS images from an input dark current ramp, and calculating
the sigma-clipped standard deviation in each of a grid of NxN pixel
boxes across the detector.

We'll also add some unique options in order to explore how to propagate
those into the jwst_reffiles world.
"""
import os
import copy

from astropy.io import fits
import numpy as np
from scipy.stats import sigmaclip

from jwst.datamodels import ReadnoiseModel


class MyReadnoise():
    def __init__(self, dark, boxsize=64, sigma_threshold=3, output_name='./'):
        """
        Parameters
        ----------
        dark : str
            Filename of dark current ramp. Assume dq_init has been run

        verbose : bool
            Print extra info while running


        NOTE: if you have 'verbose' as a parameter, it appears to collide
        with the verbose parameter in mkrefs. Armin uses integers for mkref's
        verbose, and I was using a True/False string, and mkrefs crashed bec
        it was trying to comapre a string to an integer. We need to be careful
        of collisions like this...

        """
        # Set up class variables
        self.dark = dark
        self.boxsize = int(boxsize)
        #self.verbose = True
        self.sigma_threshold = sigma_threshold
        self.output_name = output_name

        # Run
        self.compute()

    def box_grid(self):
        """Create a list of x,y coordinates to create a grid
        of boxes across the face of the detector
        """
        if 2048 % self.boxsize != 0:
            raise ValueError(("Warning: boxsize must be a factor of the detector width"))

        num_boxes = int(2048 / self.boxsize)
        xs = np.arange(0, 2049, self.boxsize)
        ys = np.arange(0, 2049, self.boxsize)

        # Exclude reference pixels
        xs[xs == 0] = 5
        ys[ys == 0] = 5
        xs[xs == 2048] = 2043
        ys[ys == 2048] = 2043
        return xs, ys

    def calculate_stdev(self, array, x, y, clipthresh):
        """Calculate the sigma-clipped stdev through the stack of array,
        within each box defined by the x and y coordinate values

        Parameters
        ----------
        array : numpy.ndarray
            3D array containing a stack of CDS images

        x : numpy.ndarray
            1D array containing x coordinate values corresponding to the
            edges of the grid of boxes

        y : numpy.ndarray
            1D array containing y coordinate values corresponding to the
            edges of the grid of boxes

        clipthresh : int
            Sigma-clipping threshold to use. Number of sigma. Default=3

        Returns
        -------
        readnoise : numpy.ndarray
            2D array containing a map of the readnoise across the detector
        """
        n_frames, n_y, n_x = array.shape

        # Create empty readnoise frame
        readnoise = np.zeros((n_y, n_x))

        # Loop over boxes and calculate stdev
        for i in range(len(x)-1):
            for j in range(len(y)-1):
                xstart = x[i]
                ystart = y[j]
                xend = x[i+1]
                yend = y[j+1]
                cube = array[:, ystart:yend, xstart:xend]
                clipped, thresh_low, thresh_high = sigmaclip(cube, low=clipthresh, high=clipthresh)
                dev = np.std(clipped)
                readnoise[ystart:yend, xstart:xend] = dev
        return readnoise

    def compute(self):
        """Main function

        Parameters
        ----------
        out_dir : str
            Directory in which to save output file
        """
        # Read in dark current ramps
        #print('Working on file: {}'.format(self.dark))
        with fits.open(self.dark) as hdu:
            dark_ramp = hdu[1].data
            self.instrument = hdu[0].header['INSTRUME']
            self.detector = hdu[0].header['DETECTOR']

        # Create a series of CDS images
        cds = self.make_cds(dark_ramp)

        # Create grid of boxes across detector
        xs, ys = self.box_grid()

        # Calcualte the stdev in each box through the stack of CDS images
        # Create an output frame where all pixels within the box are
        # set equal to this stdev value
        readnoise_frame = self.calculate_stdev(cds, xs, ys, self.sigma_threshold)

        # Save the readnoise frame in the format of a readnoise reference file
        self.save(readnoise_frame)

    def make_cds(self, ramp):
        """Given an integration (or exposure with multiple
        integrations) create a series of CDS images. Assume
        the input comes from a DMS-formatted file, in which
        the data are in a 4D array (integration, group, y, x)

        Parameters
        ----------
        ramp1 : numpy.ndarray
            4-dimensional array of dark signals

        Returns
        -------
        cds_ramp : numpy.ndarray
            3-dimensional cube of difference images
        """
        num_integrations, num_groups, num_y, num_x = ramp.shape

        for integration in range(num_integrations):
            integ_cds = ramp[integration, 1:, :, :] - ramp[integration, 0:-1, :, :]
            if integration == 0:
                cds_ramp = copy.deepcopy(integ_cds)
            else:
                cds_ramp = np.concatenate(cds_ramp, integ_cds)
        return cds_ramp

    def save(self, data):
        """Save the readnoise frame using JWST datamodel"""
        ron = ReadnoiseModel()
        ron.data = data
        ron.meta.instrument.name = self.instrument
        ron.meta.instrument.detector = self.detector
        ron.meta.subarray.name = 'GENERIC'
        ron.meta.subarray.xstart = 1
        ron.meta.subarray.ystart = 1
        ron.meta.subarray.xsize = 2048
        ron.meta.subarray.ysize = 2048
        ron.meta.subarray.fastaxis = 1
        ron.meta.subarray.slowaxis = -2
        ron.meta.exposure.readpatt = 'ANY'
        ron.meta.description = 'CDS Noise Image'
        ron.meta.pedigree = 'GROUND'
        ron.meta.useafter = '2115-10-01T00:00:00'
        ron.meta.author = 'NOBODY'
        ron.save(self.output_name)

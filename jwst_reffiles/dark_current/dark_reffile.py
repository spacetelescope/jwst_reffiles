#! /usr/env/bin python

"""
This module contains code for creating a dark current reference file.

Input:
    List of exposures at a calibration state from which you wish to create
    dark current reference files. For NIR detectors, this should be after:

    0) group_scale
    1) dq_init
    2) saturation
    3) ipc* - not currently run
    4) superbias
    5) refpix
    6) linearity
    7) persistence

    For MIRI, this should be after:

    0) group_scale
    1) dq_init
    2) saturation
    3) linearity
    4) rscd

    Also, versions of these files that have been run up through the JUMP
    step must also be present. This is because we need to know where the
    CR hits are in the data.

    All files should be from the same instrument/detector/aperture/readpattern



Ouputs:
    A dark current reference file in the proper format
    A series of images showing for each pixel in each group, how many
        measurements of the signal were used to calculate the mean
        dark value.

Function:
Get file headers from each file. Check for consistency.
Extract CR flags from the JUMP version of the files.
    From these, create arrays giving the group number of the first cosmic
    ray hit in each pixel for each integration
Create images of how many dark current measurements are considered good
    for each pixel and each group. All groups after the first CR hit in a
    pixel are considered bad and not used.
Read in the data. To avoid memory problems only read in N groups at a time,
    where N is determined from user inputs.
Flag all groups that occur between the first CR hit and the end of the integration
    as bad.
From the remaining good data, calculate the sigma-clipped mean dark signal
    in each pixel for each group (i.e. take the mean across integrations).
    Also calculate the sigma-clipped stdard deviation.
Save the sigma-clipped mean as the dark data, and the sigma-clipped stdev
    as the uncertainty in a DarkModel or MIRIDarkModel instance. The
    DQ array is currently all zeros, anticipating bad pixels being placed
    in the bad pixel reference file.
"""

from astropy.io import fits
from astropy.stats import sigma_clip
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

from jwst.datamodels import DarkModel, DarkMIRIModel

from jwst_reffiles.utils import dq_flags
from jwst_reffiles.utils.file_utils import read_ramp, get_file_header


class Dark():
    def __init__(self, file_list=[], max_equivalent_groups=60, sigma_threshold=3, reffile_metadata=None,
                 output_file=None, contribution_file='good_signals_per_group.pdf'):
        self.file_list = file_list
        self.sigma_threshold = sigma_threshold
        self.reffile_metadata = reffile_metadata  # Dict?
        self.output_file = output_file
        self.contribution_file = contribution_file

        # Check that all necessary entries are in the dictionary to be used
        # to populate the header of the final reference file
        #self.check_reffile_info()

        # Get metadata from first input dark.
        self.metadata = get_file_header(self.file_list[0], extension='PRIMARY')

        # Use array size and number of
        # groups as a guide to determine if and how to split up the data
        # into chunks, such that one chunk is read in from all input files
        # and calculations performed without using too much memory
        if self.metadata['INSTRUME'] != 'MIRI':
            x_fullframe_dim = 1024
            y_fullframe_dim = 1024
        else:
            x_fullframe_dim = 2048
            y_fullframe_dim = 2048
        pix_per_fullframe_group = x_fullframe_dim * y_fullframe_dim

        pix_per_group = self.metadata['SUBSIZE1'] * self.metadata['SUBSIZE2']
        pix_ratio = np.int(pix_per_fullframe_group / pix_per_group)

        # The number of groups that can be read in from each file is the
        # user-entered max number of groups, divided by the ratio of the
        # number of pixels in a full frame exposure to the number of pixels
        # in the input aperture, divided by the number of files, and
        # finally divided by 2, since we will also need to read in the
        # cosmic ray flags associated with each file
        delta_groups = np.int(max_equivalent_groups / pix_ratio / len(self.file_list) / 2)
        self.group_index = np.arange(0, self.metadata['NGROUPS'], delta_groups)
        self.group_index = np.append(self.group_index, self.metadata['NGROUPS'])

    def consistency_check(self):
        """Make sure that all input files have the same dimensions,
        number of groups, etc
        """
        keywords_to_check = ['INSTRUME', 'DETECTOR', 'SUBARRAY', 'READPATT', 'SUBSIZE1',
                             'SUBSIZE2', 'NGROUPS']
        for filename in self.file_list[1:]:
            metadata = get_file_header(filename)
            for kw in keywords_to_check:
                if metadata[kw] != self.metadata[kw]:
                    raise ValueError(("ERROR: inconsistent input files. The first input "
                                      "file contains: \n{}: {}, while {} contains: \n{}: {}"
                                      .format(kw, self.metadata[kw], filename, kw, metadata[kw])))

    def create_num_contributions_plots(self, clean_signals, bad_signals):
        """Assume PDF
        """
        # Make sure the output file for the plots is a PDF
        if self.contribution_file[-4:].lower() != '.pdf':
            self.contribution_file = '{}.pdf'.format(self.contribution_file)

        # Save all plots to a single PDF
        pdf = PdfPages(self.contribution_file)
        for group in range(ngroups):
            fig = plt.figure(figsize=(9, 9))
            ax1 = fig.add_subplot(1, 1, 1)
            ngroups, ysize, xsize = clean_signals.shape
            im = ax1.imshow(clean_signals[group, :, :], extent=[0, xsize, 0, ysize], interpolation='None',
                            cmap=cm.RdYlGn, origin='lower', vmin=0, vmax=image_max)
            plt.colorbar(im)

            ax1.set_title('Signal Values Used to Create Mean Dark, Group {}'.format(group))
            fig.tight_layout()
            pdf.savefig(fig)
        pdf.close()

    def create_reffile(self):
        """MAIN FUNCTION
        Creates a dark current reference file from a list of input dark
        current exposures
        """
        # Data consistency check
        self.consistency_check()

        # Read in the DQ arrays from the Jump files and create a map of
        # group numbers associated with the first CR hit in the integration
        cr_indexes = self.prepare_cr_maps()

        # Collect info on the number of remaining good pixels
        for filename in filenames:
            # Number of good groups per integration
            num_good = copy.deepcopy(cr_indexes[filename])

            # Set pixels with no hits as good in all groups
            num_good[np.isnan(num_good)] = self.metadata['NGROUPS']

            # Number of bad groups per integration
            num_bad = self.metadata['NGROUPS'] - num_good

            # Translate into the number of good signal values
            # for each pixel in each group (i.e. look across
            # integrations)
            good_counter, bad_counter = dq_flags.signals_per_group(num_good, self.metadata['NGROUPS'])

            # Save num_good and num_bad as images
            #do_that_here()

            # What we really want is an image of the number of good signal
            # values for every pixel in every group
            self.create_num_contributions_plots(good_counter, bad_counter)

        # Get a map of reference pixels
        refpix_map = dq_flags.create_refpix_map(filenames[0])

        # Create arrays to hold the final dark current data and uncertainty
        final_dark = np.zeros((self.metadata['NGROUPS'], self.metadata['SUBSIZE2'], self.metadata['SUBSIZE1']))
        final_dark_dev = np.zeros((self.metadata['NGROUPS'], self.metadata['SUBSIZE2'], self.metadata['SUBSIZE1']))
        # Read in the data
        for i in range(len(self.group_index) - 1):
            for file_index, filename in enumerate(self.file_list):

                # Get metadata so we know how many integrations are in
                # the file, so that we can read them all in.
                metadata = get_file_header(self.file_list[0], extension='PRIMARY')
                for int_index, integration in enumerate(metadata['NINTS']):

                    # Read in the dark data
                    file_data = read_ramp(filename, integration, self.group_index[i],
                                          self.group_index[i + 1], extension='SCI')

                    jump_map = cr_indexes[filename][integration, :, :]

                    # Pixels/groups that come after a CR hit are set to NaN
                    group_arr = np.zeros_like(file_data)
                    for idx, g in enumerate(groups):
                        group_arr[idx, :, :] = g
                    diff = group_arr - jump_map
                    hits = np.where(diff >= 0)
                    file_data[hits[0], hits[1], hits[2]] = np.nan

                    # Place data in full container
                    if int_index == 0:
                        data = np.zeros_like(file_data)
                        data[:, :, :, :] = file_data
                    else:
                        data = np.vstack(data, file_data)

            # Now calculate the sigma_clipped mean for each pixel through
            # all the files
            clipped_data = sigma_clip(data, sigma=self.sigma_threshold, axis=0, masked=False)
            mean_dark = np.nanmean(clipped_data, axis=0)
            dev_dark = np.nanstd(clipped_data, axis=0)

            # Place the mean and stdev of the dark signals into the correct
            # location within the final dark arrays that contain all groups
            final_dark[self.group_index[i]: self.group_index[i + 1], :, :] = mean_dark
            final_dark_dev[self.group_index[i]: self.group_index[i + 1], :, :] = dev_dark

            # DQ array associated with dark data
            # All bad pix are put in bad pixel mask, so dark reffile DQ array
            # is empty?
            dq_dark = np.zeros_like(final_dark).astype(np.int8)

            # Save reference file
            self.save_reffile(fianl_dark, final_dark_dev, dq_dark)

    def prepare_cr_maps(self):
        """Read in groupdq extension from jump files, extract only the
        JUMP_DET flags, and collapse into an array that lists the group
        number of the first CR hit for each pixel in each integration

        Returns
        -------
        hit_indexes : dict
            key values are the filenames associated with the data. Each
            value is a 3D array (integ, y, x) of group numbers of CR hits
        """
        hit_indexes = {}
        self.jump_files = []
        for filename in self.file_list:
            # Find the jump file (i.e. output from jump step of calwebb_detector1)
            suffix = filename.split('_')[-1]
            jump_file = filename.replace(suffix, 'jump.fits')
            self.jump_files.append(jump_files)

            if not os.path.isfile(jump_file):
                raise FileNotFoundError("ERROR: Jump file {} not found.".format(jump_file))

            print('Getting CR map from Jump File {}'.format(jump_file))

            # Read in the entire exposure from the GROUPDQ extension
            groupdq = read_ramp(jump_file, extension='GROUPDQ')

            # Keep only the JUMP_DET flags
            file_cr_map = dq_flags.flag_map(groupdq, 'JUMP_DET')

            # Translate CR maps per group into one array per integration.
            # The pixel value is the group number of the first CR hit
            first_hits = self.collapse_cr_map(file_cr_map)

            # Update the CR index map for this file
            hit_indexes[filename] = first_hits
        return hit_indexes

    def save_reffile(self, data, uncertainty, dq):
        """Save the given data in the official reference file format

        Parameters
        ----------
        data : numpy.ndarray
            3D array containing the mean dark current signals

        uncertainty : numpy.ndarray
            3D array containing associated uncertainties

        dq : numpy.ndarray
            2D array containing associated data quality flags
        """
        if self.metadata['instrument'].upper() == 'MIRI':
            model = DarkMIRIModel()
        else:
            model = DarkModel()

        # Insert data
        model.data = data
        model.err = uncertainty
        model.dq = dq

        # Metadata copied from input files
        model.meta.reftype = 'DARK'
        model.meta.instrument.name = self.metadata['INSTRUME'].upper()
        model.meta.instrument.detector = self.metadata['DETECTOR'].upper()
        model.meta.exposure.type = self.metadata['EXPTYPE'].upper()
        model.meta.exposure.readpatt = self.metadata['READPATT'].upper()
        model.meta.exposure.nframes = self.metadata['NFRAMES']
        model.meta.exposure.ngroups = self.metadata['NGROUPS']
        model.meta.exposure.groupgap = self.metadata['GROUPGAP']
        model.meta.subarray.name = self.metadata['SUBARRAY']
        model.meta.subarray.xstart = self.metadata['SUBSTRT1']
        model.meta.subarray.ystart = self.metadata['SUBSTRT2']
        model.meta.subarray.xsize = self.metadata['SUBSIZE1']
        model.meta.subarray.ysize = self.metadata['SUBSIZE2']
        model.meta.subarray.fastaxis = self.metadata['FASTAXIS']
        model.meta.subarray.slowaxis = self.metadata['SLOWAXIS']

        # Instrument-specific metadata
        if self.metadata['instrument'].upper() != 'MIRI':
            model.meta.exposure.gain_factor = self.metadata['GAINFACT']
        else:
            model.meta.exposure.preadpatt = self.metadata['P_READPA']

        # User-provided metadata
        model.meta.pedigree = self.reffile_metadata['pedigree']
        model.meta.descrip = self.reffile_metadata['descrip']
        model.meta.author = self.reffile_metadata['author']
        model.meta.useafter = self.reffile_metadata['useafter']

        # Add HISTORY entries
        if 'history' in self.reffile_metadata.keys():
            model.history.append(self.reffile_metadata['history'])
        if 'HISTORY' in self.reffile_metadata.keys():
            model.history.append(self.reffile_metadata['HISTORY'])
        model.history.append(("{} {} {} dark current file, created from data in files: "
                              .format(self.metadata['INSTRUME'].upper(), self.metadata['DETECTOR'].upper(),
                                      self.metadata['SUBARRAY'])))
        for filename, jump_file in zip(self.file_list, self.jump_files):
            model.history.append(filename)
            model.history.append(jump_file)

        # Create output name if necessary
        if self.output_file is None:
            self.output_file = '{}_{}_{}_dark_ref_file.fits'.format(self.metadata['INSTRUME'],
                                                                    self.metadata['DETECTOR'],
                                                                    self.metadata['SUBARRAY'])
        model.save(self.output_file)

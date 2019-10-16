#! /usr/env/bin python

"""
This module contains code for creating a dark current reference file
"""

from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np

from jwst.datamodels import DarkModel, DarkMIRIModel


class Dark():
    def __init__(self, file_list=[], max_equivalent_groups=60, sigma_threshold=3, reffile_header_info=None,
                 output_file=None):
        self.file_list = file_list
        self.sigma_threshold = sigma_threshold
        self.reffile_header_info = reffile_header_info  # Dict?
        self.output_file = output_file

        # Check that all necessary entries are in the dictionary to be used
        # to populate the header of the final reference file
        #self.check_reffile_info()

        # Get metadata from first input dark.
        self.metadata = utils.get_file_info(self.file_list[0])

        # Use array size and number of
        # groups as a guide to determine if and how to split up the data
        # into chunks, such that one chunk is read in from all input files
        # and calculations performed without using too much memory
        if self.metadata['instrument'] != 'MIRI':
            x_fullframe_dim = 1024
            y_fullframe_dim = 1024
        else:
            x_fullframe_dim = 2048
            y_fullframe_dim = 2048
        pix_per_fullframe_group = x_fullframe_dim * y_fullframe_dim

        pix_per_group = metadata['xdim'] * metadata['ydim']
        pix_ratio = np.int(pix_per_fullframe_group / pix_per_group)

        # The number of groups that can be read in from each file is the
        # user-entered max number of groups, divided by the ratio of the
        # number of pixels in a full frame exposure to the number of pixels
        # in the input aperture, divided by the number of files, and
        # finally divided by 2, since we will also need to read in the
        # cosmic ray flags associated with each file
        delta_groups = np.int(max_equivalent_groups / pix_ratio / len(self.file_list) / 2)
        self.group_index = np.arange(0, self.metadata['ngroup'], delta_groups)
        self.group_index = np.append(self.group_index, self.metadata['ngroup'])

    def consistency_check(self):
        """Make sure that all input files have the same dimensions,
        number of groups, etc
        """
        for filename in self.file_list[1:]:
            metadata = utils.get_file_info(filename)
            if metadata != self.metadata:
                raise ValueError(("ERROR: inconsistent input files. The first input "
                                  "file is:\n {}\n While {} is {}.".format(self.metadata,
                                                                           filename,
                                                                           metadata)))

    def create_reffile(self):
        """
        """
        # Data consistency check
        self.consistency_check()

        # Read in the DQ arrays from the Jump files
        cr_indexes = self.prepare_cr_maps()

        # Collect info on the number of remaining good pixels
        for filename in filenames:
            # Number of good groups per integration
            num_good = self.metadata['ngroup'] - cr_indexes[filename]

            # Set pixels with no hits as good in all groups
            num_good[np.isnan(num_good)] = self.metadata['ngroup']

            # Number of bad groups per integration
            num_bad = self.metadata['ngroup'] - num_good

            # Save num_good and num_bad as images
            do_that_here()

            # What we really want is an image of the number of good signal
            # values for every pixel in every group
            do_that_here()

        # Read in the data
        for i in range(len(self.group_index) - 1):
            for file_index, filename in enumerate(self.file_list):

                # Get metadata so we know how many integrations are in
                # the file, so that we can read them all in.
                metadata = utils.get_file_info(self.file_list[0])
                for int_index, integration in enumerate(metadata['nint']):

                    # Read in the dark data
                    file_data = read_ramp(filename, integration, self.group_index[i],
                                          self.group_index[i + 1], extension='SCI')

                    jump_map = cr_indexes[filename][integration, :, :]
                    #groups = self.group_index[i], self.group_index[i + 1]

                    #file_data = arr[3, 2048, 2048]
                    #groups = [min, max]
                    #jump_map = arr[2048, 2048]

                    #for pix1 in x:
                    #    for pix2 in y:
                    #        ramp = file_data[:, pix1, pix2]
                    #        jumps = jump_map[pix1, pix2]
                    #        groups = np.arange(self.group_index[i], self.group_index[i + 1])

                    #        groups = [3, 4, 5, 6]
                    #        jumps = 4
                    #        ramp = [d1, d2, d3, d4]
                    #        hit = np.where(groups >= jumps)
                    #        ramp[hit] = np.nan

                    # Same as above but without looping over pixels
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

            # DQ array associated with dark data
            # All bad pix are put in bad pixel mask, so dark reffile DQ array
            # is empty?
            dq_dark = np.zeros_like(mean_dark).astype(np.int8)

            # Save reference file
            self.save_reffile(mean_dark, dev_dark, dq_dark)








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
        for filename in self.file_list:
            # Find the jump file (i.e. output from jump step of calwebb_detector1)
            suffix = filename.split('_')[-1]
            jump_file = filename.replace(suffix, 'jump.fits')

            if not os.path.isfile(jump_file):
                raise FileNotFoundError("ERROR: Jump file {} not found.".format(jump_file))

            print('Getting CR map from Jump File {}'.format(jump_file))

            # Read in the entire exposure from the GROUPDQ extension
            groupdq = read_full_ramp(jump_file, extension='GROUPDQ')

            # Keep only the JUMP_DET flags
            file_cr_map = dq_flags.flag_map(groupdq, 'JUMP_DET')

            Flag map above may need to be tweaked to handle multiple groups

            pixels with no cr hits should be set to np.nan

            # Translate CR maps per group into one array per integraiton.
            # The pixel value is the group number of the first CR hit
            hit_indexes = self.collapse_cr_map(file_cr_map)

            # Update the CR index map for this file
            hit_indexes[filename] = cr_indexes
        return hit_indexes







    def read_ramp(filenames):
    """Read in the science extension from a group of slope images

    Parameters
    ----------
    filenames : list
        List of fits files containing slope values

    Returns
    -------
    instrument : str
        Name of instrument associated with the data

    slope_data : numpy.ndarray
        3D array containing slope values for science pixels only.
        Reference pixels have been stripped off.

    starting_indexes : list
        List of numbers corresponding to the index numbers within slope_data
        where each file's data begins. For example, if slope_data is an array
        of size (10, 2048, 2048), and starting_indexes = [0, 5, 7, 10] then
        we can pull apart slope_data into its constituent exposures using
        slope_data[starting_indexes[0]: starting_indexes[1]], etc.

    left_cols : int
        Number of columns of reference pixels on the left side of the array

    right_cols : int
        Number of columns of reference pixels on the right side of the array

    bottom_rows : int
        Number of rows of reference pixels on the bottom of the array

    top_rows : int
        Number of rows of reference pixels on the top of the array
    """
    slope_stack = []
    starting_indexes = []
    for i, slope_file in enumerate(filenames):
        if not os.path.isfile(slope_file):
            raise FileNotFoundError('ERROR: Input slope file {} does not exist'.format(slope_file))

        with fits.open(slope_file) as hdulist:
            slope_img = hdulist['SCI'].data
            dq_int = hdulist['DQ'].data
            header = hdulist[0].header
            instrument = header['INSTRUME']
            slope_shape = slope_img.shape
            if len(slope_shape) == 2:
                dq_img = (dq_int[:, :] & dqflags.pixel['REFERENCE_PIXEL'] == 0)
            elif len(slope_shape) == 3:
                dq_img = (dq_int[0, :, :] & dqflags.pixel['REFERENCE_PIXEL'] == 0)
            else:
                raise ValueError("Slope image should be either 2D or 3D.")

        # Create a mask where 1 indicates a science pixel and 0 indicates
        # a reference pixel

        science = np.where(dq_img == 1)
        left_edge = np.min(science[1])
        right_edge = np.max(science[1]) + 1
        bottom_edge = np.min(science[0])
        top_edge = np.max(science[0]) + 1

        left_cols = left_edge
        right_cols = dq_img.shape[1] - right_edge
        bottom_rows = bottom_edge
        top_rows = dq_img.shape[0] - top_edge

        # Add to the list of starting indexes
        starting_indexes.append(len(slope_stack))

        # loop over integrations and pull out slope for int
        # Crop the reference pixels from the array.

        if len(slope_shape) == 2:
            slopes = slope_img[bottom_edge:top_edge, left_edge:right_edge]
            slope_stack.append(slopes)
        elif len(slope_shape) == 3:
            num_int = slope_shape[0]
            for i in range(num_int):
                slopes = slope_img[i, bottom_edge:top_edge, left_edge:right_edge]
                slope_stack.append(slopes)
    starting_indexes.append(len(slope_stack))
    slope_data = np.array(slope_stack)
    return instrument, slope_data, starting_indexes, (left_cols, right_cols, bottom_rows, top_rows)

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

        model.data = data
        model.err = uncertainty
        model.dq = dq

        model.meta.instrument.name = self.metadata['instrument'].upper()
        model.meta.otherstuff = self.reffile_header_info['otherstuff']

        if self.output_file is None:
            self.output_file = '{}_{}_dark_ref_file.fits'.format(self.metadata['instrument'],
                                                                 self.metadata['aperture'])
        model.save(self.output_file)

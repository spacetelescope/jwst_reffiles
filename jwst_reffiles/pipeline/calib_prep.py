#! /usr/bin/env python

'''
Organize science files to be used by mkrefs.py to create reference files.
Convert a list of input files and required pipeline steps into
a list of strun commands that can be farmed out to condor.

Inputs:

inputs -
         Astropy table of filenames, reference file types, and required
         calibration pipeline steps.

input table example:

index    reftype     fitsfile             ssbsteps
1         gain      something_uncal.fits     bpm,sat
2         gain      another_uncal.fits       bpm,sat
3        readnoise  something_uncal.fits     bpm,sat,superbias


colnames for input table, as defined in mkrefs.py
['cmdID', 'reftype', 'imlabel', 'imtype', 'detector', 'ssbsteps', 'imindex', 'imID', 'MJD', 'fitsfile']



search_dir -
             Directory to search for science files (from the 'filenames'
             column of the input table. Search looks into any
             subdirectories as well, to allow the user flexibility when
             organizing input files

output_dir -
             Directory in which the pipeline-processed science files will
             be placed. This code does not actually run the pipeline, but
             only generates the commands needed to run the pipeline. The
             output_dir will be added to these commands.

use_only_given -
                 Use exactly the file listed in the 'filenames'
                 column of the input table. If this is False,
                 the code will search the search_dir to see
                 what versions of the input file are present,
                 and choose the most appropriate version for the
                 requested pipeline processing state

overwrite_existing_files -
                           If True, any existing files in the output
                           directory that have names matching the output
                           names will be over-written. This is essentially
                           a way to force the calibration pipeline to be
                           re-run on a given input file

Outputs:

proc_table -
             Copy of the input table with several columns added.

             real_input_file:
             This added column lists the name of the file to be run through
             the pipeline. This may differ from the filename originally given
             as input, in the case where a version of the input file that has
             already been partially processed is found and will be used to
             save time.

             output_name:
             This column lists the name of the output file that will be
             produced by the pipeline

             steps_to_run:
             The column lists the pipeline steps that need to be run in order
             to transform the real_input_file into the output_name file.


strun -
        This is a list of calls to the calibration pipeline. In order to use
        condor, we must use the command-line, strun, calls to the pipeline.


HISTORY:

Version 0.1 - Created by Hilbert, 11 Dec 2017.
'''

import copy
from collections import OrderedDict
from glob import glob
import itertools
import os
import re
import subprocess
import sys
import time

from astropy.io import ascii, fits
from astropy.table import Column, Table
import numpy as np
from jwst import datamodels

from jwst_reffiles.pipeline.pipeline_steps import get_pipeline_steps
from jwst_reffiles.utils.definitions import PIPE_STEPS, PIPE_KEYWORDS, INSTRUMENTS


class CalibPrep:
    def __init__(self, instrument):
        self.__version__ = 0.1
        self.verbose = True
        self.inputs = None
        self.search_dir = ''
        self.use_only_given = False
        self.overwrite_existing_files = True
        self.output_dir = ''
        #self.pipe_step_dict = OrderedDict(PIPE_STEPS)

        instrument = instrument.lower()
        if instrument not in INSTRUMENTS:
            raise ValueError("WARNING: {} is not a valid instrument name.".format(instrument))
        self.pipe_step_list = get_pipeline_steps(instrument)

    def activate_encompassing_entries(self):
        """Once the search for ssb commands contained within other ssb commands has
        been completed, go through the list of sets and activate the top level rows
        that contain other commands. Remember that there can be nested values here.
        (e.g. row 4 contains [0], row 3 contains [0,4], row 5 contains [0,3,4]) In this
        example, we want to set the contained_within value to -1 for row 5, as that
        single run of the pipeline will output products needed by rows 0, 3, and 4 also.

        Rows that do not contain the command for any other rows will have an empty
        set. These need to be set to -1 as well, so that they are run.
        """
        empty_entries = self.inputs['index_contained_within'] == set([])
        self.inputs['index_contained_within'][empty_entries] = set([-1])

        # Change the entries in the 'index_contained_within' to be lists, for easy
        # indexing later.
        new_contained_data = [sorted(list(entry)) for entry in self.inputs['index_contained_within']]

        # To ensure that each entry in the Column is a list (rather than a scalar)
        # add a dummy value that is a 2-element list.
        new_contained_data.append([-42, -43])
        new_contained_column = Column(new_contained_data, name='index_contained_within')

        # Now remove dummy value
        new_contained_column = new_contained_column[0:-1]

        self.inputs.remove_column('index_contained_within')
        self.inputs.add_column(new_contained_column)

    def activate_one_from_repeats(self):
        """Once the repeat search has been completed, activate only one copy of each
        repeated entry.
        """
        repeats = self.inputs['repeat_of_index_number']

        # Convert to a list so we can find unique elements
        repeats_list = [sorted(list(s)) for s in repeats]

        # It seems that when all elements of repeats_list are one element lists,
        #np.unique returns a list, whereas if some of the lists within repeat_list
        #are multiple-element lists, then np.unique returns a list of lists. Try adding
        # this dummy entry to force unique to always return a list of lists
        repeats_list.append([-42, -42])
        unique_repeats = np.unique(repeats_list)

        # Now pop the dummy entry out of unique_repeats
        unique_repeats = unique_repeats[1:]

        for repeat_values in unique_repeats:
            # Set the first instance of the repeated group to -1 so it will be run.
            # Leave the other instances as-is so they won't be run
            self.inputs['repeat_of_index_number'][repeat_values[0]] = set([-1])

        # Change the entries in the 'repeat_of_index_number' to be lists, for easy
        # indexing later.
        new_repeat_data = [sorted(list(entry)) for entry in self.inputs['repeat_of_index_number']]

        # Repeat the dummy value trick here so that the column will contain lists
        new_repeat_data.append([-42, -42])
        new_repeat_column = Column(new_repeat_data, name='repeat_of_index_number')
        new_repeat_column = new_repeat_column[0:-1]

        self.inputs.remove_column('repeat_of_index_number')
        self.inputs.add_column(new_repeat_column)

    def build_search_generator(self):
        '''Create a generator to use for searching for files

        Returns
        -------
        generator : generator
            Generator from os.walk
        '''
        if isinstance(self.search_dir, str):
            # Does the directory exist?
            if not os.path.isdir(self.search_dir):
                print('WARNING: directory {} does not exist!'.format(self.search_dir))
                print('WARNING: nothing found!')
                return([])
            generator = os.walk(self.search_dir, topdown=True)
        elif isinstance(self.search_dir, list):
            # If the search directory is a list, then construct a generator
            # for each directory, and chain them together into a single generator
            generators = []
            for searchdir in self.search_dir:
                # Does the directory exist?
                if not os.path.isdir(searchdir):
                    print('WARNING: directory {} does not exist!'.format(searchdir))
                    continue
                generators.append(os.walk(searchdir, topdown=True))

            # Any entries? If not return None
            if len(generators) == 0:
                print('WARNING: nothing found!')
                return([])

            # Chain together
            generator = itertools.chain()
            for gen in generators:
                generator = itertools.chain(generator, gen)
        return generator

    def choose_file(self, filelist, req, current):
        '''Given a list of files, their calibration state
        and the required calibration state, choose the best file
        with which to proceed

        Paramters
        ---------
        filelist : list
            List of fits filenames

        req : dict
            Dictionary of pipeline steps containing boolean entries to describe
            the required calibration state of the file

        current : dict
            Nested dictionary: for each filename (key), the value is a dictionary
            of pipeline steps containing boolean entires to describe the current
            calibration state of the file

        Returns
        -------
        usefile : str
            Name of the fits file to begin further processing with
        '''
        usefile = None
        first_false = -1
        for file in filelist:
            comp = []
            for key in req:
                comp.append(req[key] == current[file][key])
            comp = np.array(comp).astype('int')
            f = np.where(comp == 0)[0]
            if ((len(f) > 0) and (f[0] > first_false)):
                first_false = f[0]
                usefile = file
            elif len(f) == 0:
                return file
        return usefile

    def combine_repeats(self):
        '''Search for repeated entries in the table input to calib_prep, and remove them.
        In this case, a repeated entry means that the filename is the same, and that the
        ssbsteps are the same OR one set are contained exactly within another (with no
        extra intervening steps). In this case, a single call to strun can be used to
        create the requested data for both cases. The only thing that needs to be done is
        to provide an intermediate output name to capture both versions of the calibrated
        file.'''

        # Create two columns to add to the table. One indicates whether an entry
        # is a duplicate of another row of the table. The other indicates if the command
        # in one row is contained within the command of another row. Both columns list the
        # index number of the row that it is a duplicate of or contained within.
        num_rows = len(self.inputs['real_input_file'])
        repeat_column = Column(data=[set([i]) for i in range(num_rows)], name='repeat_of_index_number')
        self.inputs.add_column(repeat_column)

        within_column = Column(data=[set([]) for i in range(num_rows)], name='index_contained_within')
        self.inputs.add_column(within_column)

        # We also need to make a new column to hold the strun commands. After this step,
        # command lengths can potentially become much longer than they are at the time of
        # input, but the maximum string length for the input column has already been set.
        # So we need to make a replacement column with the updated commands so that the new
        # commands aren't truncated at the previous maximum string length
        # new_commands = [''] * len(self.inputs['strun_command'])
        new_commands = list(self.inputs['strun_command'].data)

        for row in self.inputs:
            # Find all entries that have the same INPUT filename
            repeated_filename = self.inputs['real_input_file'] == row['real_input_file']

            # Don't count the row itself when looking for matches
            current_row_index = row['index']
            repeated_filename[current_row_index] = False

            # If the same starting file is present in more than one row:
            if np.sum(repeated_filename.astype('int')) > 0:
                matching_rows = self.inputs[repeated_filename]

                # Strip all whitespace from the list of ssb steps. Use this as
                # a comparison to see if there are repeats among those with matching
                # filenames.
                row_ssb = re.sub(r'\s+', '', row['steps_to_run'])
                outname = row['output_name']
                for matching_row in matching_rows:
                    ssbsteps = re.sub(r'\s+', '', matching_row['steps_to_run'])
                    row_index = matching_row['index']

                    if ssbsteps in row_ssb:
                        # If the list of ssb steps is contained within the list of
                        # steps from the comparison file, then this entry is duplicate.
                        # Save the output name from this row if it is different from the
                        # comparison row.
                        additional_output = matching_row['output_name']

                        if (additional_output == outname):
                            # Same output filename means the two rows are exact copies
                            self.inputs[row_index]['repeat_of_index_number'].add(current_row_index)
                        elif (additional_output != outname):
                            # Different output names means that to combine these rows we need
                            # to have two outputs. Identify which ssb step the intermediate
                            # output is from, and add it to the strun command
                            self.inputs[row_index]['index_contained_within'].add(current_row_index)
                            final_step = ssbsteps.split(',')[-1]
                            additional_out_str = (" --steps.{}.output_file={}"
                                                  .format(final_step, additional_output))
                            new_command = copy.deepcopy(new_commands[current_row_index]) + additional_out_str
                            new_commands[current_row_index] = new_command
                        else:
                            print('input names match but no repeat nor subcommand. or already flagged')
                    else:
                        print('ssb steps is not in or matching other ssbsteps.')
            else:
                print("No repeat nor sub-command.")

        # Now remove the old strun_command column and insert the new one
        self.inputs.remove_column('strun_command')
        new_command_column = Column(data=new_commands, name='strun_command')

        self.inputs.add_column(new_command_column)

        # Activate only one copy of each matching repeat set
        self.activate_one_from_repeats()

        # Check for subsets in the index_contained_within column, and keep only the entries
        # that contain all others
        self.activate_encompassing_entries()

    def completed_steps(self, input_file):
        '''Identify and return the pipeline steps completed
        for the input file

        Parameters
        ----------
        input_file : str
            Name of input fits file

        Returns
        -------
        completed : dict
        '''
        completed = OrderedDict({})
        for key in self.pipe_step_list:
            completed[key] = False
        header = fits.getheader(input_file)
        for key in header.keys():
            if key in PIPE_KEYWORDS.keys():
                value = header.get(key)
                if value == 'COMPLETE':
                    completed[PIPE_KEYWORDS[key]] = True
        return completed

    def copy_config_to_output_dir(self, filename):
        """Check to see if the given pipeline configuration file exists
        in the output directory. If not, make a copy from the files in
        the repo

        Parameters
        ----------
        filename : str
            Name of the cfg file to check for (e.g. calwebb_detector1.cfg)
        """
        full_filename = os.path.join(self.output_dir, filename)
        if not os.path.isfile(full_filename):
            print('INFO: {} file does not exist in {}. Creating.'.format(filename,
                                                                         self.output_dir))
            cfg_reference = os.path.join(os.path.dirname(__file__), 'config_files/{}'.format(filename))
            subprocess.call(['cp', cfg_reference, filename])

    def create_output(self, base, req):
        '''Create the output name of the pipeline-processed
        input file given the list of required steps

        Parameters
        ----------
        base : str
            Base name of the fits file. Should be the filename minus '.fits'

        req : dict
            Dictionary of pipeline steps, with boolean entries describing whether
            those steps should be completed or not.

        Returns
        -------
        output_filename : str
            Full path + filename of the output once all the required pipeline steps
            have been run. Should be the true base name + the suffix associated with
            the latest pipeline step to be run.

        true_base : str
            Base name of the fits file. Should be the filename minus '.fits' as
            well as any pipeline-generated suffix
        '''
        # base (e.g.) something_dq_init
        # Find the latest pipeline step that needs to be done
        # Keep track of any intemediate skipped steps (we
        # probably want to differentiate between files that have
        # had all pipeline steps done, compared to those which
        # share the latest completed step but may have skipped
        # intervening steps
        step = None
        suffix = 'uncal'
        final_suffix_piece = 'uncal'
        skip = copy.deepcopy(self.pipe_step_list)
        baseend = len(base)
        for key in self.pipe_step_list:
            if req[key]:
                step = key
                suffix = "{}_{}".format(suffix, key)
                final_suffix_piece = key
                skip.remove(key)
                # In the case where the filenamne has multiple pipeline step names attached,
                # walk back until we find the actual basename of the file
                #if self.pipe_step_dict[key] in base:
                #    idx = base.index(self.pipe_step_dict[key])
                #    if ((idx < baseend) and (idx != -1)):
                #        baseend = copy.deepcopy(idx) - 1

        # If steps were added, then remove the uncal from the suffix
        if suffix != 'uncal':
            suffix = suffix.replace('uncal_', '')

        # Suffix automatically added by the pipeline
        pipeline_suffix = ''
        if 'rate' not in skip:
            pipeline_suffix = '_rate'

        # Remove the entries in skip that are after the last
        # required pipeline step
        stepvals = np.array(self.pipe_step_list)
        if final_suffix_piece in stepvals:
            lastmatch = np.where(stepvals == final_suffix_piece)[0][0]
        elif final_suffix_piece == 'uncal':
            lastmatch = -1
        else:
            raise IndexError("No entry {} in pipeline step dictionary.".format(final_suffix_piece))
        for sval in stepvals[lastmatch+1:]:
            skip.remove(sval)

        # Find the true basename of the file
        true_base = base[0:baseend]

        # Remove 'uncal' if it's in the file base
        true_base = true_base.replace('_uncal', '')
        if self.verbose:
            print("True base name is {}".format(true_base))
            print("Suffix is {}".format(suffix))
            print("Skip vals are {}".format(skip))
        true_base = os.path.split(true_base)[1]

        # Create the output filename by adding the name of the latest
        # pipeline step to be run
        ofile = true_base + '_' + suffix + pipeline_suffix
        ofile = ofile + '.fits'
        output_filename = os.path.join(self.output_dir, ofile)
        return output_filename, true_base

    def file_search(self, base, generator_object):
        '''Search for all versions of a particular
        file in the given directory

        Parameters
        ----------
        base : str
            Basename of the file to search for. This should be the filename minus
            the .fits as well as any pipeline-created suffix.

        generator_object : obj
            Generator constructed from os.walk on all input search directories

        Returns
        -------
        files : list
            List of found files containing the input base string
        '''
        files = []
        for dirpath, dirnames, fnames in generator_object:
            mch = [f for f in fnames if base in os.path.join(dirpath, f) and f[-4:] == 'fits']

            for m in mch:
                files.append(os.path.join(dirpath, m))
        return files

    def ORIGINAL_file_search_DELME_PROBABLY(self, base, dir):
        '''Search for all versions of a particular
        file in the given directory

        Parameters
        ----------
        base : str
            Basename of the file to search for. This should be the filename minus
            the .fits as well as any pipeline-created suffix.

        dir : str
            Directory to search. The search will look down the directory tree
            into all subdirectories below this level.

        Returns
        -------
        files : list
            List of found files containing the input base string
        '''
        # This version searches all subdirectories also...
        files = []
        for dirpath, dirnames, fnames in os.walk(dir, topdown=True):
            mch = [f for f in fnames if base in os.path.join(dirpath, f)]
            # or could try:
            # fnmatch.filter(fnames,base)
            for m in mch:
                files.append(os.path.join(dirpath, m))
        return files

    def get_file_basename(self, input_file):
        '''Determine a given file's basename.
        Currently this is the full name minus
        the .fits at the end

        Parameters
        ----------
        input_file : str
            Input filename

        Returns
        -------
        base : str
            Base of the filename
        '''
        df = input_file.rfind('.fits')
        if df != -1:
            base = input_file[0:df]
        else:
            raise ValueError(("WARNING: {} is not a valid fits file name."
                              .format(input_file)))
        return base

    def output_exist_check(self, filename):
        '''Check to see if the given file already
        exists. Remove if the user allows it.

        Parameters
        ----------
        filename : str
            Name of file to look for
        '''
        if os.path.isfile(filename):
            if self.overwrite_existing_files:
                try:
                    os.remove(filename)
                except (FileNotFoundError, PermissionError) as error:
                    print(error)
            elif ((self.overwrite_existing_files is False) and
                  (self.output_dir == os.path.dirname(filename))):
                raise ValueError(("WARNING: Output file {} already exists, "
                                  "and overwrite_existing_files set to False."
                                  .format(filename)))

    def prepare(self):
        '''Main function'''

        # Column of output names to add to table
        outfiles = []
        self.strun = []
        realinput = []
        all_to_run = []

        # Create generators of files in self.search_dir, to be searched later for
        # matching filenames
        self.search_dir_decompose()
        search_generator = self.build_search_generator()

        for line in self.inputs:
            file = line['fitsfile']

            req_steps = self.step_dict(line['ssbsteps'])

            if self.verbose:
                print("")
                print("File: {}".format(file))
                print("Input required steps: {}".format(line['ssbsteps']))
                print("Required steps: {}".format(req_steps))

            # In order to search for other
            # versions of the file we need
            # to know the basename
            basename = self.get_file_basename(file)

            if self.verbose:
                print("Basename: {}".format(basename))

            # Create the output filename based on
            # the required pipeline steps
            outname, true_base = self.create_output(basename, req_steps)
            outfiles.append(outname)

            if self.verbose:
                print("Output name: {}".format(outname))

            # Check to see if a file with the same name as the
            # output filename already exists in the output directory.
            # If so and the user allows it, remove the file. If the user
            # has not allowed the removal of files, throw an error.
            # Similarly, if permissions prevent you from successfully
            # removing the file, throw an error
            # AR: I remove this! We don't want to remove files that already exist! Only if --force_redo_ssb
            #self.output_exist_check(os.path.join(self.output_dir, outname))

            # Search the self.search_dir directories for partially
            # processed files.
            if not self.use_only_given:
                files = self.file_search(true_base, search_generator)
            else:
                files = [os.path.join(self.search_dir, file)]

            if len(files) == 0:
                print("No matching files found in {} for {}".format(self.search_dir, true_base))
                print("Falling back to the input file: {}".format(file))
                files = [file]

            if self.verbose:
                print("Found files: {}".format(files))

            # Determine the completed calibration steps
            # for each file
            current_state = {}
            for f in files:
                state = self.completed_steps(f)
                current_state[f] = state

                if self.verbose:
                    print("    file {}, current state: {}".format(f, current_state[f]))

            # Select the file to use
            if self.use_only_given:
                input_file = file
            else:
                input_file = self.choose_file(files, req_steps, current_state)

            realinput.append(input_file)
            if self.verbose:
                print("File to use: {}".format(input_file))

            # Create list of pipeline steps that must be
            # run on the file
            to_run = self.steps_to_run(input_file, req_steps,
                                       current_state[input_file])

            # Add column to input table
            tr = [key for key in to_run if to_run[key]]
            trstr = ''
            for s in tr:
                trstr = trstr + ',' + s
            if len(trstr) > 0:
                trstr = trstr[1:]
            else:
                trstr = 'None'
            all_to_run.append(trstr)

            if self.verbose:
                print("Steps that need to be run: {}".format(to_run))

        # Add the output filename column to the input table
        realcol = Column(data=realinput, name='real_input_file')
        self.inputs.add_column(realcol)
        outcol = Column(data=outfiles, name='output_name')
        self.inputs.add_column(outcol)
        toruncol = Column(all_to_run, name='steps_to_run')
        self.inputs.add_column(toruncol)

        # Add an overall index column (to be used as a reference when flagging
        # repeated rows)
        indexcol = Column(data=np.arange(len(realinput)), name='index')
        self.inputs.add_column(indexcol, index=0)

        # Turn the table into a series of strun commands
        self.strun = self.strun_command(realcol, toruncol, outcol)  # ,reffiles=??)

        # Add the list of strun commands to the input table
        cmds = Column(data=self.strun, name='strun_command')
        self.inputs.add_column(cmds)

        # Find any groups of strun commands that can be combined
        self.combine_repeats()

        self.proc_table = copy.deepcopy(self.inputs)

        if self.verbose:
            print("Input table updated.")

        # Done
        # calling function can now use self.proc_table and
        # self.strun to access results

    def search_dir_decompose(self):
        '''Turn the self.search_dir into a list if necessary'''
        if ',' in self.search_dir:
            self.search_dir = [element.strip() for element in self.search_dir.split(',')]

    def step_dict(self, stepstr):
        '''Translate the input list of required pipeline steps
        into a dictionary with boolean entries

        Parameters
        ----------
        stepstr : str
            Comma-separated list of pipeline step abbreviations. These come from
            the origina config file.

        Returns
        -------
        req : dict
            Dictionary of pipeline steps with boolean values designating whether
            each step was requested in stepstr
        '''
        # Make a copy of pipeline step dictionary and
        # intialize all steps to False
        req = OrderedDict({})
        for key in self.pipe_step_list:
            req[key] = False
        if stepstr is None:
            return req

        # Strip all whitespace from the list of steps
        # and split into a list
        stepslist = [element.strip() for element in stepstr.split(',')]
        for ele in stepslist:
            if ele not in self.pipe_step_list:
                raise ValueError(("WARNING: unrecognized pipeline step: {}"
                                  .format(ele)))
            try:
                req[ele] = True
            except KeyError as error:
                print(error)
        return req

    def steps_to_run(self, infile, req, current):
        '''Return a list of the pipeline steps that need
        to be run to bring the input file up to the state
        described by req. Return a dictionary of boolean
        values just like the req and current inputs

        Parameters
        ----------
        infile : str
            Input fits filename

        req : dict
            Dictionary of pipeline steps with boolean entries descibing whether
            each step should be complete on the final output file

        current : dict
            Dictionary of pipeline steps with boolean entries describing whether
            each step has already been completed on infile

        Returns
        -------
        torun : dict
            Dictionary of pipeline steps with boolean entires descrbing which
            required steps must still be run on infile
        '''
        torun = copy.deepcopy(current)
        for key in req:
            if req[key] == current[key]:
                torun[key] = False
            elif ((req[key] is True) & (current[key] is False)):
                torun[key] = True
            elif ((req[key] is False) & (current[key] is True)):
                print(("WARNING: Input file {} has had {} step run, "
                       "but the requirements say that it should not "
                       "be. Need a new input file.".format(infile, key)))
        return torun

    def strun_command(self, input, steps_to_run, outfile_name, overrides=[], instrument='nircam'):
        '''Create the necessary strun command to run the
        appropriate JWST calibration pipeline steps

        Parameters
        ----------
        input : astropy.table.Column
            astropy table column object listing input fits filenames

        steps_to_run : astropy.table.Column
            astropy table column object giving a comma-separated list of pipeline
            steps to run for each file

        outfile_name : astropy.table.Column
            astropy table column object listing pipeline ouput fits file name

        Returns
        -------
        cmds : list
            List of command line calls to the JWST calibration pipeline with appropriate
            flags to run/skip requested steps
        '''

        # Quick fix for ramp fitting, which uses a different
        # output suffix than step name, unlike the other steps

        #step_names = copy.deepcopy(self.pipe_step_dict)
        # Create a dictionary of step names and corresponding step class names
        # In all cases except ramp-fitting, the key and value will be identical
        step_names = OrderedDict({})
        for key in self.pipe_step_list:
            step_names[key] = key
        # Update the ramp-fitting suffix
        step_names['rate'] = 'ramp_fit'

        cmds = []

        # Determine appropriate reference file overrides
        # something

        # Assume the pipeline config files are in self.output_directory. If not,
        # copy from the repo.
        cfg_file_list = copy.deepcopy(PIPE_STEPS)[0: -1] + ['calwebb_detector1']
        cfg_fiile_list = [entry+'.cfg' for entry in cfg_file_list]
        for cfile in cfg_file_list:
            cfg_file = os.path.join(cfile)
            self.copy_config_to_output_dir(cfg_file)
        pipeline_cfg_file = os.path.join(self.output_dir, 'calwebb_detector1.cfg')

        # Create the strun command
        initial = 'strun {} '.format(pipeline_cfg_file)
        for infile, steps, outfile in zip(input, steps_to_run, outfile_name):
            with_file = initial + infile

            if steps == 'None':
                cmd = 'None'
            else:
                skip_text = ''
                out_text = ''
                for key in step_names:

                    # Add skip statements for the steps to be skipped
                    if key not in steps:
                        skip_text += ' --steps.{}.skip=True'.format(step_names[key])

                    # Add override reference files
                    print("Add ability to override reference files here")

                # Get the suffix from the difference between the input and output names
                infile_base = os.path.split(infile)[1]
                infile_base = infile_base.split('.fits')[0]
                if '_uncal' in infile_base:
                    infile_base = infile_base.split('_uncal')[0]
                outfile_base = os.path.split(outfile)[1]
                outfile_base = outfile_base.split('.fits')[0]
                idx = len(infile_base)
                added_suffix = outfile_base[idx+1:]

                # Add step-specific options for saving
                finstep = steps.split(',')[-1]
                final_step = step_names[finstep]
                out_text += (' --steps.{}.suffix={} --steps.{}.output_dir={} --steps.{}.save_results=True'
                             .format(final_step, added_suffix, final_step, self.output_dir, final_step))

                # Put the whole command together
                cmd = with_file + skip_text + out_text  # +override_text

                # Since we are getting the output from the final step directly
                # we can turn off the output from the pipeline
                cmd = cmd + ' --save_results=False'

                # For NIRCam we skip the odd even rows in refpix.
                if instrument == 'nircam':
                    cmd = cmd + ' --steps.refpix.odd_even_rows=False'
            cmds.append(cmd)
        return cmds

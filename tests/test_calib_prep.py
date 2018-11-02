'''Define unit tests for calib_prep with pytest.

Authors
-------
    - Bryan Hilbert

Use
---
    Ensure you have pytest installed. Then, simply run pytest in any
    parent directory of jwst_reffiles/tests/:
    >>> pytest
'''
from collections import OrderedDict
import copy
import os

from astropy.io import ascii
from jwst.pipeline import Detector1Pipeline

import jwst_reffiles
from jwst_reffiles.pipeline.calib_prep import CalibPrep
from jwst_reffiles.utils.definitions import PIPE_STEPS


def choose_file_test(cp_object):
    """Test function for choosing the proper filename to use
    as input to pipeline. This should be the file that has
    been run through the pipeline the farthest in a way that matches
    the requested final calibration state
    """
    files = ['file1.fits', 'file2.fits', 'file3.fits', 'file4.fits', 'file5.fits']

    # Create a step dictionary for each input file that lists its
    # current calibration state
    file1_state = create_step_dictionary('dq, sat')
    file2_state = create_step_dictionary('dq, sat, super')
    file3_state = create_step_dictionary('dq, sat, super, ref')
    file4_state = create_step_dictionary('dq, sat, super, ref, ipc')
    file5_state = create_step_dictionary('dq, sat, super, ref, ipc, jump')
    all_current = [file1_state, file2_state, file3_state, file4_state, file5_state]
    current_state = {}
    for i, filename in enumerate(files):
        current_state[filename] = all_current[i]

    # Create dictionary giving the required pipeline steps to compare
    required = []
    required.append(create_step_dictionary('dq, sat, super, ref, ipc, dark, jump'))
    required.append(create_step_dictionary('dq, sat'))
    required.append(create_step_dictionary('dq, sat, ref'))
    required.append(create_step_dictionary('dq, sat, super, ref, ipc, jump, rampfit'))
    # required.append(create_step_dictionary('dq'))

    # Select the best files
    chosen = []
    for req in required:
        chosen.append(cp_object.choose_file(files, req, current_state))

    truth = ['file4.fits', 'file1.fits', 'file1.fits', 'file5.fits']
    assert chosen == truth


def completed_steps_test(cp_object):
    """Test that header keywords describing the calibration state are
    correctly translated into a boolean dictionary"""
    filename = os.path.join(os.path.split(os.path.abspath(__file__))[0], 'test_data',
                            'file_for_completed_steps_test.fits')
    truth = create_step_dictionary('dq, sat, super, ref, lin')
    test_case = cp_object.completed_steps(filename)
    assert truth == test_case


def create_output_test(cp_object):
    """Test that the correct output filename is generated for a given
    requested calibration state and input file
    """
    requested = create_step_dictionary('dq, sat, ref, super, jump')
    basefilename = 'dark_current_file_number_42_uncal'
    base_and_suffix = basefilename.replace('uncal', 'dq_init_saturation_superbias_refpix_jump.fits')
    truth = os.path.join(cp_object.output_dir, base_and_suffix)
    constructed_name, true_basename = cp_object.create_output(basefilename, requested)
    assert constructed_name == truth


def create_step_dictionary(steplist):
    """Create a dictionary of pipeline steps with boolean entries
    given a list of steps to be set to True

    Parameters
    ----------
    steplist : list
        Comma separated list of pipeline steps to set to True

    Returns
    -------
    step_dict : dict
        Dictionary with requested steps set to True
    """
    all_false = copy.deepcopy(OrderedDict(PIPE_STEPS))
    for key in all_false:
        all_false[key] = False

    step_dict = copy.deepcopy(all_false)
    pipesteps = [element.strip() for element in steplist.split(',')]
    for stepname in pipesteps:
        if stepname in all_false.keys():
            step_dict[stepname] = True
        else:
            raise ValueError("{} is not a valid step dictionary entry.".format(stepname))
    return step_dict


def get_basename_test(cp_object):
    """Test of simple script to remove fits from filename"""
    filename = "jw00000999888_test_file_cal.fits"
    basename = cp_object.get_file_basename(filename)
    assert basename == filename.strip('.fits')


def initialize_calib_prep():
    """Create instance of calibPrep and prepare for later tests
    """
    instance = CalibPrep()
    instance.verbose = False
    # instance.inputs = ascii.read('test_calib_input_table.txt')
    instance.inputs = 'nothing'
    instance.search_dir = jwst_reffiles.__path__
    instance.output_dir = '/dummy/output/dir/'
    return instance


def steps_to_run_test(cp_object):
    """Test that the correct dictionary for calibration steps to be run
    is generated given an input file in some current state of calibration
    as well as a list of total calibration steps required.
    """
    required = create_step_dictionary('dq, sat, super, ref, dark, jump')
    current1 = create_step_dictionary('sat, ref, dark')
    current2 = create_step_dictionary('dq')
    truth1 = create_step_dictionary('dq, super, jump')
    truth2 = create_step_dictionary('sat, super, ref, dark, jump')
    torun1 = cp_object.steps_to_run('dummy.fits', required, current1)
    torun2 = cp_object.steps_to_run('dummy.fits', required, current2)
    assert torun1 == truth1
    assert torun2 == truth2


def strun_command_test(cp_object):
    """Test that the correct strun commands are generated for a given
    input file and required calibration state
    """
    from astropy.table import Column

    in_names = ['dummy.fits', 'dummy.fits', 'dummy.fits']
    input_names = Column(data=in_names)
    torun1 = 'sat,super,dark,jump'
    torun2 = 'dq,super,ref'
    torun3 = 'dark,lin,rampfit'
    steps_to_run = Column(data=[torun1, torun2, torun3])
    out_names = ['dummy_sat_superbias_dark_jump.fits', 'dummy_dq_init_superbias_refpix.fits',
                 'dummy_dark_linearity_rate.fits']
    suffixes = ['sat_superbias_dark_jump', 'dq_init_superbias_refpix', 'dark_linearity_rate']
    output_filename = Column(data=out_names)
    constructed_commands = cp_object.strun_command(input_names, steps_to_run, output_filename)
    truth1 = ('strun calwebb_detector1.cfg {} --steps.dq_init.skip=True --steps.refpix.skip=True '
              '--steps.ipc.skip=True --steps.linearity.skip=True --steps.persistence.skip=True '
              '--steps.ramp_fit.skip=True --steps.jump.suffix={} --steps.jump.output_dir={}'
              ' --steps.jump.save_results=True --save_results=False --steps.refpix.odd_even_rows=False'
              .format(in_names[0], suffixes[0], cp_object.output_dir))
    truth2 = ('strun calwebb_detector1.cfg {} --steps.saturation.skip=True --steps.ipc.skip=True '
              '--steps.linearity.skip=True --steps.persistence.skip=True --steps.dark_current.skip=True '
              '--steps.jump.skip=True --steps.ramp_fit.skip=True --steps.refpix.suffix={}'
              ' --steps.refpix.output_dir={} --steps.refpix.save_results=True --save_results=False '
              '--steps.refpix.odd_even_rows=False'
              .format(in_names[1], suffixes[1], cp_object.output_dir))
    truth3 = ('strun calwebb_detector1.cfg {} --steps.dq_init.skip=True --steps.saturation.skip=True '
              '--steps.superbias.skip=True --steps.refpix.skip=True --steps.ipc.skip=True '
              '--steps.persistence.skip=True --steps.jump.skip=True --steps.ramp_fit.suffix={}'
              ' --steps.ramp_fit.output_dir={} --steps.ramp_fit.save_results=True --save_results=False '
              '--steps.refpix.odd_even_rows=False'
              .format(in_names[2], suffixes[2], cp_object.output_dir))
    truths = [truth1, truth2, truth3]
    assert truths == constructed_commands


def test_calib_prep_steps():
    """Wrapper around tests of individual functions in calib.prep.py"""
    obj = initialize_calib_prep()
    choose_file_test(obj)
    get_basename_test(obj)
    choose_file_test(obj)
    completed_steps_test(obj)
    create_output_test(obj)
    steps_to_run_test(obj)
    strun_command_test(obj)


def test_pipeline_step_list():
    """Basic test that the number of steps in CALWEBB_DETECTOR1 is what we expect"""
    pipeline = Detector1Pipeline()
    assert len(pipeline.step_defs) == 15

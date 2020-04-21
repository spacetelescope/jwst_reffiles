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
from jwst_reffiles.pipeline.pipeline_steps import get_pipeline_steps
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
    file1_state = create_step_dictionary('dq_init, saturation')
    file2_state = create_step_dictionary('dq_init, saturation, superbias')
    file3_state = create_step_dictionary('dq_init, saturation, superbias, refpix')
    file4_state = create_step_dictionary('dq_init, saturation, superbias, refpix, ipc')
    file5_state = create_step_dictionary('dq_init, saturation, superbias, refpix, ipc, jump')
    all_current = [file1_state, file2_state, file3_state, file4_state, file5_state]
    current_state = {}
    for i, filename in enumerate(files):
        current_state[filename] = all_current[i]

    # Create dictionary giving the required pipeline steps to compare
    required = []
    required.append(create_step_dictionary('dq_init, saturation, superbias, refpix, ipc, dark_current, jump'))
    required.append(create_step_dictionary('dq_init, saturation'))
    required.append(create_step_dictionary('dq_init, saturation, refpix'))
    required.append(create_step_dictionary('dq_init, saturation, superbias, refpix, ipc, jump, rate'))
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
    truth = create_step_dictionary('dq_init, saturation, superbias, refpix, linearity')
    test_case = cp_object.completed_steps(filename)
    assert truth == test_case


def create_output_test(cp_object):
    """Test that the correct output filename is generated for a given
    requested calibration state and input file
    """
    requested = create_step_dictionary('dq_init, saturation, superbias, refpix, jump')
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
    steplist : str
        Comma separated list of pipeline steps to set to True

    Returns
    -------
    step_dict : dict
        Dictionary with requested steps set to True
    """
    all_false = OrderedDict()
    for key in get_pipeline_steps('nircam'):
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
    instance = CalibPrep('nircam')
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
    required = create_step_dictionary('dq_init, saturation, superbias, refpix, dark_current, jump')
    current1 = create_step_dictionary('saturation, refpix, dark_current')
    current2 = create_step_dictionary('dq_init')
    truth1 = create_step_dictionary('dq_init, superbias, jump')
    truth2 = create_step_dictionary('saturation, superbias, refpix, dark_current, jump')
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
    torun1 = 'saturation,superbias,dark_current,jump'
    torun2 = 'dq_init,superbias,refpix'
    torun3 = 'dark_current,linearity,rate'
    tosave1 = 'superbias,dark_current'
    tosave2 = 'dq_init'
    tosave3 = 'linearity'
    steps_to_run = Column(data=[torun1, torun2, torun3])
    steps_to_save = Column(data=[tosave1, tosave2, tosave3])
    out_names = ['dummy_sat_superbias_dark_jump', 'dummy_dq_init_superbias_refpix',
                 'dummy_dark_linearity_rate']
    suffixes = ['sat_superbias_dark_jump', 'dq_init_superbias_refpix', 'dark_linearity_rate']
    output_filename = Column(data=out_names)
    constructed_commands = cp_object.strun_command(input_names, steps_to_run, output_filename, steps_to_save, testing=True)
    pipeline_cfg_file = os.path.join(cp_object.output_dir, 'calwebb_detector1.cfg')
    truth1 = ('strun {} {} --steps.group_scale.skip=True --steps.dq_init.skip=True --steps.ipc.skip=True '
              '--steps.refpix.skip=True --steps.linearity.skip=True --steps.persistence.skip=True '
              '--steps.ramp_fit.skip=True --steps.jump.output_file={} --steps.jump.save_results=True '
              '--output_dir={} --steps.superbias.save_results=True --steps.dark_current.save_results=True '
              '--steps.superbias.output_file=/dummy/output/dir/dummy_saturation_superbias.fits '
              '--steps.dark_current.output_file=/dummy/output/dir/dummy_saturation_superbias_dark_current.fits '
              '--save_results=False --steps.jump.rejection_threshold=100 --steps.refpix.odd_even_rows=False'
              .format(pipeline_cfg_file, in_names[0], out_names[0], cp_object.output_dir))
    truth2 = ('strun {} {} --steps.group_scale.skip=True --steps.saturation.skip=True --steps.ipc.skip=True '
              '--steps.linearity.skip=True --steps.persistence.skip=True --steps.dark_current.skip=True '
              '--steps.jump.skip=True --steps.ramp_fit.skip=True --steps.refpix.output_file={} '
              '--steps.refpix.save_results=True --output_dir={} --steps.dq_init.save_results=True '
              '--steps.dq_init.output_file=/dummy/output/dir/dummy_dq_init.fits --save_results=False '
              '--steps.jump.rejection_threshold=100 --steps.refpix.odd_even_rows=False'
              .format(pipeline_cfg_file, in_names[1], out_names[1], cp_object.output_dir))
    truth3 = ('strun {} {} --steps.group_scale.skip=True --steps.dq_init.skip=True --steps.saturation.skip=True '
              '--steps.ipc.skip=True --steps.superbias.skip=True --steps.refpix.skip=True '
              '--steps.persistence.skip=True --steps.jump.skip=True --steps.ramp_fit.output_file={} '
              '--steps.ramp_fit.save_results=True --output_dir={} --steps.linearity.save_results=True '
              '--steps.linearity.output_file=/dummy/output/dir/dummy_linearity_dark_current.fits --save_results=False '
              '--steps.jump.rejection_threshold=100 --steps.refpix.odd_even_rows=False'
              .format(pipeline_cfg_file, in_names[2], out_names[2], cp_object.output_dir))
    truths = [truth1, truth2, truth3]


    print('')
    print(constructed_commands[0])
    print(constructed_commands[1])
    print(constructed_commands[2])
    print('')
    print(truth1)
    print(truth2)
    print(truth3)
    print('')


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

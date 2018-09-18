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

from jwst.pipeline import Detector1Pipeline

from jwst_reffiles.pipeline.calib_prep import calibPrep
from jwst_reffiles.utils.definitions import PIPE_STEPS


def initialize_calib_prep():
    """Create instance of calibPrep and prepare for later tests
    """
    instance = calibPrep()
    instance.verbose = False
    instance.inputs =
    instance.search_dir = ''
    instance.output_dir = ''
    return instance


def test_pipeline_step_list():
    """Basic test that the number of steps in CALWEBB_DETECTOR1 is what we expect"""
    pipeline = Detector1Pipeline()
    assert len(pipeline.step_defs) == 15


def test_calib_prep_steps():
    """Wrapper around tests of individual functions in calib.prep.py"""
    obj = initialize_calib_prep()
    choose_file_test(obj)
    get_basename_test(obj)


def choose_file_test(cp_object):
    """Test function for choosing the proper filename to use
    as input to pipeline. This should be the file that has
    been run through the pipeline the farthest in a way that matches
    the requested final calibration state
    """
    files = ['file1.fits', 'file2.fits', 'file3.fits', 'file4.fits', 'file5.fits']

    # Create a basic copy of the pipeline step dictionary
    all_false = copy.deepcopy(PIPE_STEPS)
    for key in all_false:
        all_false[key] = False

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
    required.append(create_step_dictionary('dq, sat, super, ref, ipc, jump, rate'))

    # Select the best files
    chosen = []
    for req in required:
        chosen.append(cp_object.choose_file(files, req, current_state))

    truth = ['file4.fits', 'file1.fits', 'file1.fits', 'file5.fits']
    assert chosen == truth


def get_basename_test(cp_object):
    """Test of simple script to remove fits from filename"""
    filename = "jw00000999888_test_file_cal.fits"

    basename = cp_object.get_file_basename(filename)
    assert basename == filename.strip('.fits')


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
    all_false = copy.deepcopy(PIPE_STEPS)
    for key in all_false:
        all_false[key] = False

    step_dict = {}
    pipesteps = [PIPE_STEPS[element.strip()] for element in steplist.split(',')]
    for stepname in pipesteps:
        step_dict[stepname] = True
    return step_dict


def test_trivial():
    assert 1 == 1

#! /usr/bin/env python

"""This module contains code for retrieving a list of JWST calibration pipeline steps
to help `jwst_reffiles` define the names and order of the default pipeline steps, such that
we can make shortcuts (e.g. "slope-" = run all steps in the nominal order up to ramp-fitting)
"""

# from jwst.pipeline import Detector1Pipeline
from jwst_reffiles.utils.definitions import INSTRUMENTS


def get_pipeline_steps(instrument):
    """Get the names and order of the calwebb_detector1
    pipeline steps for a given instrument. Use values that match up with the values in the
    PIPE_STEP defintion in definitions.py

    Parameters
    ----------

    instrument : str
        Name of JWST instrument

    Returns
    -------

    steps : dict
        Dictionary of step names (and modules? do we care?)
    """
    instrument = instrument.lower()
    if instrument not in INSTRUMENTS:
        raise ValueError("WARNING: {} is not a valid instrument name.".format(instrument))
    # all_steps = Detector1Pipeline.step_defs

    # Order is important in 'steps' lists below!!
    if instrument == 'miri':
        steps = ['group_scale', 'dq_init', 'saturation', 'ipc', 'firstframe', 'lastframe',
                 'linearity', 'rscd', 'dark_current', 'refpix', 'persistence', 'jump', 'rate']
        # No persistence correction for MIRI
        steps.remove('persistence')
    else:
        steps = ['group_scale', 'dq_init', 'saturation', 'ipc', 'superbias', 'refpix', 'linearity',
                 'persistence', 'dark_current', 'jump', 'rate']

        # No persistence correction for NIRSpec
        if instrument == 'nirspec':
            steps.remove('persistence')
    return steps


def step_minus(step_name, instrument):
    """Enable "-" shorthand when talking about pipeline steps. For example,
    "rate-" will indicate all pipeline steps up to (but not including) ramp-fitting.
    For a given step_name-, return the list of step names to be run.

    Parameters
    ----------

    step_name : str
        Name of step name (including the trailing '-')

    instrument : str
        Name of JWST instrument

    Returns
    -------

    step_list : list
        List of step names up to but not including the input step
    """
    step_name = step_name.split('-')[0].lower()
    all_steps = get_pipeline_steps(instrument)

    # Basic error check
    if step_name not in all_steps:
        raise ValueError("WARNING: {} is not a valid pipeline step name.".format(step_name))

    step_list = []
    for step in all_steps:
        if step == step_name:
            break
        else:
            step_list.append(step)

    # Assuming this will be called by mkrefs, we need to return a comma-separated list of steps
    step_list_csv = str(step_list).strip('[]').replace("'", "")
    return step_list_csv


def step_plus(step_name, instrument):
    """Enable "+" shorthand when talking about pipeline steps. For example,
    "rate+" will indicate all pipeline steps up to and including ramp-fitting.
    For a given step_name+, return the list of step names to be run.

    Parameters
    ----------

    step_name : str
        Name of step name (including the trailing '+')

    instrument : str
        Name of JWST instrument

    Returns
    -------

    step_list : list
        List of step names up to and including the input step
    """
    step_name = step_name.split('+')[0].lower()
    all_steps = get_pipeline_steps(instrument)

    # Basic error check
    if step_name not in all_steps:
        raise ValueError("WARNING: {} is not a valid pipeline step name.".format(step_name))

    location = all_steps.index(step_name)
    step_list = all_steps[0:location + 1]

    # Assuming this will be called by mkrefs, we need to return a comma-separated list of steps
    step_list_csv = str(step_list).strip('[]').replace("'", "")
    return step_list_csv

#! /usr/bin/env python

"""File containing various necessary definitions, such as abbreviations for each
pipeline step
"""

# Abbreviations for pipeline steps. The abbreviations (keys) here, must be used
# in the configuration files when listing pipeline steps
#PIPE_STEPS = [('dq', 'dq_init'), ('sat', 'saturation'), ('super', 'superbias'),
#              ('ref', 'refpix'), ('ipc', 'ipc'), ('lin', 'linearity'),
#              ('persistence', 'persistence'), ('dark', 'dark_current'),
#              ('jump', 'jump'), ('rampfit', 'rate')]

# Abbreviations for pipeline steps. These are the values used in the JWST calibration
# pipeline, and are also the values that must be used in jwst_reffiles input configuration
# files
PIPE_STEPS = ['group_scale', 'dq_init', 'saturation', 'ipc', 'firstframe', 'lastframe',
              'superbias', 'refpix', 'linearity', 'persistence', 'rscd', 'dark_current',
              'jump', 'rate']

# Define header keyword that accompanies each step
# This needs to be made more instrument-agnostic. The list is going to
# vary depending on instrument.
#PIPE_KEYWORDS = {'S_DQINIT': 'dq', 'S_SATURA': 'sat', 'S_REFPIX': 'ref', 'S_SUPERB': 'super',
#                 'S_IPC': 'ipc', 'S_PERSIS': 'persistence', 'S_DARK': 'dark', 'S_LINEAR': 'lin',
#                 'S_JUMP': 'jump',  'S_RAMP': 'rampfit'}

PIPE_KEYWORDS = {'S_DQINIT': 'dq_init', 'S_SATURA': 'saturation', 'S_REFPIX': 'refpix', 'S_SUPERB': 'superbias',
                 'S_IPC': 'ipc', 'S_PERSIS': 'persistence', 'S_DARK': 'dark_current', 'S_LINEAR': 'linearlity',
                 'S_JUMP': 'jump',  'S_RAMP': 'rate'}

INSTRUMENTS = ['nircam', 'niriss', 'nirspec', 'miri', 'fgs']

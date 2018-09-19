#! /usr/bin/env python

"""File containing various necessary definitions, such as abbreviations for each
pipeline step
"""

# Abbreviations for pipeline steps. The abbreviations (keys) here, must be used
# in the configuration files when listing pipeline steps
PIPE_STEPS = [('dq', 'dq_init'), ('sat', 'saturation'), ('super', 'superbias'),
              ('ref', 'refpix'), ('ipc', 'ipc'), ('lin', 'linearity'),
              ('persistence', 'persistence'), ('dark', 'dark_current'),
              ('jump', 'jump'), ('rampfit', 'rate')]

# Define header keyword that accompanies each step
# This needs to be made more instrument-agnostic. The list is going to
# vary depending on the instrument.
PIPE_KEYWORDS = {'S_DQINIT': 'dq', 'S_SATURA': 'sat', 'S_REFPIX': 'ref', 'S_SUPERB': 'super',
                 'S_IPC': 'ipc', 'S_PERSIS': 'persistence', 'S_DARK': 'dark', 'S_LINEAR': 'lin',
                 'S_JUMP': 'jump',  'S_RAMP': 'rampfit'}

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


def test_pipeline_step_list():
    """Basic test that the number of steps in CALWEBB_DETECTOR1 is what we expect"""
    pipeline = Detector1Pipeline()
    assert len(pipeline.step_defs) == 15


def test_trivial():
    assert 1 == 1

#!/usr/bin/env python
'''
example script that is loaded by mkref_example_bpm.py
A. Rest
'''

import numpy as np
from jwst.datamodels import MaskModel

# convenience function
def make_maskmodel(ysize, xsize):

    # create a mask model for the dq_init step
    dq = np.zeros((ysize, xsize), dtype=int)

    # define a dq_def extension
    mask = MaskModel()

    dqdef = [(0, 1, 'DO_NOT_USE', 'Bad Pixel do not use'),
             (1, 2, 'DEAD', 'Dead Pixel'),
             (2, 4, 'HOT', 'Hot pixel'),
             (3, 8, 'UNRELIABLE_SLOPE', 'Large slope variance'),
             (4, 16, 'RC', 'RC pixel'),
             (5, 32, 'REFERENCE_PIXEL', 'Reference Pixel')]

    dq_def = np.array((dqdef), dtype=mask.dq_def.dtype)

    return dq, dq_def

class example_mkbpmclass:
    
    def mkbpm(self,dark1filename,dark2filename,flat1filename,flat2filename,
              whateveroption1=False,whateveroption2=2.0):

        print('HELLO mkbpm in myexamplescript.py!!!')

        # create an empty MaskModel for the bad pixel mask
        xstart = 1
        ystart = 1
        xsize = 2048
        ysize = 2048
        
        dq, dq_def = make_maskmodel(ysize, xsize)
        
        # write mask model and fix the meta data
        mask_model = MaskModel(dq=dq, dq_def=dq_def)

        return(mask_model)
        

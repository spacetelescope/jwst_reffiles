#! /usr/bin/env python

'''
Create gain reference file using SSB datamodels
given the final gain array and error maps
'''

import os
from jwst.datamodels import GainModel

class GainFile:
    def __init__(self):
        self.detector = ''
        self.outfile = ''
        self.author = ''
        self.descrip = ''
        self.useafter = ''
        self.pedigree = ''
        self.history = ''
        self.outdir = '.'
        self.fastaxis = 0
        self.slowaxis = 0
        
    def save(self,array,error):
        # Save using datamodel
        mod = GainModel()
        mod.data = array
        #mod.err = error
        mod.meta.telescope = 'JWST'
        mod.meta.author = self.author
        mod.meta.description = self.descrip
        mod.meta.useafter = self.useafter
        mod.meta.instrument.name = 'NIRCAM'
        mod.meta.instrument.detector = self.detector
        mod.meta.pedigree = self.pedigree
        mod.meta.reftype = 'GAIN'
        mod.meta.subarray.fastaxis = self.fastaxis
        mod.meta.subarray.slowaxis = self.slowaxis
        mod.meta.subarray.name = 'GENERIC'
        mod.meta.subarray.xsize = 2048
        mod.meta.subarray.ysize = 2048
        mod.meta.subarray.xstart = 1
        mod.meta.subarray.ystart = 1
        #mod = self.add_history(mod,self.history)
        mod.history.append(self.history)
        mod.save(os.path.join(self.outdir,self.outfile))

    def add_history(model,hist):
        # Break history entry into a list of small
        # strings and add individually
        pass

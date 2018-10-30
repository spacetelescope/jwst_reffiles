#!/usr/bin/env python
'''
create read noise
A. Rest
'''

import argparse,os,re,sys,types

from mkref_template import mkrefclass

from jwst_reffiles.utils.tools import astrotableclass,yamlcfgclass

class mkrefclassX(mkrefclass):
    def __init__(self,*args, **kwargs):
        mkrefclass.__init__(self,*args, **kwargs)
        self.reflabel='bpm'
        self.reftype='bpm'

    def inputfileoptions(self, parser):
        parser.add_argument('--D1',
                            help='specify the dark exposure (default=%(default)s)')
        return(0)
        
    def extraoptions(self, parser):
        parser.add_argument('-x','--dummy_bpm_option1',
                            help='testoption1 bpm')
        return(0)


if __name__ == '__main__':
    print('bbb',sys.argv)
    mkref = mkrefclassX()
    (parser) = mkref.refoptions()
    args = parser.parse_args()
    print('bbb',args)
    sys.exit(0)
    #mkref.mkref(sys.argv[1:])

#!/usr/bin/env python
'''
create gain with gain_armin method
A. Rest
'''

import argparse,os,re,sys,types

from mkref_template import mkrefclass

from jwst_reffiles.utils.tools import astrotableclass,yamlcfgclass

class mkrefclassX(mkrefclass):
    def __init__(self,*args, **kwargs):
        mkrefclass.__init__(self,*args, **kwargs)
        self.reflabel='gain_armin'
        self.reftype='gain'
        
    def inputfileoptions(self, parser):
        parser.add_argument('--D1',required=True,
                            help='specify the first dark exposure (default=%(default)s)')
        parser.add_argument('--D2',required=True,
                            help='specify the second dark exposure (default=%(default)s)')
        parser.add_argument('--F1',required=True,
                            help='specify the first flat exposure (default=%(default)s)')
        parser.add_argument('--F2',required=True,
                            help='specify the second flat exposure (default=%(default)s)')
        return(0)
        
    def extraoptions(self, parser):
        parser.add_argument('-d', '--dummy_gain_option1', default=None, nargs=2,
                            help='testoption1 gainarmin')
        parser.add_argument('--dummy_gain_option2', default=None,
                            help='testoption2 gainarmin')
        return(0)

    def callagorithm(self,args):
        return(0)
        #import myscript from bla
        #myscript.myroutine(args.D1,args.F1,newgain=args.dummy_gain_option1)
    
if __name__ == '__main__':
    mkref = mkrefclassX()
    (parser) = mkref.refoptions()
    args = parser.parse_args()
    mkref.callagorithm(args)

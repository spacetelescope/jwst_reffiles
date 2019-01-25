#!/usr/bin/env python
'''
create gain with gain_armin method
A. Rest
'''

import argparse,os,re,sys,types

from plugin_wrapper import mkrefclass_template

from jwst_reffiles.utils.tools import astrotableclass,yamlcfgclass

class mkrefclass(mkrefclass_template):
    def __init__(self,*args, **kwargs):
        mkrefclass_template.__init__(self,*args, **kwargs)
        self.reflabel='gain_armin'
        self.reftype='gain'

    def extra_optional_arguments(self, parser):
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
    mkref = mkrefclass()
    parser = mkref.allargs()
    args = parser.parse_args()
    mkref.callagorithm(args)

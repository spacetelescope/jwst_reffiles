#!/usr/bin/env python
'''
create read noise
A. Rest
'''

import argparse,os,re,sys,types

from plugin_wrapper import mkrefclass_template

from jwst_reffiles.utils.tools import astrotableclass,yamlcfgclass

class mkrefclass(mkrefclass_template):
    def __init__(self,*args, **kwargs):
        mkrefclass_template.__init__(self,*args, **kwargs)
        self.reflabel='bpm'
        self.reftype='bpm'

    def extra_optional_arguments(self, parser):
        parser.add_argument('-x','--dummy_bpm_option1',nargs=2,
                            help='testoption1 bpm')
        return(0)


if __name__ == '__main__':
    mkref = mkrefclass()
    parser = mkref.allargs()
    args = parser.parse_args()
    mkref.callagorithm(args)

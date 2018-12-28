#!/usr/bin/env python
'''
create read noise
A. Rest
'''

import argparse,os,re,sys,types

from mkref_template import mkrefclass_template

from jwst_reffiles.utils.tools import astrotableclass,yamlcfgclass

class mkrefclass(mkrefclass_template):
    def __init__(self,*args, **kwargs):
        mkrefclass_template.__init__(self,*args, **kwargs)
        self.reflabel='rdnoisenircam'
        self.reftype='rdnoise'

    def extra_optional_arguments(self, parser):
        parser.add_argument('-x','--dummy_rdnoise_option1',
                            help='testoption1 readnoise')
        parser.add_argument('--dummy_rdnoise_option2',
                            help='testoption2 readnoise')
        return(0)


if __name__ == '__main__':
    mkref = mkrefclass()
    parser = mkref.allargs()
    args = parser.parse_args()
    mkref.callagorithm(args)

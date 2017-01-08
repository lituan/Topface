#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage: python get_hotspot.py *.wdsp
output hotspots in following format
pro xxx xxx xxx xxx xxx xxx
"""

import os
import sys
import lt
from wdsp import Wdsp

with open(sys.argv[-1]) as wdsp_f:
    w = Wdsp(wdsp_f)
    with lt.open_file(file_suffix='hotspots') as w_f:
        for pro, hots in w.hotspots.iteritems():
            print >> w_f, '{0:<25}{1:<}'.format(pro,' '.join(hots))


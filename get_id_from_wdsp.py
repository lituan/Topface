#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
output protein id from a wdsp file
"""
import os
import sys

with open(sys.argv[-1]) as o_f:
    lines = o_f.readlines()
    ids = [line.split()[1] for line in lines if '>' in line]

with open('wdsp_id.txt', 'w') as w_f:
    for i in ids:
        print >> w_f, i

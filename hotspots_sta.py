#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
each blade has three hotspots, in theory, there will be 20**3 possible combinations
so, in reality, is there any prefer for some combination
"""

import itertools
import os
import sys

from wdsp import Wdsp


def classify_blade(wdsp_f):
    with open(wdsp_f) as o_f:
        wdsp = Wdsp(o_f)
        hotspots = ' '.join([' '.join(v)
                             for k, v in wdsp.hotspots.iteritems()]).split()

        aa = {'K', 'R', 'H', 'D', 'E', 'F', 'W', 'Y', 'S', 'T',
              'N', 'Q', 'V', 'L', 'I', 'M', 'A', 'C', 'P', 'G', '*', 'X', 'B'}
        aa_combi = itertools.product(aa, repeat=3)
        blades = {}
        for c in aa_combi:
            blades[''.join(c)] = 0

        for hot in hotspots:
            blades[hot] += 1

        blades = [(v, k) for k, v in blades.iteritems()]
        blades = sorted(blades, reverse=True)

        return blades


def main():
    wdsp_f = sya.argv[-1]
    blades = classify_blade(wdsp_f)

    with open('hotspots_sta.txt', 'w') as w_f:
        for v, k in blades:
            print >> w_f, '{0:<10}{1:<}'.format(k, v)

if __name__ == "__main__":
    main()

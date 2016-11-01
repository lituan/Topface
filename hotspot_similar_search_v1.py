#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import itertools
import operator
from collections import defaultdict

import lt
from wdsp import Wdsp

"""
find similar hotspots
similarity cutoff is set by using blades identity
b values is 1 2 3, which means at least b residues are exactly the same
h valuse is [1,2,3,4,5,6,7,8,9] means at least h blases are similar
"""


def check_hot_similar(hot1, hot2, b, h):
    if len(hot1) < h or len(hot2) < h:
        return 0
    hot1h = [(hot1 + hot1)[i:i + h] for i in xrange(len(hot1))]
    hot2h = [(hot2 + hot2)[i:i + h] for i in xrange(len(hot2))]

    for h1 in hot1h:
        for h2 in hot2h:
            success = 1
            for h11, h22 in zip(h1, h2):
                common = [1 for r1, r2 in zip(h11, h22) if r1 == r2]
                if len(common) < b:
                    success = 0
                    break
            if success:
                return 1
    return 0


def get_similar_hots(tem_hots, all_hots, b, h):
    tem_all_hots = {}
    for pro, hot in tem_hots.iteritems():
        tem_all_hots[pro] = []
        for p, h in all_hots.iteritems():
            if check_hot_similar(hot, h, b, h):
                tem_all_hots[pro].append(p)

    return tem_all_hots


def write_result(tem_all_hots, tem_hots, all_hots, b, h, cutoff=1):
    for tem_pro, all_pros in tem_all_hots.iteritems():
        f_name = tem_pro + '_similar_hotspots_b' + str(b) + 'h' + str(h)
        with lt.open_file(f_name) as w_f:
            print >> w_f, '{0:<20}{1:<}'.format(
                tem_pro, ' '.join(tem_hots[tem_pro]))
            for pro in all_pros:
                print >> w_f, '{0:<20}{1:<}'.format(
                    pro, ' '.join(all_hots[pro]))


def read_hots(wdsp_f):
    with open(sys.argv[-2]) as wdsp_f:
        wdsp = Wdsp(wdsp_f)
        hots = wdsp.hotspots

        return hots


@lt.run_time
def main():
    tem_hots = read_hots(sys.argv[-2])
    all_hots = read_hots(sys.argv[-1])
    for b in 2, 3:
        for h in 2, 3, 4, 5, 6, 7:
            tem_all_hots = get_similar_hots(tem_hots, all_hots, b, h)
            write_result(tem_all_hots, tem_hots, all_hots, b, h)
    for b in [1]:
        for h in 3, 4, 5, 6, 7:
            tem_all_hots = get_similar_hots(tem_hots, all_hots, b, h)
            write_result(tem_all_hots, tem_hots, all_hots, b, h)


if __name__ == "__main__":
    main()

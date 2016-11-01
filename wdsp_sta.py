#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
calculate statistics for WDSP output file
usage: python wdsp_sta.py wdsp_f
"""
import lt
import sys
import os

from wdsp import Wdsp

with open(sys.argv[-1]) as o_f:
    wdsp = Wdsp(o_f)

    scores_sta = lt.lis_sta(wdsp.scores.values())
    with lt.open_file(file_suffix='total_score_sta') as w_f:
        for num, freq in scores_sta:
            print >> w_f, '{0:<10}{1}'.format(num, freq)

    tetrad_sta = [len([vi for vi in v if vi >= 44.0])
                  for k, v in wdsp.blade_scores.iteritems()]
    tetrad_sta = lt.lis_sta(tetrad_sta)
    with lt.open_file(file_suffix='tetrad_num_sta') as w_f:
        for num, freq in tetrad_sta:
            print >> w_f, '{0:<5}{1}'.format(num, freq)

    blades_sta = [len(blades) for pro, blades in wdsp.blades.iteritems()]
    blades_sta = lt.lis_sta(blades_sta)
    with lt.open_file(file_suffix='blades_sta') as w_f:
        for num, freq in blades_sta:
            print >> w_f, '{0:<5}{1}'.format(num, freq)

    R1 = []
    R1_2 = []
    D_1 = []
    for pro, hotspot in wdsp.hotspots.iteritems():
        for hot in hotspot:
            R1.append(hot[0])
            R1_2.append(hot[1])
            D_1.append(hot[2])
    R1_sta = lt.lis_sta(R1)
    R1_2_sta = lt.lis_sta(R1_2)
    D_1_sta = lt.lis_sta(D_1)
    with lt.open_file(file_suffix='R1_AA_sta') as w_f:
        for num, freq in R1_sta:
            print >> w_f, '{0:<5}{1}'.format(num, freq)
    with lt.open_file(file_suffix='R1_2_AA_sta') as w_f:
        for num, freq in R1_2_sta:
            print >> w_f, '{0:<5}{1}'.format(num, freq)
    with lt.open_file(file_suffix='D_1_AA_sta') as w_f:
        for num, freq in D_1_sta:
            print >> w_f, '{0:<5}{1}'.format(num, freq)

    # repeat length distribution between proteins
    repeat_sta = [len(''.join(blade))
                  for pro, blades in wdsp.blades.iteritems() for blade in blades]
    with lt.open_file(file_suffix='repeat_num') as w_f:
        for i in repeat_sta:
            print >> w_f, i
    repeat_sta = lt.lis_sta(repeat_sta)
    with lt.open_file(file_suffix='repeat_sta') as w_f:
        for num, freq in repeat_sta:
            print >> w_f, '{0:<5}{1}'.format(num, freq)
    # repeat length distribution within proteins
    repeat = []
    for pro, blades in wdsp.blades.iteritems():
        pro_repeat = []
        for blade in blades:
            pro_repeat.append(len(''.join(blade)))
        repeat.append(pro_repeat)
    # protein sequence length
    pro_len = [len(seq) for pro, seq in wdsp.seqs.iteritems()]
    with lt.open_file(file_suffix='pro_len') as w_f:
        for i in pro_len:
            print >> w_f, i
    # blade length
    blade_sta = [[len(b) for b in blade]
                 for pro, blades in wdsp.blades.iteritems() for blade in blades]
    d = [i[0] for i in blade_sta]
    da = [i[1] for i in blade_sta]
    a = [i[2] for i in blade_sta]
    ab = [i[3] for i in blade_sta]
    b = [i[4] for i in blade_sta]
    bc = [i[5] for i in blade_sta if len(i) >= 6]
    c = [i[6] for i in blade_sta if len(i) >= 7]
    cd = [i[7] for i in blade_sta if len(i) >= 8]

    print b

    with lt.open_file(file_suffix='d') as w_f:
        for i in d:
            print >> w_f, i
    with lt.open_file(file_suffix='da') as w_f:
        for i in da:
            print >> w_f, i
    with lt.open_file(file_suffix='a') as w_f:
        for i in a:
            print >> w_f, i
    with lt.open_file(file_suffix='ab') as w_f:
        for i in ab:
            print >> w_f, i
    with lt.open_file(file_suffix='b') as w_f:
        for i in b:
            print >> w_f, i
    with lt.open_file(file_suffix='bc') as w_f:
        for i in bc:
            print >> w_f, i
    with lt.open_file(file_suffix='c') as w_f:
        for i in c:
            print >> w_f, i
    with lt.open_file(file_suffix='cd') as w_f:
        for i in cd:
            print >> w_f, i

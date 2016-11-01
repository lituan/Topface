#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
calculate repeats similaritites for WDSP output file
usage: python repeat_similarity.py wdsp_f
"""
import sys
from collections import OrderedDict
from wdsp import Wdsp


def align(seq1, seq2):
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62
    gap_open = -10  # usual value
    gap_extend = -0.5  # usual value

    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)

    seq1 = alns[0][0]
    seq2 = alns[0][1]
    identity = [1 for i, s in enumerate(seq1) if s == seq2[i]]
    identity = int(100 * len(identity) / len(seq1))

    return identity


def repeat_similarity(repeats):
    lens = len(repeats)
    sims = []
    for i in xrange(lens):
        sim_i = []
        for j in xrange(lens):
            if j < i:
                sim_i.append(sims[j][i])
            elif j >= i:
                sim_i.append(align(repeats[i], repeats[j]))
        sims.append(sim_i)
    average = (sum([sum(i) for i in sims]) - lens * 100) / (lens * (lens - 1))
    return average, sims


def main():
    with open(sys.argv[-1]) as o_f:
        w = Wdsp(o_f)
        sims = OrderedDict()
        for pro, repeats in w.repeats.iteritems():
            sims[pro] = repeat_similarity(repeats)

        with open('sims.txt', 'w') as w_f:
            for k, v in sims.iteritems():
                print >> w_f, '{0:<20}{1:<}'.format(k, v[0])

if __name__ == "__main__":
    main()

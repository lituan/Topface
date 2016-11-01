#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
calculate sequence similarities for two WDSP output file
usage: python seq_similarity.py wdsp_f1 wdsp_f2
"""
import sys
import operator
import lt
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


@lt.run_time
def main():
    with open(sys.argv[-2]) as o_f:
        tem = Wdsp(o_f)
        tem_seq = tem.seqs
    with open(sys.argv[-1]) as o_f:
        all1 = Wdsp(o_f)
        all_seq = all1.seqs

    similarity = OrderedDict()
    for t_name, t_seq in tem_seq.iteritems():
        sim = []
        for a_name, a_seq in all_seq.iteritems():
            sim.append((a_name, align(t_seq, a_seq)))
        # sim = sorted(sim, key=operator.itemgetter(1),reverse=True)
        similarity[t_name] = sim

    for k, v in similarity.iteritems():
        with lt.open_file(k) as w_f:
            for a_name, a_identity in v:
                print >> w_f, '{0:<15}{1:<}'.format(a_name, a_identity)


if __name__ == "__main__":
    main()

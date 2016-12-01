#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
reform sequence so that each line contain 80 residues
and simplize uniprot-type id
"""

import os
import sys
from cPickle import load


def readfa(fa_f):
    # readin seqs in fasta format
    # seqs foramt:[(pro,seq),...]
    with open(fa_f,'r') as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        pro_line_num = [i for i, line in enumerate(
            lines) if '>' in line] + [len(lines)]
        seqs = [lines[n:pro_line_num[i + 1]]
                for i, n in enumerate(pro_line_num[:-1])]
        seqs = [(seq[0][1:], ''.join(seq[1:])) for seq in seqs]
        return seqs

def map_to_uniprot(oma_id):
    with open('oma-uniprot.txt') as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n').split() for line in lines if line and (not '#' in line)]
        map_dic = {}
        for k,v in lines:
            if not map_dic.has_key(k):
                map_dic[k] = [v]
            else:
                map_dic[k].append(v)
        lt.pickle_dump(map_dic,'oma_map_to_uniprot')
        return map_dic[oma_id]

def map_to_uniprot(oma_id):
    dic = load(open('oma_map_to_uniprot.pickle','r'))
    if dic.has_key(oma_id):
        return dic[oma_id]
    else:
        return False

def writefa(seqs,filename):
    with open(filename+'_map_to_uniprot.fa', 'w') as w_f:
        seqs = [(pro.split()[2], seq) for pro, seq in seqs] # oma-fasta
        # seqs = [(pro, seq) for pro, seq in seqs]
        # seqs = sorted(seqs)
        dic = load(open('oma_map_to_uniprot.pickle','r'))
        for pro, seq in seqs:
            if dic.has_key(pro):
                pro = '-'.join(dic[pro])
                print >> w_f, '>{}'.format(pro)
            else:
                print >> w_f, '>{}'.format(pro)
            for i in [seq[i:i + 80] for i in range(0, len(seq), 80)]:
                print >> w_f, i

def main():
    fa_f = sys.argv[-1]
    filename = os.path.splitext(os.path.split(fa_f)[1])[0]
    seqs = readfa(fa_f)
    writefa(seqs,filename)

main()

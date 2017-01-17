#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
transform fasta file downloaded from uniprot to wdsp input format
wdsp input format:
each seq has two lines, one for title (there is space between > and pro_name), one for sequences
> H0X5L8
GIDEPLHIKRRKVIKPGFIHSPWK...

usage: python make_wdsp_input,py 50 sequences.fasta
"""
import sys
import os

def readfa(fa_f):
    with open(fa_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        pro_line_num = [i for i,line in enumerate(lines) if '>' in line] + [len(lines)]
        print 'sequences num: ',len(pro_line_num)-1
        seqs = [lines[n:pro_line_num[i+1]] for i,n in enumerate(pro_line_num[:-1])]
        seqs = [(seq[0][1:].split('|')[1],''.join(seq[1:])) for seq in seqs]
        return seqs

def writefa(seqs,seg_len):
    f_path,f_name = os.path.split(sys.argv[-1])
    f_name,f_ext = os.path.splitext(f_name)
    seqgroups = [seqs[i:i+seg_len] for i in range(0,len(seqs),seg_len)]
    print 'sequences groups num: ',len(seqgroups)
    for i,seqs in enumerate(seqgroups):
        with open(f_name+'_wdsp_input_'+str(i)+'.fa','w') as w_f:
            for pro,seq in seqs:
                print >> w_f,'> {0:<}'.format(pro)
                print >> w_f,seq

def main():
    seqs = readfa(sys.argv[-1])
    seg_len = int(sys.argv[-2])
    writefa(seqs,seg_len)

if __name__ == "__main__":
    main()

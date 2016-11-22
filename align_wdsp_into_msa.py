#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
input wdsp out file, output msa by aligning  by structures
"""
import os
import sys

def trans_lis_lis(lis_lis):
    """align and trans nested list to print a table"""
    lis_lis = [[str(l) for l in lis]
               for lis in lis_lis]  # trans every element to str
    lis_lis_len = len(lis_lis)
    lis_max_len = max(len(lis) for lis in lis_lis)
    # construct a matrix
    lis_lis = [lis + (lis_max_len - len(lis)) * [''] for lis in lis_lis]
    # transform the matrix to its T matrix
    trans = [[lis_lis[i][j]
              for i in range(lis_lis_len)] for j in range(lis_max_len)]
    # aligned = [[l + (max([len(l) for l in lis]) - len(l)) * " " for l in lis] for lis in trans]
    aligned = []
    for lis in trans:
        width = max([len(l) for l in lis])
        lis = [l + (width - len(l)) * '-' for l in lis]
        aligned.append(lis)
    # transform the T matrix to the original  matrix
    trans = [[aligned[i][j]
              for i in range(lis_max_len)] for j in range(lis_lis_len)]
    return trans

def main():
    with open(sys.argv[-1]) as o_f:
        lines = o_f.readlines()
        # strip end of line symbol and empty lines
        lines = [line.strip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        # cut into entries
        begin = [i for i, line in enumerate(lines) if '>' in line]
        end = begin[1:] + [len(lines)]
        entries = [lines[begin[i]:end[i]] for i in range(len(begin))]
        entries = [[l.split() for l in entry] for entry in entries]
        # check last line of each entry
        for entry in entries:
            if len(entry[-1]) < len(entry[1]):
                complement = len(entry[1]) - len(entry[-1])
                entry[-1] = entry[-1][0:-1] + \
                    complement * [' '] + [entry[-1][-1]]
        # align each entry
        entries_new = []
        for entry in entries:
            combine = ['{0:<40}'.format(entry[0][1])]
            for b in [blade[3:-1] for blade in entry[1:]]:
                combine +=  b
            entries_new.append(combine)
        entries_new = trans_lis_lis(entries_new)

        # write out new aligned wdsp files
        with open('align_wdsp_into_msa.meg', 'w') as w_f:
            print >> w_f, '#mega'
            print >> w_f,'Title:WDSP Alignment'
            for entry in entries_new:
                print >> w_f,'#{0}'.format(entry[0])
                seq = ''.join(entry[1:])
                seqs = [seq[i:i+80] for i in range(0,len(seq),80)]
                for seq in seqs:
                    print >> w_f,seq


if __name__ == "__main__":
    main()

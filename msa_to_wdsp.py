#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
given a msa(fasta format) file and a wdsp file, change the msa file to wdsp file accordint to alignment

standard wdsp format is as follows:
> FBW1A_MOUSEQ3ULA2 1086.98 7
122 289 330 GRHSLQ RIHCRSETSKG VYCLQY DDQ KIVSGL RDN TIKIWD KSTL 44.0
145 334 370 ECKRIL TGHTGS VLCLQY DER VIITGS SDS TVRVWD VNAG 44.0
140 374 410 EMLNTL IHHCEA VLHLRF NNG MMVTCS KDR SIAVWD MASPTDI 44.0
144 457 493 TLRRVL VGHRAA VNVVDF DDK YIVSAS GDR TIKVWN TSTC 44.0
144 457 493 EFVRTL NGHKRG IACLQY RDR LVVSGS SDN TIRLWD IECG 44.0
129 497 533 ACLRVL EGHEEL VRCIRF DNK RIVSGA YDG KIKVWD LMAALDPRAPAGT 28.0
141 546 582 LCLRTL VEHSGR VFRLQF DEF QIVSSS HDD TILIWD 44.0

usage: python msa_to_wdsp *.wdsp *.fasta
"""

import os
import sys

def read_wdsp(wdsp_f):
    with open(wdsp_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        blades = [line.split() for line in lines]
        print blades
        return blades

def read_msa(msa_f):
    # readin seqs in fasta format
    # seqs foramt:[(pro,seq),...]
    with open(msa_f,'r') as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        pro_line_num = [i for i, line in enumerate(
            lines) if '>' in line] + [len(lines)]
        seqs = [lines[n:pro_line_num[i + 1]]
                for i, n in enumerate(pro_line_num[:-1])]
        seqs = [(seq[0][1:], ''.join(seq[1:])) for seq in seqs]
        return seqs

def match(blades,seq):
    blades_index = []
    pointer = 0
    for blade in blades:
        blade_index = blade[0:3]+[]
        for b in blade[3:-1]:
            b_index = []
            for bi in b:
                index = seq[pointer:].index(bi)
                b_index.append(pointer+index)
                pointer = pointer + index+1
            blade_index.append(b_index)

        try:
            blade_index.append(blade[-1])
        except:
            print blade
        blades_index.append(blade_index)

    # move extra aa to loop instead of sheet
    for blade_index in blades_index:
        for i in range(1,len(blade_index[3:-2]),2):
            if blade_index[3+i][0] - blade_index[3+i-1][-1] > 1:
                blade_index[3+i][0] = blade_index[3+i-1][-1] + 1
            if blade_index[3+i][-1] - blade_index[3+i+1][0] < -1:
                blade_index[3+i][-1] = blade_index[3+i+1][0] - 1
    for i in range(len(blades_index)-1):
        if blades_index[i][-2][-1] - blades_index[i+1][3][0] > 1:
            blades_index[i][-2][-1] = blades_index[i+1][3][0] + 1

    return blades_index


def write_wdsp(blades,seqs,filename):
    proid,score,repeat_num = blades[0][1:]
    match_seq = ''
    for pro,seq in seqs:
        if proid in pro:
            match_seq = seq
            break
    blades_index = match(blades[1:],match_seq)

    with open(filename+'_adjuted.wdsp','w') as w_f:
        for pro,seq in seqs:
            pro = pro.split()[0]
            print >> w_f,'> {0} {1} {2}'.format(pro,score,repeat_num)
            for blade in blades_index:
                new_blade = blade[0:3] + [seq[b[0]:b[-1]+1] for b in blade[3:-1]] + [blade[-1]]
                new_blade = ' '.join(new_blade).replace('-','')
                print >> w_f,new_blade

def main():
    blades = read_wdsp(sys.argv[-2])
    seqs = read_msa(sys.argv[-1])
    filename = os.path.splitext(os.path.split(sys.argv[-1])[1])[0]
    write_wdsp(blades,seqs,filename)

main()

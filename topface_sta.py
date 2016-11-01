#!/usr/bin/env python
# -*- coding: utf-8 -*-

import lt
import sys

def trans_topface():
    topface = lt.pickle_load(t_f)
    basic = ['K','R','H']
    acid = ['D','E']
    aromatic = ['F','W','Y']
    polar =['S','T','C','P','N','Q']
    nonpolar = ['G','V','L','I','M','A']
    res_hash = {'K':'b','R':'b','H':'b','D':'a','E':'a','S':'p','T':'p','C':'p','P':'p',\
            'N':'p','Q':'p','G':'n','V':'n','L':'n','I':'n','M':'n','A':'n','*':'*'}
    transed_topface = {}
    for p,h in topface:
        h_new = [res_hash[i] for bla in h for i in bla]
        h_new = [''.join(h_new[i:i+3]) for i in range(0,len(h_new),3)]
        transed_topface[p] = h_new


def hotsta():
    hot = lt.pickle_load(sys.argv[-1])
    r1 =        {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0 }
    r2 =        {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0 }
    d1 =        {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0 }
    r1_r2 =     {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0 }
    r1_d1 =     {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0 }
    r2_d1 =     {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0 }
    r1_r2_d1 =  {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0 }
    for pro_name,pro_hot in hot.iteritems():
        for h in pro_hot:
            if h[0] in r1.keys():
                r1[h[0]] += 1
            if h[1] in r2.keys():
                r2[h[1]] += 1
            if h[2] in d1.keys():
                d1[h[2]] += 1
    for k,v in r1.items():
        r1_r2[k] = r1[k] + r2[k]
        r2_d1[k] = r2[k] + d1[k]
        r1_d1[k] = r1[k] + d1[k]
        r1_r2_d1[k] = r1[k] + r2[k] +d1[k]
    fname = lt.fname(sys.argv[-1])
    lt.write_sta_dict(r1,fname + 'R1_hotsta')
    lt.write_sta_dict(r2,fname + 'R1_2_hotsta')
    lt.write_sta_dict(d1,fname + 'D_1_hotsta')
    lt.write_sta_dict(r1_d1,fname + 'R1_D_1_hotsta')
    lt.write_sta_dict(r2_d1,fname + 'R1_2_D_1_hotsta')
    lt.write_sta_dict(r1_r2,fname + 'R1_R1_2_hotsta')
    lt.write_sta_dict(r1_r2_d1,fname + 'R1_R1_2_D_1_hotsta')

def hotsta_classify():
    seq =  {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0 }
    hot = lt.pickle_load(sys.argv[-1])
    basic = ['K','R','H']
    acid = ['D','E']
    aromatic = ['F','W','Y']
    polar =['S','T','C','P','N','Q']
    nonpolar = ['G','V','L','I','M','A']
    r1 = {'basic':0,'acid':0,'aromatic':0,'polar':0,'nonpolar':0}
    r2 = {'basic':0,'acid':0,'aromatic':0,'polar':0,'nonpolar':0}
    d1 = {'basic':0,'acid':0,'aromatic':0,'polar':0,'nonpolar':0}
    sta = {'basic':0,'acid':0,'aromatic':0,'polar':0,'nonpolar':0}
    for pro_name,pro_hot in hot.iteritems():
        for h in pro_hot:
            if h[0] in basic:
                r1['basic'] += 1
                sta['basic'] += 1
            elif h[0] in acid:
                r1['acid'] += 1
                sta['acid'] += 1
            elif h[0] in aromatic:
                r1['aromatic'] += 1
                sta['aromatic'] += 1
            elif h[0] in polar:
                r1['polar'] += 1
                sta['polar'] += 1
            elif h[0] in nonpolar:
                r1['nonpolar'] += 1
                sta['nonpolar'] += 1
            if h[1] in basic:
                r2['basic'] += 1
                sta['basic'] += 1
            elif h[1] in acid:
                r2['acid'] += 1
                sta['acid'] += 1
            elif h[1] in aromatic:
                r2['aromatic'] += 1
                sta['aromatic'] += 1
            elif h[1] in polar:
                r2['polar'] += 1
                sta['polar'] += 1
            elif h[1] in nonpolar:
                r2['nonpolar'] += 1
                sta['nonpolar'] += 1
            if h[2] in basic:
                d1['basic'] += 1
                sta['basic'] += 1
            elif h[2] in acid:
                d1['acid'] += 1
                sta['acid'] += 1
            elif h[2] in aromatic:
                d1['aromatic'] += 1
                sta['aromatic'] += 1
            elif h[2] in polar:
                d1['polar'] += 1
                sta['polar'] += 1
            elif h[2] in nonpolar:
                d1['nonpolar'] += 1
                sta['nonpolar'] += 1
    fname = lt.fname(sys.argv[-1])
    lt.write_sta_dict(r1,fname + 'R1_hot_rc_sta')
    lt.write_sta_dict(r2,fname + 'R1_2_hot_rc_sta')
    lt.write_sta_dict(d1,fname + 'D_1_hot_rc_sta')
    lt.write_sta_dict(sta,fname + 'R1_R1_2_D_1_hot_rc_sta')


hotsta_classify()


#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import cPickle as pickle
from functools import wraps
from multiprocessing import Pool
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


gap_open = -11
gap_extend = -1.5

def gap_function(begin,length):
    if length == 0:
        return 0
    elif length == 1:
        return gap_open
    else:
        return gap_open + gap_extend*(length-1)

def wrap_match(blosum):
    def match_function(triple1,triple2):
        matrix = blosum
        score = []
        for a1,a2 in zip(triple1[0],triple2[0]):
            if matrix.has_key((a1,a2)):
                score.append(matrix[(a1,a2)])
            else:
                score.append(matrix[(a2,a1)])
        return sum(score)
    return match_function


def local_align_triple(p):
    blosums = [matlist.blosum30,matlist.blosum35,matlist.blosum40,matlist.blosum45, \
               matlist.blosum50,matlist.blosum55,matlist.blosum60,matlist.blosum62, \
               matlist.blosum65,matlist.blosum70,matlist.blosum75,matlist.blosum80, \
               matlist.blosum85,matlist.blosum90,matlist.blosum95,matlist.blosum100]

    i1,i2,seq1,seq2 = p
    seq1 = [[seqi] for seqi in seq1]
    seq2 = [[seqi] for seqi in seq2]

    max_identity = 0
    max_align = [[],[]]
    for blosum in blosums:
        alns = pairwise2.align.localcc(seq1, seq2,wrap_match(blosum),gap_function,gap_function,gap_char=['---'])
        inner_identity = 0
        inner_align = []
        for aln in alns:
            align1 = [a if not isinstance(a,list) else a[0] for a in aln[0] ]
            align2 = [a if not isinstance(a,list) else a[0] for a in aln[1] ]
            align1 = ''.join(align1)
            align2 = ''.join(align2)
            identical_res = [1 for si,s in enumerate(align1) if s == align2[si]]
            identity = 1.0*len(identical_res)/len(align1)
            align1 = [align1[i:i+3] for i in range(0,len(align1),3) ]
            align2 = [align2[i:i+3] for i in range(0,len(align2),3) ]
            if identity > inner_identity:
                inner_identity = identity
                inner_align = [align1,align2]
        if inner_identity > max_identity:
            max_identity = inner_identity
            max_align = inner_align

    return i1,i2,max_identity,max_align[0],max_align[1]


# p = [1,2,['TYT','HGY'],['TYT']]
# local_align_triple(p)
# sys.exit()

def global_align_triple(p):
    blosums = [matlist.blosum30,matlist.blosum35,matlist.blosum40,matlist.blosum45, \
               matlist.blosum50,matlist.blosum55,matlist.blosum60,matlist.blosum62, \
               matlist.blosum65,matlist.blosum70,matlist.blosum75,matlist.blosum80, \
               matlist.blosum85,matlist.blosum90,matlist.blosum95,matlist.blosum100]

    i1,i2,seq1,seq2 = p
    seq1 = [[seqi] for seqi in seq1]
    seq2 = [[seqi] for seqi in seq2]

    max_identity = 0
    max_align = [[],[]]
    for blosum in blosums:
        alns = pairwise2.align.globalcc(seq1, seq2,wrap_match(blosum),gap_function,gap_function,gap_char=['---'])
        inner_identity = 0
        inner_align = []
        for aln in alns:
            align1 = [a if not isinstance(a,list) else a[0] for a in aln[0] ]
            align2 = [a if not isinstance(a,list) else a[0] for a in aln[1] ]
            align1 = ''.join(align1)
            align2 = ''.join(align2)
            identical_res = [1 for si,s in enumerate(align1) if s == align2[si]]
            identity = 1.0*len(identical_res)/len(align1)
            align1 = [align1[i:i+3] for i in range(0,len(align1),3) ]
            align2 = [align2[i:i+3] for i in range(0,len(align2),3) ]
            if identity > inner_identity:
                inner_identity = identity
                inner_align = [align1,align2]
        if inner_identity > max_identity:
            max_identity = inner_identity
            max_align = inner_align

    return i1,i2,max_identity,max_align[0],max_align[1]


def local_align_hot(p):
    i1,i2,hot1,hot2 = p
    h1 = ''.join(hot1)
    h2 = ''.join(hot2)
    try:
        alns = pairwise2.align.localxx(h1, h2)
        align1 = alns[0][0]
        align2 = alns[0][1]
        identical_res = [1 for hi,h in enumerate(align1) if h == align2[hi]]
        identity = 1.0*len(identical_res)/len(align1)
        align1 = [align1[i:i+3] for i in range(0,len(align1),3)]
        align2 = [align2[i:i+3] for i in range(0,len(align2),3)]
        return i1,i2,identity,align1,align2
    except:
        return i1,i2,0,[],[]



def global_align_hot(p):
    i1,i2,hot1,hot2 = p
    h1 = ''.join(hot1)
    h2 = ''.join(hot2)
    try:
        alns = pairwise2.align.globalxx(h1, h2)
        align1 = alns[0][0]
        align2 = alns[0][1]
        identical_res = [1 for hi,h in enumerate(align1) if h == align2[hi]]
        identity = 1.0*len(identical_res)/len(align1)
        align1 = [align1[i:i+3] for i in range(0,len(align1),3)]
        align2 = [align2[i:i+3] for i in range(0,len(align2),3)]
        return i1,i2,identity,align1,align2
    except:
        return i1,i2,0,[],[]

def full_align(p):
    label,name,h1,h2 = p
    if label == 0:
        _,_,local_identity,local_align1,local_align2 = local_align_triple([0,1,h1,h2])
        return 0,local_identity,local_align1,local_align2
    elif label == 1:
        _,_,global_identity,global_align1,global_align2 = global_align_triple([0,1,h1,h2])
        return 1,global_identity,global_align1,global_align2
    elif label == 2:
        _,_,xx_local_identity,xx_local_align1,xx_local_align2 = local_align_hot([0,1,h1,h2])
        return 2,xx_local_identity,xx_local_align1,xx_local_align2
    elif label == 3:
        _,_,xx_global_identity,xx_global_align1,xx_global_align2 = global_align_hot([0,1,h1,h2])
        return 3,xx_global_identity,xx_global_align1,xx_global_align2



import lt
@lt.run_time
def main():
    fname = os.path.split(sys.argv[-1])[1].split('.')[0]
    with open(sys.argv[-1]) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        hotgroup = [lines[i:i+3] for i in range(0,len(lines),3)]
        hotgroup = [[h[0][1:],h[1].split(),h[2].split()] for h in hotgroup]

    with open(fname+'_compare_align.txt','w') as w_f:
        compare = [[],[],[],[]]
        for name,h1,h2 in hotgroup:
            # parameters = [[0,name,h1,h2],[1,name,h1,h2],[2,name,h1,h2],[3,name,h1,h2]]
            # p = Pool(4)
            # result = p.map(full_align,parameters)
            # p.close()
            # for r in result:
                # if r[0] == 0:
                    # _,local_identity,local_align1,local_align2 = r
                # elif r[0] == 1:
                    # _,global_identity,global_align1,global_align2 = r
                # elif r[0] == 2:
                    # _,xx_local_identity,xx_local_align1,xx_local_align2 = r
                # elif r[0] == 3:
                    # _,xx_global_identity,xx_global_align1,xx_global_align2 = r

            _,_,local_identity,local_align1,local_align2 = local_align_triple([0,1,h1,h2])
            _,_,global_identity,global_align1,global_align2 = global_align_triple([0,1,h1,h2])
            _,_,xx_local_identity,xx_local_align1,xx_local_align2 = local_align_hot([0,1,h1,h2])
            _,_,xx_global_identity,xx_global_align1,xx_global_align2 = global_align_hot([0,1,h1,h2])

            max_identity = max([local_identity,global_identity,xx_local_identity,xx_global_identity])

            good_align = []
            if max_identity == local_identity:
                compare[0].append(local_identity)
                good_align.append('localcc')
            if max_identity == global_identity:
                compare[1].append(global_identity)
                good_align.append('globalcc')

            if max_identity == xx_local_identity:
                compare[2].append(xx_local_identity)
                good_align.append('localxx')
            if max_identity == xx_global_identity:
                compare[3].append(xx_global_identity)
                good_align.append('globalxx')

            print >> w_f,name,good_align
            print >> w_f,name,local_identity,global_identity,xx_local_identity,xx_global_identity
            print >> w_f,h1
            print >> w_f,h2
            print >> w_f,'localcc'
            print >> w_f,local_align1
            print >> w_f,local_align2
            print >> w_f,'globalcc'
            print >> w_f,global_align1
            print >> w_f,global_align2
            print >> w_f,'localxx'
            print >> w_f,xx_local_align1
            print >> w_f,xx_local_align2
            print >> w_f,'globalxx'
            print >> w_f,xx_global_align1
            print >> w_f,xx_global_align2
            print >> w_f,''

            print >> w_f,''
            print >> w_f,''

        pickle.dump(compare,open(fname+'_compare.pickle','w'))
        compare = pickle.load(open(fname+'_compare.pickle'))

        sns.set_context('paper',font_scale=2)
        sns.set_style('white')
        plt.rc('text',usetex=False)
        fig,ax = plt.subplots(figsize=(6,6))
        sns.despine(left=True)
        labels = ['localcc','globalcc','localxx','globalxx']
        ax.hist(compare,histtype='bar',align='mid',label=labels)

        # sns.distplot(compare[0],hist=False,label='localcc')
        # sns.distplot(compare[1],hist=False,label='globalcc')
        # sns.distplot(compare[2],hist=False,label='localxx')
        # sns.distplot(compare[3],hist=False,label='globalxx')
        ax.set_xlabel('Sequence Identity')
        ax.get_yaxis().set_visible(False)
        ax.legend()
        plt.savefig(fname+'compare_align_dist.png',dpi=300)

        print >> w_f, map(len,compare)
        print >> w_f, 'total',len(hotgroup)
        print map(len,compare)
        print 'total',len(hotgroup)

if __name__ == "__main__":
    main()


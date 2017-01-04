#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
simulate a set of WD40s according to a template wd40
mutate a residue according to a similarity cutoff, indels not in consideration
calculate top_seq linearity
"""

"""
> Q969H0
KSPKVLKGHDDHVITCLQFCGNRIVSGSDDNTLKVWSAVTGKCLRTLVGHTGGVWSSQMRDNIIISGSTDRTLKVWNAETGECIHTLYGHTSTVRCMHLHEKRVVSGSRDATLRVWDIETGQCLHVLMGHVAAVRCVQYDGRRVVSGAYDFMVKVWDPETETCLHTLQGHTNRVYSLQFDGIHVVSGSLDTSIRVWDVETGNCIHTLTGHQSLTSGMELKDNILVSGNADSTVKIWDIKTGQCLQTLQGPNKHQSAVTCLQFNKNFVITSSDDGTVKLWDLKTGEFIRNLVTLESGGSGGVVWRIRASNTKLVCAVGSRNGTEETKLLVLD
"""
import sys
import os
import numpy as np
from numpy.random import randint
from copy import deepcopy
from scipy.stats import linregress
import matplotlib.pyplot as plt
import seaborn as sns

AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K','M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
TEM = 'KSPKVLKGHDDHVITCLQFCGNRIVSGSDDNTLKVWSAVTGKCLRTLVGHTGGVWSSQMRDNIIISGSTDRTLKVWNAETGECIHTLYGHTSTVRCMHLHEKRVVSGSRDATLRVWDIETGQCLHVLMGHVAAVRCVQYDGRRVVSGAYDFMVKVWDPETETCLHTLQGHTNRVYSLQFDGIHVVSGSLDTSIRVWDVETGNCIHTLTGHQSLTSGMELKDNILVSGNADSTVKIWDIKTGQCLQTLQGPNKHQSAVTCLQFNKNFVITSSDDGTVKLWDLKTGEFIRNLVTLESGGSGGVVWRIRASNTKLVCAVGSRNGTEETKLLVLD'
INDEX = [13,15,29,53,55,69,93,95,109,133,135,149,173,175,189,213,215,229,256,258,272,301,303,319]

def sim_one(similarity):
    new = list(deepcopy(TEM))
    for i in range(len(TEM)):
        if randint(1,101) > similarity:
            for j in range(200):
                mut = AA[randint(0,20)]
                if not mut == new[i]:
                    new[i] = mut
                    break
    return new

def calculate_similarity(sims):
    sim_len = len(sims)
    top_len = len(INDEX)
    tops = [[s[i-1] for i in INDEX] for s in sims]
    seq_score = []
    top_score = []
    for i in range(sim_len):
        for j in range(sim_len):
            if j > i:
                seq_score.append(sum(1 for i,s in enumerate(sims[i]) if s == sims[j][i])*1.0/sim_len)
                top_score.append(sum(1 for i,s in enumerate(tops[i]) if s == tops[j][i])*1.0/top_len)
    return top_score,seq_score


def sim_many(num,similarity,times=10000):
    regressions = []
    for i in range(times):
        sims = [sim_one(similarity) for i in range(num)]
        top_score,seq_score = calculate_similarity(sims)
        regressions.append(linregress(seq_score,top_score))
    slops = [r[0] for r in regressions]
    var = np.var(slops)
    mean = np.mean(slops)
    sns.distplot(slops)
    plt.savefig('sim_top_seq_'+str(num)+'_'+str(similarity)+'_'+str(mean)+'_'+str(var)+'.png',dpi=300)
    plt.close('all')
    return regressions

for num in range(10,200,10):
    for similarity in np.arange(0.3,1.0,0.1):
        sim_many(num,similarity)











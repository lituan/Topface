#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import operator
import itertools
from random import randint
from collections import OrderedDict

import numpy as np
# from matplotlib import pyplot as plt
# from sklearn import manifold
# from sklearn.decomposition import PCA
# from multiprocessing import freeze_support

import lt
from blosum import BLOSUM45, BLOSUM62
from wdsp import Wdsp

seq_matrix = BLOSUM62
UP = 1
LEFT = 2
DIAG = 0
INDEL = -7


def gap_penalty(i, j):
    return -3 * abs(i - j)


def score_hot(h1, h2):
    return sum([seq_matrix[h11][h22] for h11, h22 in zip(h1, h2)])


def needleman_wunsch_matrix(hot1, hot2, score_fun):
    """
    calculate the matrix
    """

    m, n = len(hot1), len(hot2)
    score_matrix = np.zeros((m + 1, n + 1))
    pointer_matrix = np.zeros((m + 1, n + 1), dtype=int)

    for i in range(1, m + 1):
        score_matrix[i, 0] = INDEL * i + (-3) * i * (i + 1) * 0.5
    for j in range(1, n + 1):
        score_matrix[0, j] = INDEL * j + (-3) * j * (j + 1) * 0.5

    pointer_matrix[0, 1:] = LEFT
    pointer_matrix[1:, 0] = UP

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            up = score_matrix[i - 1, j] + INDEL + gap_penalty(i, j)
            left = score_matrix[i, j - 1] + INDEL + gap_penalty(i, j)
            diag = score_matrix[i - 1, j - 1] + \
                score_fun(hot1[i - 1], hot2[j - 1])

            max_score = max(up, left, diag)
            if max_score == diag:
                pointer_matrix[i, j] = DIAG
            elif max_score == up:
                pointer_matrix[i, j] = UP
            elif max_score == left:
                pointer_matrix[i, j] = LEFT

            score_matrix[i, j] = max_score
    return score_matrix, pointer_matrix


def needleman_wunsch_trace(hot1, hot2, score_matrix, pointer_matrix):
    align1 = []
    align2 = []
    m, n = len(hot1), len(hot2)

    curr = pointer_matrix[m, n]
    i, j = m, n
    while (i > 0 or j > 0):
        if curr == DIAG:
            align1.append(hot1[i - 1])
            align2.append(hot2[j - 1])
            i -= 1
            j -= 1
        elif curr == LEFT:
            align1.append('---')
            align2.append(hot2[j - 1])
            j -= 1
        elif curr == UP:
            align1.append(hot1[i - 1])
            align2.append('---')
            i -= 1
        curr = pointer_matrix[i, j]

    return align1[::-1], align2[::-1]


def needleman_wunsch(hot1, hot2):
    score_matrix, pointer_matrix = needleman_wunsch_matrix(
        hot1, hot2, score_hot)
    hot1, hot2 = needleman_wunsch_trace(
        hot1, hot2, score_matrix, pointer_matrix)

    hot11 = ''.join(hot1)
    hot22 = ''.join(hot2)
    identity = int(200 * sum([1 for h1, h2 in zip(hot11, hot22)
                          if h1 == h2]) / (len(hot11) + len(hot22)))

    return identity

def check_hot_similar(hot1,hot2):
    return needleman_wunsch(hot1,hot2)


def get_similar_hots(tem_hots, all_hots, cutoff=0):
    tem_all_hots = {}
    for pro, hot in tem_hots.iteritems():
        similar = []
        i = 0
        for p, h in all_hots.iteritems():
            if not p == pro:
                identity = check_hot_similar(hot,h)
                if identity >= cutoff:
                    similar.append((p,identity))
                i +=1
                print i
        similar = sorted(similar,key=operator.itemgetter(1),reverse=True)
        tem_all_hots[pro] = similar

    return tem_all_hots


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
def seq_similarity(tem_wdsp,all_wdsp):
    tem_similarity = OrderedDict()
    for tem,tem_seq in tem_wdsp.seqs.iteritems():
        sim = OrderedDict()
        for pro,pro_seq in all_wdsp.seqs.iteritems():
            sim[pro] = align(tem_seq,pro_seq)
        tem_similarity[tem] = sim
    return tem_similarity


def repeat_similarity(repeats):
    lens = len(repeats)
    sims = []
    for i in xrange(lens):
        sim_i = []
        for j in xrange(lens):
            if j < i:
                sim_i.append(sims[j][i])
            elif j >= i:
                sim_i.append(align(repeats[i],repeats[j]))
        sims.append(sim_i)
    average = (sum([sum(i) for i in sims])-lens*100)/(lens*(lens-1))
    return average,sims


@lt.run_time
def wdsp_repeat_similarity(repeat_dic):
    sims = OrderedDict()
    for pro,repeats in repeat_dic.iteritems():
        sims[pro] = repeat_similarity(repeats)
    return sims


@lt.run_time
def write_result(tem_all_hots, tem_hots, tem_wdsp, all_hots, all_wdsp,tem_repeats_similarity,all_repeats_similarity,tem_all_seq_similarity,cutoff=0):
    for tem_pro, all_pros in tem_all_hots.iteritems():

        f_name = tem_pro + '_similar_hotspots_' + str(cutoff)
        with lt.open_file(f_name) as w_f:
            print >> w_f, '{0:<20}{1:<15}{2:<18}{3:<15}{4:<15}{5:<15}{6:<15}{7:<}'.format(
                'protein_id', 'identity','seq_similarity', "hotspot_num", 'repeat_length',"repeats_sim",'tetrad_num','hotspots' )
            print >> w_f, '{0:<20}{1:<15}{2:<18}{3:<15}{4:<15}{5:<15}{6:<15}{7:<}'.format(
                tem_pro, '100', '100', len(tem_hots[tem_pro]), tem_wdsp.repeat_num[tem_pro],tem_repeats_similarity[tem_pro][0],tem_wdsp.tetrad_num[tem_pro], ' '.join(tem_hots[tem_pro]))
            for pro,identity in all_pros:
                print >> w_f, '{0:<20}{1:<15}{2:<18}{3:<15}{4:<15}{5:<15}{6:<15}{7:<}'.format(
                    pro, identity, tem_all_seq_similarity[tem_pro][pro], len(all_hots[pro]), all_wdsp.repeat_num[pro],all_repeats_similarity[pro][0],all_wdsp.tetrad_num[pro], ' '.join(all_hots[pro]))

        f_name = tem_pro + '_similar_hotspots_wdsp' + str(cutoff)
        with lt.open_file(f_name,file_extension='.wdsp') as w_f:
            for line in tem_wdsp.wdsps[tem_pro]:
                print >> w_f, line
            for pro,identity in all_pros:
                for line in all_wdsp.wdsps[pro]:
                    print >> w_f,line

        f_name = tem_pro + '_similar_hotspots_seq' + str(cutoff)
        with lt.open_file(f_name) as w_f:
            print >> w_f,'>',tem_pro
            seq = tem_wdsp.seqs[tem_pro]
            for s in [seq[i:i+80] for i in range(0,len(seq),80)]:
                print >> w_f,s
            for pro in all_pros:
                print >> w_f,'>',pro[0]
                seq = all_wdsp.seqs[pro[0]]
                for s in [seq[i:i+80] for i in range(0,len(seq),80)]:
                    print >> w_f,s

@lt.run_time
def main():

    with open(sys.argv[-2]) as wdsp_f:
        tem_wdsp = Wdsp(wdsp_f)
        tem_hots = tem_wdsp.hotspots
        tem_repeats_similarity = wdsp_repeat_similarity(tem_wdsp.repeats)
        with open(sys.argv[-1]) as wdsp_f:
            all_wdsp = Wdsp(wdsp_f)
            all_hots = all_wdsp.hotspots
            all_repeats_similarity = wdsp_repeat_similarity(all_wdsp.repeats)
            tem_all_seq_similarity = seq_similarity(tem_wdsp,all_wdsp)


    cutoff = 70
    for cutoff in [30,40,50,60,70,80,90]:
        tem_all_hots = get_similar_hots(tem_hots, all_hots,cutoff)
        write_result(tem_all_hots, tem_hots, tem_wdsp, all_hots,all_wdsp,tem_repeats_similarity,all_repeats_similarity,tem_all_seq_similarity,cutoff)


if __name__ == "__main__":
    main()


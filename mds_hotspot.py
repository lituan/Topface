#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot mds plot for a a list of hotspots
input can be WDSP outfile or hotspots file
"""
import itertools
import sys
from random import randint

import numpy as np
from matplotlib import pyplot as plt
from sklearn import manifold
from sklearn.decomposition import PCA
from multiprocessing import freeze_support

import lt
from wdsp import Wdsp
from Bio.SubsMat import MatrixInfo

SEQ_MATRIX = MatrixInfo.blosum62
UP = 1
LEFT = 2
DIAG = 0
INDEL = -7


def gap_penalty(i, j):
    return -3 * abs(i - j)


def score_hot(h1, h2):
    score = []
    for h11, h22 in zip(h1, h2):
        if (h11, h22) in SEQ_MATRIX.keys():
            score.append(SEQ_MATRIX[(h11, h22)])
        elif (h22, h11) in SEQ_MATRIX.keys():
            score.append(SEQ_MATRIX[(h22, h11)])
        else:
            score.append(-5)

    return sum(score)


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

    return ho1,hot2

def align_hotspot(hot1,hot2):
    hot1,hot2 = needleman_wunsch(hot1,hot2)
    hot11 = ''.join(hot1)
    hot22 = ''.join(hot2)
    identity =  int(100*1.0*sum([1 for h1, h2 in zip(hot11, hot22) if h1 == h2]) / len(hot11))

    return identity, hot1, hot2

# a = ['ATT', 'BRR', 'GGY']
# b = ['ATT', 'BRR', 'BTT', 'GGG']
# c = ['FHN', 'GRR', 'CEN', 'KAA', 'CQY', 'LRG']
# d = ['SN*', 'GRC', 'CEN', 'KAA', 'CQY', 'LRA']
# e = ['QYG', 'SCD', 'GRR', 'CEN', 'KAT', 'CQY', 'LRG']
# print needleman_wunsch(d, b)

# input a list of hotspots, output distance matrix
def hot_scores_matrix(hotspots_list):
    len_hot = len(hotspots_list)
    scores = []
    for i in xrange(len_hot):
        score_i = []
        for j in xrange(len_hot):
            if j < i:
                score_i.append(scores[j][i])
            elif j == i:
                score_i.append(0)
            elif j > i:
                identity, _, _ = align_hotspot(
                    hotspots_list[i], hotspots_list[j])
                score_i.append(100 - identity)
        scores.append(score_i)
    return scores


def mds(scores_matrix):
    seed = np.random.RandomState(seed=3)
    mds = manifold.MDS(n_components=2, max_iter=3000000, eps=1e-9,
                        dissimilarity='precomputed', n_jobs=1)
    pos = mds.fit(scores_matrix).embedding_
    x = [p[0]*100 for p in pos]
    y = [p[1]*100 for p in pos]

    lt.pickle_dump(x,'pos_x')
    lt.pickle_dump(y,'pos_y')
    fig = plt.figure(1)
    plt.scatter(x, y)
    plt.savefig('cluster_mds.png')
    plt.close()

def read_hotspot(hot_f):
    lines = hot_f.readlines()
    hotspots = {}
    for line in lines:
        words = line.split()
        hotspots[words[0]] = words[1:]
    return hotspots.values()

def main():
    with open(sys.argv[-1]) as o_f:
        # input is WDSP output file
        # wdsp = Wdsp(o_f)
        # hotspots = wdsp.hotspots.values()

        # input is hotspots file
        hotspots = read_hotspot(o_f)

        scores = hot_scores_matrix(hotspots)
        mds(scores)

if __name__ == '__main__':
    freeze_support()
main()

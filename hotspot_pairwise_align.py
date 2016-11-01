#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

import numpy as np
from Bio.SubsMat import MatrixInfo

"""
input hots in the following format, output similarity matrix
pro xxx xxx xxx xxx xxx xxx

usuage example
python pairwise_align_hotspot.py example.hotspot
"""


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

    return hot1, hot2


def align_hotspot(hot1, hot2):
    return needleman_wunsch(hot1, hot2)

# d = ['SN*', 'GRC', 'CEN', 'KAA', 'CQY', 'LRA']
# e = ['QYG', 'SCD', 'GRR', 'CEN', 'KAT', 'CQY', 'LRG']
# print align_hotspot(d, e)


def align_lis_lis(lis_lis):
    """align and  nested list to print a table"""
    lis_lis = [[str(l) for l in lis]
               for lis in lis_lis]  # trans every element to str
    # make all inner lists of the same length
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [lis + (inner_lis_max_len - len(lis)) * [''] for lis in lis_lis]
    # trans list, so that the elements of the same column are in one list
    lis_lis = [[lis[i] for lis in lis_lis] for i in range(inner_lis_max_len)]
    # make element in the same list have the same length
    aligned = []
    for lis in lis_lis:
        width = max([len(l) for l in lis])
        lis = [l + (width - len(l)) * ' ' for l in lis]
        aligned.append(lis)
    # trans list_list to the original list_list
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [[lis[i] for lis in aligned] for i in range(inner_lis_max_len)]
    return lis_lis


def align(hot1, hot2):

    hot1, hot2 = align_hotspot(hot1, hot2)

    hot1 = ''.join(hot1)
    hot2 = ''.join(hot2)
    identity = [1 for i, s in enumerate(hot1) if s == hot2[i]]
    identity = 1.0 * len(identity) / len(hot1)

    return '{0:<12.10f}'.format(identity)

# d = ['SN*', 'GRC', 'CEN', 'KAA', 'CQY', 'LRA']
# e = ['QYG', 'SCD', 'GRR', 'CEN', 'KAT', 'CQY', 'LRG']
# print align(d,e)


def read_hot(hot_f):
    # readin hotspots
    # foramt:proname xxx xxx xxx xxx
    lines = hot_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    hots = [(line.split()[0], line.split()[1:]) for line in lines]

    return hots


def align_hots(hots):
    # hots format
    # [(proname,['xxx','xxx',...]),...]
    lens = len(hots)
    scores = []
    for i in xrange(lens):
        score_i = []
        for j in xrange(lens):
            # print i, '\t', j
            if j < i:
                score_i.append(scores[j][i])
            elif j >= i:
                score = align(hots[i][1], hots[j][1])
                score_i.append(score)
        scores.append(score_i)
    return scores


def write_resutls(hots, scores, file_path, file_name):
    result = [[hot[0]] + score for hot, score in zip(hots, scores)]
    header = [['ID'] + [hot[0] for hot in hots]]
    result = header + result

    # filename = os.path.join(file_path,file_name+'scores_tab.txt')
    # with open(filename,'w') as w_f:
    # for r in result:
    # print >> w_f, '\t'.join([str(ri)for ri in r])

    result = align_lis_lis(result)
    filename = os.path.join(file_path, file_name + 'scores_align.txt')
    with open(filename, 'w') as w_f:
        for r in result:
            print >> w_f, '\t'.join([str(ri)for ri in r])


def plot_heatmap(hots, scores):
    import matplotlib.pyplot as plt
    from numpy import array

    column_labels = [s[0] for s in hots]
    row_labels = column_labels
    scores = array(scores)

    fig, ax = plt.subplots()
    ax.axis('off')
    heatmap = ax.pcolor(scores, cmap=plt.cm.Blues)
    cb = plt.colorbar(heatmap)
    # ax.set_xticklabels(row_labels,minor=False)
    # ax.set_yticklabels(column_labels,minor=False)
    fig.savefig('test')


def main():
    with open(sys.argv[-1]) as hot_f:
        hots = read_hot(hot_f)
        scores = align_hots(hots)

    file_path, file_name = os.path.split(sys.argv[-1])
    file_name, file_extention = os.path.splitext(file_name)
    write_resutls(hots, scores, file_path, file_name)
    # plot_heatmap(hots,scores)

if __name__ == "__main__":
    main()

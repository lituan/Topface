#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
workflow
1. randomly select n pros from wdsp
2. compute top_seq
3. repeat step 1 and 2 10000 times
usage: python sim_top_seq.py all_nr.wdsp
"""

import sys
import os
import operator
import itertools
from numpy.random import randint
import numpy as np
import scipy
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as spd
from collections import OrderedDict
from scipy.stats import linregress
from wdsp import Wdsp

import igraph
from multiprocessing import Pool


def needleman_wunsch(hot1, hot2):
    from blosum import BLOSUM45, BLOSUM62
    seq_matrix = BLOSUM62
    UP,LEFT,DIAG,INDEL = 1,2,0,-7

    def gap_penalty(i, j):
        return -3 * abs(i - j)

    def score_hot(h1, h2):
        return sum([seq_matrix[h11][h22] for h11, h22 in zip(h1, h2)])

    def needleman_wunsch_matrix(hot1, hot2, score_fun):

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

    score_matrix, pointer_matrix = needleman_wunsch_matrix(
        hot1, hot2, score_hot)
    aligned_hot1, aligned_hot2 = needleman_wunsch_trace(
        hot1, hot2, score_matrix, pointer_matrix)

    return aligned_hot1,aligned_hot2

def align_hot(hot1,hot2):
    aligned_hot1,aligned_hot2 = needleman_wunsch(hot1,hot2)
    hot11 = ''.join(aligned_hot1)
    hot22 = ''.join(aligned_hot2)
    identity = 1.0 * sum([1 for h1, h2 in zip(hot11, hot22)
                          if h1 == h2]) / len(hot11)
    return identity

def align_hots(hots):
    score = []
    hots_len = len(hots)
    for i in range(hots_len):
        for j in range(hots_len):
            if j > i:
                score.append(align_hot(hots[i],hots[j]))
    return score

def align_seq(seq1, seq2):
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62
    gap_open = -10  # usual value
    gap_extend = -0.5  # usual value

    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)

    seq1 = alns[0][0]
    seq2 = alns[0][1]
    identity = [1 for i, s in enumerate(seq1) if s == seq2[i]]
    identity = 1.0 * len(identity) / len(seq1)

    return identity

def align_seqs(seqs):
    score = []
    seqs_len = len(seqs)
    for i in range(seqs_len):
        for j in range(seqs_len):
            if j > i:
                score.append(align_seq(seqs[i],seqs[j]))
    return score

def pair_to_matrix(pair):
    n = int((1+np.sqrt(len(pair)*8+1))/2.0)
    matrix = np.ones((n,n))
    for i in range(n):
        i_shift = i*n-i*(i+1)/2
        for j in range(n):
            j_shift = j*n-j*(j+1)/2
            if j > i:
                matrix[i][j] = pair[j-i-1+i_shift]
            if j < i:
                matrix[i][j] = pair[i-j-1+j_shift]
    return matrix

def matrix_to_pair(matrix):
    pair = []
    n = len(matrix)
    for i in range(n):
        for j in range(n):
            if j > i:
                pair.append(matrix[i][j])
    return pair


def single_fun(cluster):
    pros,label = cluster
    hots = [[pro,all_hots[pro]] for pro in pros]
    seqs = [[pro,all_wdsp.seqs[pro]] for pro in pros]
    hots_score = align_hots([hot[1] for hot in hots])
    seqs_score = align_seqs([seq[1] for seq in seqs])

    hots_score_matrix = pair_to_matrix(hots_score)
    seqs_score_matrix = pair_to_matrix(seqs_score)
    top_seq_matrix = np.array(hots_score_matrix) - np.array(seqs_score_matrix)
    top_seq_matrix = [[1 if e > 0.01 else 0 for e in line ] for line in top_seq_matrix]
    for i in range(len(top_seq_matrix)):
        top_seq_matrix[i][i]  = 0
    graph = igraph.Graph.Adjacency(top_seq_matrix,mode='undirected')
    # mc = graph.maximal_cliques()
    # mcm = sorted(mc,key=lambda x:len(x),reverse=True)[0]
    visual_style = {}
    # visual_style['vertex_label'] = labels
    # visual_style['vertex_label_size'] = 2
    visual_style['layout'] = graph.layout('kk')
    filename = 'sim_top_seq_'+str(cluster_num)+'_'+str(label)
    igraph.plot(graph,filename+'_graph.png',**visual_style)
    # slope,intercept,rvalue,pvalue,stderr = linregress(seqs_score,hots_score)
    # return [slope,intercept,rvalue,pvalue,stderr]
    # plot_scatter(seqs_score,hots_score,filename+'_scatter')

    # hots = adjust_hots(nr_hots)
    # hots = [(pro,''.join(hot)) for pro,hot in hots]
    # plotlogo(hots,filename+'_logo')

with open(sys.argv[-1]) as wdsp_f:
    all_wdsp = Wdsp(wdsp_f)
    all_hots = all_wdsp.hotspots
    all_seqs = all_wdsp.seqs
    pros = all_wdsp.pros

import lt
@lt.run_time
def main():



    for cluster_num in range(10,100,50):
        clusters = []
        for i in range(10):
            clusters.append([[pros[randint(0,len(pros))] for i in range(cluster_num)],i])

            p = Pool(10)
            p.map(single_fun,clusters)
            p.close()

if __name__ == "__main__":
    main()



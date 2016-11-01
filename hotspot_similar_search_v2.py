#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import itertools
import operator
from collections import defaultdict

import lt
from wdsp import Wdsp

"""
find similar hotspots
similarity cutoff is set by using blades identity
return the best score of all permutations:w
"""


def check_hot_similar(hot1, hot2, b, h):
    if len(hot1) < h or len(hot2) < h:
        return 0
    hot1h = [(hot1 + hot1)[i:i + h] for i in xrange(len(hot1))]
    hot2h = [(hot2 + hot2)[i:i + h] for i in xrange(len(hot2))]

    for h1 in hot1h:
        for h2 in hot2h:
            success = 1
            for h11, h22 in zip(h1, h2):
                common = [1 for r1, r2 in zip(h11, h22) if r1 == r2]
                if len(common) < b:
                    success = 0
                    break
            if success:
                return 1
    return 0


def get_similar_hots(tem_hots, all_hots, b, h):
    tem_all_hots = {}
    for pro, hot in tem_hots.iteritems():
        tem_all_hots[pro] = []
        for p, h in all_hots.iteritems():
            if check_hot_similar(hot, h, b, h):
                tem_all_hots[pro].append(p)

    return tem_all_hots


def write_result(tem_all_hots, tem_hots, all_hots, b, h, cutoff=1):
    for tem_pro, all_pros in tem_all_hots.iteritems():
        f_name = tem_pro + '_similar_hotspots_b' + str(b) + 'h' + str(h)
        with lt.open_file(f_name) as w_f:
            print >> w_f, '{0:<20}{1:<}'.format(
                tem_pro, ' '.join(tem_hots[tem_pro]))
            for pro in all_pros:
                print >> w_f, '{0:<20}{1:<}'.format(
                    pro, ' '.join(all_hots[pro]))


def read_hots(wdsp_f):
    with open(sys.argv[-2]) as wdsp_f:
        wdsp = Wdsp(wdsp_f)
        hots = wdsp.hotspots

        return hots


@lt.run_time
def main():
    tem_hots = read_hots(sys.argv[-2])
    all_hots = read_hots(sys.argv[-1])
    for b in 2, 3:
        for h in 2, 3, 4, 5, 6, 7:
            tem_all_hots = get_similar_hots(tem_hots, all_hots, b, h)
            write_result(tem_all_hots, tem_hots, all_hots, b, h)
    for b in [1]:

def score_hot(h1, h2):
    common = [1 for h11, h22 in zip(h1, h2) if h11 == h22]
    return len(common)


def best_fit_score(lis1, lis2, score_fun):
    """
    find bet_pairs between element in lis1 and lis2
    the idea is to calculate the score for each pair and
    get the best the pairs from the matrix
    """
    short_len = min(len(lis1), len(lis2))

    # calculate scores for each pair and sort pairs
    score_matrix = []
    for index_l1, l1 in enumerate(lis1):
        score_l1_lis2 = [(score_fun(l1, l2), index_l1, index_l2)
                         for index_l2, l2 in enumerate(lis2)]
        score_l1_lis2 = sorted(
            score_l1_lis2, key=operator.itemgetter(0), reverse=True)
        score_matrix.append(score_l1_lis2)

    # get best scores  from the matrix
    fit = []
    for i in range(short_len):
        max_score = [(score_l1_lis2[0]) for score_l1_lis2 in score_matrix]
        max_score = sorted(max_score, key=operator.itemgetter(0), reverse=True)
        fit.append(max_score[0])

        new_score_matrix = []
        for score_l1_lis2 in score_matrix:
            if score_l1_lis2[0][1] != fit[-1][1]:
                for l2_index, score in enumerate(score_l1_lis2):
                    if score[2] == fit[-1][2]:
                        score_l1_lis2.pop(l2_index)
                        new_score_matrix.append(score_l1_lis2)
        score_matrix = new_score_matrix

    pairs = [(f[0], lis1[f[1]], lis2[f[2]]) for f in fit]
    score = sum([f[0] for f in fit])
    return score, pairs

# a = ['ATT', 'BRR', 'GUU']
# b = ['ATT', 'BRR', 'BTT', 'GGG']
# print best_fit_score(b, a, score_hot)


def check_hot_similar(hot1, hot2):
    score, pairs = best_fit_score(hot1, hot2, score_hot)
    identity = score * 1.0 / ((len(hot1) + len(hot2)) * 0.5 * 3)
    return identity


def get_similar_hots(tem_hots, all_hots, cutoff=0.3):
    tem_all_hots = {}
    for pro, hot in tem_hots.iteritems():
        similar = []
        for p, h in all_hots.iteritems():
            identity = check_hot_similar(hot,h)
            if identity >= cutoff:
                similar.append((p,identity))
        similar = sorted(similar,key=operator.itemgetter(1),reverse=True)
        tem_all_hots[pro] = similar

    return tem_all_hots


def write_result(tem_all_hots, tem_hots, all_hots, cutoff=0.3):
    for tem_pro, all_pros in tem_all_hots.iteritems():
        f_name = tem_pro + '_similar_hotspots_' + str(cutoff)
        with lt.open_file(f_name) as w_f:
            print >> w_f, '{0:<20}{1:<10}{2:<}'.format(
                tem_pro, " ",' '.join(tem_hots[tem_pro]))
            for pro,identity in all_pros:
                print >> w_f, '{0:<20}{1:<10.2f}{2:<}'.format(
                    pro, identity, ' '.join(all_hots[pro]))


@lt.run_time
def main():
    with open(sys.argv[-2]) as wdsp_f:
        tem_wdsp = Wdsp(wdsp_f)
        tem_hots = tem_wdsp.hotspots
    with open(sys.argv[-1]) as wdsp_f:
        all_wdsp = Wdsp(wdsp_f)
        all_hots = all_wdsp.hotspots

    tem_all_hots = get_similar_hots(tem_hots, all_hots)
    write_result(tem_all_hots, tem_hots, all_hots)

if __name__ == "__main__":
    main()

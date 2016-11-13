# -*- coding:utf-8 -*-
import itertools
import operator
import os
import sys

import lt
from wdsp import Wdsp
from patch_search_base import PatchSearch,write_results


class PatchSearchSpecific(PatchSearch):

    def __init__(self, pros, seqs, wdsps, hotspots, cutoff=1):
        PatchSearch.__init__(self, pros, seqs, wdsps, hotspots, cutoff)

        self.shapes = {
            'shape_3_1_1_a':	((1, 2), (2, 0), (2, 1)),
            'shape_3_1_1_b':	((1, 0), (1, 2), (2, 1)),
            'shape_3_1_2_a':	((1, 2), (2, 0), (2, 2)),
            'shape_3_1_2_b':	((1, 0), (1, 1), (2, 1)),
            'shape_3_1_3_a':	((1, 0), (2, 1), (2, 2)),
            'shape_3_1_3_b':	((1, 1), (1, 2), (2, 0)),
            'shape_3_2_1_a':	((1, 0), (2, 0), (2, 1)),
            'shape_3_2_1_b':	((1, 0), (1, 2), (2, 0)),
            'shape_3_3':        ((1, 0), (1, 1), (1, 2)),
         }

        # self.phos_res = ['R','K','Y','W','Q','H','N','S','T']
        self.phos_res = ['R', 'K', 'Y', 'H', 'S', 'T']

    # def check_patches(self, shape_k, shape_v, patches, surround_patches):
        # #do not distinguish between 'S' and 'T'
        # patches = [''.join([p.replace('S','T') for p in patch]) for patch in patches]
        # # def validify(patch):
            # # for p in patch:
                # # if not p in self.phos_res:
                    # # return 0
            # # for i, p in enumerate(patch):
                # # if shape_v[i][1] == 0 and p in ['T', 'S']:
                    # # return 0
            # # return 1
        # # def validify_surround(surround_patch):
            # # if 'D' in surround_patch or 'E' in surround_patch:
                # # return 0
            # # return 1

        # def validify(patch):
            # return 1
        # def validify_surround(surround_patch):
            # return 1

        # good_patches = []
        # for patch,surround_patch in zip(patches,surround_patches):
            # if validify(patch) and validify_surround(surround_patch):
                # good_patches.append(patch)

        # return good_patches


@lt.run_time
def get_out_search():
    with open(sys.argv[-1]) as wdsp_f:

        CUTOFF = 20

        wdsp = Wdsp(wdsp_f)
        pros,seqs,wdsps,hotspots = wdsp.pros,wdsp.seqs,wdsp.wdsps,wdsp.hotspots

        best = {}
        pro_num = 1000000
        while pro_num > CUTOFF:
            a = PatchSearchSpecific(pros,seqs,wdsps,hotspots,CUTOFF)
            a.get_patches()
            a.classify_patches()
            shape,patch,pro_list,pro_num = a.get_best()
            if not shape in best.keys():
                best[shape] = {}
                best[shape][patch] = pro_list
            else:
                best[shape][patch] = pro_list

            for pro in pro_list:
                if pro in seqs.keys():
                    pros.pop(pros.index(pro))
                    seqs.pop(pro)
                    wdsps.pop(pro)
                    hotspots.pop(pro)

    with open(sys.argv[-1]) as wdsp_f:
        wdsp = Wdsp(wdsp_f)
        write_results(best,wdsp.pros,wdsp.seqs,wdsp.wdsps,wdsp.hotspots,CUTOFF)


@lt.run_time
def main():
    with open(sys.argv[1]) as wdsp_f:
        CUTOFF = 20
        wdsp = Wdsp(wdsp_f)
        a = PatchSearch(wdsp.pros,wdsp.seqs,wdsp.wdsps,wdsp.hotspots,CUTOFF)
        a.get_patches()
        a.classify_patches()
        write_results(a.shape_patch_pros,wdsp.pros,wdsp.seqs,wdsp.wdsps,wdsp.hotspots,CUTOFF)

get_out_search()

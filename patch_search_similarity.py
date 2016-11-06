#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
search patches on topface of wd40 proteins based on WDSP output file
"""
import itertools
import operator
import os
import sys

import lt
from wdsp import Wdsp
from patch_search_base import PatchSearch


class PatchSearchSpecific(PatchSearch):

    def __init__(self, pros, seqs, wdsps, hotspots, cutoff=1):
        PatchSearch.__init__(self, pros, seqs, wdsps, hotspots, cutoff)

        # self.phos_res = ['R','K','Y','W','Q','H','N','S','T']


    def check_patches(self, shape_k, shape_v, patches, surround_patches):
        patterns = ['shape_4_2_2_b','RRRY'] # FBXW7 patch
        def validify(patch):
            a = ['R','K']
            if shape_k == patterns[0] and patch[0] in a and patch[1] in a and patch[2] in a:
                return True

        def validify_surround(surround_patch):
            return 1

        good_patches = []
        for patch,surround_patch in zip(patches,surround_patches):
            if validify(patch) and validify_surround(surround_patch):
                good_patches.append(patch)

        return good_patches


@lt.run_time
def main():

    with open(sys.argv[-1]) as wdsp_f:
        wdsp = Wdsp(wdsp_f)
        a = PatchSearchSpecific(wdsp.pros,wdsp.seqs,wdsp.wdsps,wdsp.hotspots,cutoff=1)
        a.get_patches()
        a.classify_patches()
        a.write_results()

main()

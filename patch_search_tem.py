#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
find proteins having similar patch as template proteins

usage: python patch_search_tem.py tem_wdsp_f wdsp_f
"""
import os
import sys
import operator
from patch_search_base import PatchSearch
import lt


class PatchSearchTem(PatchSearch):

    def __init__(self, tem_f, wdsp_f, cutoff=1):
        PatchSearch.__init__(self, wdsp_f, cutoff)
        self.tem = []
        self.tem_pros = []

        lines = tem_f.readlines()
        lines = [line for line in lines if len(line.split()) > 0]
        lines = [line.strip('\n') for line in lines]
        self.tem = [(line.split()[0], line.split()[1]) for line in lines]

    def check_patches(self, shape_k, shape_v, patches, surround_patches):
        patches = [''.join([p.replace('S', 'T') for p in patch])
                   for patch in patches]

        def validify(patch):
            for p in patch:
                if not p in self.phos_res:
                    return 0
            for i, p in enumerate(patch):
                if shape_v[i][1] == 0 and p in ['T', 'S']:
                    return 0
            if patch.count('R') + patch.count('K') < 1:
                return 0
            return 1

        def validify_surround(surround_patch):
            if 'D' in surround_patch or 'E' in surround_patch:
                return 0
            return 1

        good_patches = []
        for patch, surround_patch in zip(patches, surround_patches):
            if validify(patch) and validify_surround(surround_patch):
                good_patches.append(patch)

        return good_patches

    def get_tem_pros(self):
        for shape_tem, patch_tem in self.tem:
            for shape, patch_pros in self.shape_patch_pros.iteritems():
                if shape == shape_tem:
                    for patch, pros in patch_pros.iteritems():
                        if patch == patch_tem:
                            self.tem_pros.append((shape_tem, patch_tem, pros))

    def write_tem_pros(self):
        with lt.open_file(file_suffix='tem_pros') as w_f:
            for shape, patch, pros in self.tem_pros:
                print >> w_f, '{0:<10}{1:<20}{2:<10}{3:<}'.format(
                    len(pros), shape, patch, ','.join(pros))


def main():
    tem_f = open(sys.argv[-2])
    wdsp_f = open(sys.argv[-1])
    a = PatchSearchTem(tem_f, wdsp_f, cutoff=1)
    a.get_patches()
    a.classify_patches()
    a.get_tem_pros()
    a.write_tem_pros()

if __name__ == "__main__":
    main()

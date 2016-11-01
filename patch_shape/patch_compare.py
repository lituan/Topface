#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
比较patch之间的包含关系
usage: python patch_f1 patch_f2
"""
import os
import sys
import lt
from collections import OrderedDict

SHAPES = {
    'shape_3_1_1_a':	((1, 2), (2, 0), (2, 1)),
    'shape_3_1_1_b':	((1, 0), (1, 2), (2, 1)),
    'shape_3_1_2_a':	((1, 2), (2, 0), (2, 2)),
    'shape_3_1_2_b':	((1, 0), (1, 1), (2, 1)),
    'shape_3_1_3_a':	((1, 0), (2, 1), (2, 2)),
    'shape_3_1_3_b':	((1, 1), (1, 2), (2, 0)),
    'shape_3_2_1_a':	((1, 0), (2, 0), (2, 1)),
    'shape_3_2_1_b':	((1, 0), (1, 2), (2, 0)),
    'shape_3_3':        ((1, 0), (1, 1), (1, 2)),

    'shape_4_1_1_a':	((1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_4_1_1_b':	((1, 0), (1, 1), (1, 2), (2, 1)),
    'shape_4_1_2_a':	((1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_4_1_2_b':	((1, 0), (1, 2), (2, 1), (2, 2)),
    'shape_4_2_1_a':	((1, 0), (1, 1), (1, 2), (2, 0)),
    'shape_4_2_1_b':	((1, 0), (2, 0), (2, 1), (2, 2)),
    'shape_4_2_2_a':	((1, 0), (1, 1), (2, 0), (2, 1)),
    'shape_4_2_2_b':	((1, 0), (1, 2), (2, 0), (2, 2)),
    'shape_4_4':	((1, 0), (1, 2), (2, 0), (2, 1)),

    'shape_5_2_1_a':	((1, 0), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_5_2_1_b':	((1, 0), (1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_5_2_2_a':	((1, 0), (1, 2), (2, 1), (2, 2), (3, 0)),
    'shape_5_2_2_b':	((1, 0), (2, 1), (2, 2), (3, 0), (3, 1)),
}

def read_patch(patch_f):
    with open(patch_f) as p_f:
        lines = p_f.readlines()
        return tuple([line.split()[:3] for line in lines])

def compare_patch(patch1,patch2):
    # patch1 is smaller than patch2
    shape_name_1 = patch1[1]
    shape_seq_1 = patch1[2]
    shape1 = SHAPES[shape_name_1]
    shape_name_2 = patch2[1]
    shape_seq_2 = patch2[2]
    shape2 = SHAPES[shape_name_2]

    if len(shape1) <= len(shape2):
        shape1 = [(s,shape_seq_1[i]) for i,s in enumerate(shape1)]
        shape2_1 = [(s,shape_seq_2[i]) for i,s in enumerate(shape2)]
        shape2_2 = [((s[0]-1,s[1]),shape_seq_2[i]) for i,s in enumerate(shape2)]

        if shape1 == [s for s in shape1 if s in shape2_1] or shape1 == [s for s in shape1 if s in shape2_2]:
            return True
        else:
            return False
    else:
        return False

# patch1 = [339,'shape_5_2_2_a','RRTYY']
# patch2 = [339,'shape_4_2_2_a','RRTY']
# compare_patch(patch2,patch1)

@lt.run_time
def main():
    patch3p = read_patch(sys.argv[-2])
    patch4p = read_patch(sys.argv[-1])


    dic_4_3 = OrderedDict()
    for p4 in patch4p:
        dic_4_3[tuple(p4)] = []
        for p3 in patch3p:
            if compare_patch(p3,p4):
                dic_4_3[tuple(p4)].append(p3)

    dic_3_4 = OrderedDict()
    for p3 in patch3p:
        dic_3_4[tuple(p3)] = []
        for p4 in patch4p:
            if compare_patch(p3,p4):
                dic_3_4[tuple(p3)].append(p4)

    with lt.open_file('dic_4_3') as w_f:
        lt.print_dic(dic_4_3,nest=0,output=w_f)

    with lt.open_file('dic_3_4') as w_f:
        lt.print_dic(dic_3_4,nest=0,output=w_f)


    non_unique_3 = set([tuple(s) for ss in dic_4_3.values() for s in ss])
    unique_3 = set([tuple(p) for p in patch3p]).difference(non_unique_3)
    unique_3 = [(int(u[0]),u[1],u[2]) for u in unique_3]
    unique_3 = sorted(unique_3,reverse=True)

    with lt.open_file('unique_3') as w_f:
        for p in unique_3:
            print >> w_f, p


if __name__ == "__main__":
    main()


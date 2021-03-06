#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
find good shapes
"""
import sys
import os
import itertools
import time
import random
import operator
import lt
import numpy as np
from collections import defaultdict

# all possible shapes
SHAPES = {
    #(0,0),first number indicates the blade,second number indicates position(0 for R1,
    # 1 for R1_1, 2 for D_1)
    # 3 positions
    # 1 R1
    # 1 2(0) or 2(0) 1
    'shape_3_1_a1_a': ((0, 1), (1, 0), (1, 1)),
    'shape_3_1_b2_a': ((0, 1), (1, 0), (1, 2)),
    'shape_3_1_c3_a': ((0, 2), (1, 0), (1, 1)),
    'shape_3_1_d4_a': ((0, 2), (1, 0), (1, 2)),
    'shape_3_1_d4_b': ((0, 0), (0, 1), (1, 1)),
    'shape_3_1_c3_b': ((0, 0), (0, 2), (1, 1)),
    'shape_3_1_b2_b': ((0, 0), (0, 1), (1, 2)),
    'shape_3_1_a1_b': ((0, 0), (0, 2), (1, 2)),
    'shape_3_1_e5_a': ((0, 0), (1, 1), (1, 2)),
    'shape_3_1_e5_b': ((0, 1), (0, 2), (1, 0)),
    # 3(0)
    'shape_3_1_f6': ((0, 0), (0, 1), (0, 2)),
    # 2 R1
    # 1(0) 2(0) or 2(0) 1(0)
    'shape_3_2_a1_a': ((0, 0), (1, 0), (1, 1)),
    'shape_3_2_b2_a': ((0, 0), (1, 0), (1, 2)),
    'shape_3_2_b2_b': ((0, 0), (0, 1), (1, 0)),
    'shape_3_2_a1_b': ((0, 0), (0, 2), (1, 0)),
    # 4 positions
    # 1 R1
    # 1 3(0) or 3(0) 1
    'shape_4_1_a1_a': ((0, 1), (1, 0), (1, 1), (1, 2)),
    'shape_4_1_b2_a': ((0, 2), (1, 0), (1, 1), (1, 2)),
    'shape_4_1_b2_b': ((0, 0), (0, 1), (0, 2), (1, 1)),
    'shape_4_1_a1_b': ((0, 0), (0, 1), (0, 2), (1, 2)),

    # 2 2(0) or 2(0) 2
    'shape_4_1_c3_a': ((0, 1), (0, 2), (1, 0), (1, 1)),
    'shape_4_1_d4_a': ((0, 1), (0, 2), (1, 0), (1, 2)),
    'shape_4_1_d4_b': ((0, 0), (0, 1), (1, 1), (1, 2)),
    'shape_4_1_c3_b': ((0, 0), (0, 2), (1, 1), (1, 2)),

    # 1 2(0) 1
    'shape_4_1_e5_a': ((0, 1), (1, 0), (1, 1), (2, 1)),
    'shape_4_1_f6_a': ((0, 1), (1, 0), (1, 1), (2, 2)),
    'shape_4_1_g7_a': ((0, 1), (1, 0), (1, 2), (2, 1)),
    'shape_4_1_h8_a': ((0, 2), (1, 0), (1, 1), (2, 1)),
    'shape_4_1_h8_b': ((0, 2), (1, 0), (1, 2), (2, 1)),
    'shape_4_1_g7_b': ((0, 2), (1, 0), (1, 1), (2, 2)),
    'shape_4_1_f6_b': ((0, 1), (1, 0), (1, 2), (2, 2)),
    'shape_4_1_e5_b': ((0, 2), (1, 0), (1, 2), (2, 2)),

    # 2(0) 1 1 or 1 1 2(0)
    'shape_4_1_i9_a': ((0, 0), (0, 1), (1, 1), (2, 1)),
    'shape_4_1_j10_a': ((0, 0), (0, 2), (1, 1), (2, 1)),
    'shape_4_1_k11_a': ((0, 0), (0, 1), (1, 2), (2, 1)),
    'shape_4_1_l12_a': ((0, 0), (0, 2), (1, 2), (2, 1)),
    'shape_4_1_m13_a': ((0, 0), (0, 1), (1, 1), (2, 2)),
    'shape_4_1_n14_a': ((0, 0), (0, 2), (1, 1), (2, 2)),
    'shape_4_1_o15_a': ((0, 0), (0, 1), (1, 2), (2, 2)),
    'shape_4_1_p16_a': ((0, 0), (0, 2), (1, 2), (2, 2)),

    'shape_4_1_p16_b': ((0, 1), (1, 1), (2, 0), (2, 1)),
    'shape_4_1_o15_b': ((0, 1), (1, 1), (2, 0), (2, 2)),
    'shape_4_1_n14_b': ((0, 1), (1, 2), (2, 0), (2, 1)),
    'shape_4_1_m13_b': ((0, 1), (1, 2), (2, 0), (2, 2)),
    'shape_4_1_l12_b': ((0, 2), (1, 1), (2, 0), (2, 1)),
    'shape_4_1_k11_b': ((0, 2), (1, 1), (2, 0), (2, 2)),
    'shape_4_1_j10_b': ((0, 2), (1, 2), (2, 0), (2, 1)),
    'shape_4_1_i9_b': ((0, 2), (1, 2), (2, 0), (2, 2)),

    # 2 1(0) 1 or 1 1(0) 2
    'shape_4_1_q17_a': ((0, 1), (0, 2), (1, 0), (2, 1)),
    'shape_4_1_r18_a': ((0, 1), (0, 2), (1, 0), (2, 2)),
    'shape_4_1_r18_b': ((0, 1), (1, 0), (2, 1), (2, 2)),
    'shape_4_1_q17_b': ((0, 2), (1, 0), (2, 1), (2, 2)),

    # 2 1 1(0) or 1(0) 1 2
    'shape_4_1_s19_a': ((0, 1), (0, 2), (1, 1), (2, 0)),
    'shape_4_1_t20_a': ((0, 1), (0, 2), (1, 2), (2, 0)),
    'shape_4_1_t20_b': ((0, 0), (1, 1), (2, 1), (2, 2)),
    'shape_4_1_s19_b': ((0, 0), (1, 2), (2, 1), (2, 2)),

    # 2 R1
    # 1(0) 3(0) or 3(0) 1(0)
    'shape_4_2_a1_a': ((0, 0), (0, 1), (0, 2), (1, 0)),
    'shape_4_2_a1_b': ((0, 0), (1, 0), (1, 1), (1, 2)),

    # 2(0) 2(0)
    'shape_4_2_b2_a': ((0, 0), (0, 1), (1, 0), (1, 1)),
    'shape_4_2_c3': ((0, 0), (0, 2), (1, 0), (1, 1)),
    'shape_4_2_d4': ((0, 0), (0, 1), (1, 0), (1, 2)),
    'shape_4_2_b2_b': ((0, 0), (0, 2), (1, 0), (1, 2)),

    # 1(0) 2(0) 1 or 1 2(0) 1(0)
    'shape_4_2_e5_a': ((0, 0), (1, 0), (1, 1), (2, 1)),
    'shape_4_2_f6_a': ((0, 0), (1, 0), (1, 2), (2, 1)),
    'shape_4_2_g7_a': ((0, 0), (1, 0), (1, 1), (2, 2)),
    'shape_4_2_h8_a': ((0, 0), (1, 0), (1, 2), (2, 2)),
    'shape_4_2_h8_b': ((0, 1), (1, 0), (1, 1), (2, 0)),
    'shape_4_2_g7_b': ((0, 1), (1, 0), (1, 2), (2, 0)),
    'shape_4_2_f6_b': ((0, 2), (1, 0), (1, 1), (2, 0)),
    'shape_4_2_e5_b': ((0, 2), (1, 0), (1, 2), (2, 0)),

    # 2(0) 1(0) 1 or 1 1(0) 2(0)
    'shape_4_2_i9_a': ((0, 0), (0, 1), (1, 0), (2, 1)),
    'shape_4_2_j10_a': ((0, 0), (0, 2), (1, 0), (2, 1)),
    'shape_4_2_k11_a': ((0, 0), (0, 1), (1, 0), (2, 2)),
    'shape_4_2_l12_a': ((0, 0), (0, 2), (1, 0), (2, 2)),
    'shape_4_2_l12_b': ((0, 1), (1, 0), (2, 0), (2, 1)),
    'shape_4_2_k11_b': ((0, 1), (1, 0), (2, 0), (2, 2)),
    'shape_4_2_j10_b': ((0, 2), (1, 0), (2, 0), (2, 1)),
    'shape_4_2_i9_b': ((0, 2), (1, 0), (2, 0), (2, 2)),

    # 2(0) 1 1(0) or 1(0) 1 2(0)
    'shape_4_2_m13_a': ((0, 0), (0, 1), (1, 1), (2, 0)),
    'shape_4_2_n14_a': ((0, 0), (0, 2), (1, 1), (2, 0)),
    'shape_4_2_o15_a': ((0, 0), (0, 1), (1, 2), (2, 0)),
    'shape_4_2_p16_a': ((0, 0), (0, 2), (1, 2), (2, 0)),
    'shape_4_2_p16_b': ((0, 0), (1, 1), (2, 0), (2, 1)),
    'shape_4_2_o15_b': ((0, 0), (1, 1), (2, 0), (2, 2)),
    'shape_4_2_n14_b': ((0, 0), (1, 2), (2, 0), (2, 1)),
    'shape_4_2_m13_b': ((0, 0), (1, 2), (2, 0), (2, 2)),

    # 2 1(0) 1(0) or 1(0) 1(0) 2
    'shape_4_2_q17_a': ((0, 1), (0, 2), (1, 0), (2, 0)),
    'shape_4_2_q17_b': ((0, 0), (1, 0), (2, 1), (2, 2)),

    # 3 R1
    # 2(0) 1(0) 1(0) or 1(0) 1(0) 2(0)
    'shape_4_3_a1_a': ((0, 0), (0, 1), (1, 0), (2, 0)),
    'shape_4_3_b2_a': ((0, 0), (0, 2), (1, 0), (2, 0)),
    'shape_4_3_b2_b': ((0, 0), (1, 0), (2, 0), (2, 1)),
    'shape_4_3_a1_b': ((0, 0), (1, 0), (2, 0), (2, 2)),

    # 1(0) 2(0) 1(0)
    'shape_4_3_c3_a': ((0, 0), (1, 0), (1, 1), (2, 0)),
    'shape_4_3_c3_b': ((0, 0), (1, 0), (1, 2), (2, 0)),

    # 5 positions
    # 1 R1
    # 1 3(0) 1
    'shape_5_1_a1_a': ((0, 1), (1, 0), (1, 1), (1, 2), (2, 1)),
    'shape_5_1_a1_b': ((0, 2), (1, 0), (1, 1), (1, 2), (2, 2)),
    'shape_5_1_b2': ((0, 2), (1, 0), (1, 1), (1, 2), (2, 1)),
    'shape_5_1_c3': ((0, 1), (1, 0), (1, 1), (1, 2), (2, 2)),

    # 2 3(0) or 3(0) 2
    'shape_5_1_d4_a': ((0, 1), (0, 2), (1, 0), (1, 1), (1, 2)),
    'shape_5_1_d4_b': ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2)),

    # 2 2(0) 1 or 1 2(0) 2
    'shape_5_1_g7_a': ((0, 1), (0, 2), (1, 0), (1, 1), (2, 1)),
    'shape_5_1_h8_a': ((0, 1), (0, 2), (1, 0), (1, 2), (2, 1)),
    'shape_5_1_i9_a': ((0, 1), (0, 2), (1, 0), (1, 1), (2, 2)),
    'shape_5_1_j10_a': ((0, 1), (0, 2), (1, 0), (1, 2), (2, 2)),
    'shape_5_1_j10_b': ((0, 1), (1, 0), (1, 1), (2, 1), (2, 2)),
    'shape_5_1_i9_b': ((0, 1), (1, 0), (1, 2), (2, 1), (2, 2)),
    'shape_5_1_h8_b': ((0, 2), (1, 0), (1, 1), (2, 1), (2, 2)),
    'shape_5_1_g7_b': ((0, 2), (1, 0), (1, 2), (2, 1), (2, 2)),

    # 2(0) 2 1 or 1 2 2(0)
    'shape_5_1_k11_a': ((0, 0), (0, 1), (1, 1), (1, 2), (2, 1)),
    'shape_5_1_l12_a': ((0, 0), (0, 2), (1, 1), (1, 2), (2, 1)),
    'shape_5_1_m13_a': ((0, 0), (0, 1), (1, 1), (1, 2), (2, 2)),
    'shape_5_1_n14_a': ((0, 0), (0, 2), (1, 1), (1, 2), (2, 2)),
    'shape_5_1_n14_b': ((0, 1), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_5_1_l12_b': ((0, 2), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_5_1_m13_b': ((0, 1), (1, 1), (1, 2), (2, 0), (2, 2)),
    'shape_5_1_k11_b': ((0, 2), (1, 1), (1, 2), (2, 0), (2, 2)),

    # 2 2 1(0) or 1(0) 2 2
    'shape_5_1_o15_a': ((0, 1), (0, 2), (1, 1), (1, 2), (2, 0)),
    'shape_5_1_o15_b': ((0, 0), (1, 1), (1, 2), (2, 1), (2, 2)),

    # 2 1(0) 2
    'shape_5_1_p16': ((0, 1), (0, 2), (1, 0), (2, 1), (2, 2)),
    # 2(0) 1 2 or 2 1 2(0)
    'shape_5_1_q17_a': ((0, 0), (0, 1), (1, 2), (2, 1), (2, 2)),
    'shape_5_1_r18_a': ((0, 0), (0, 2), (1, 1), (2, 1), (2, 2)),
    'shape_5_1_s19_a': ((0, 0), (0, 2), (1, 2), (2, 1), (2, 2)),
    'shape_5_1_t20_a': ((0, 0), (0, 1), (1, 1), (2, 1), (2, 2)),
    'shape_5_1_t20_b': ((0, 1), (0, 2), (1, 2), (2, 0), (2, 2)),
    'shape_5_1_s19_b': ((0, 1), (0, 2), (1, 1), (2, 0), (2, 1)),
    'shape_5_1_r18_b': ((0, 1), (0, 2), (1, 2), (2, 0), (2, 1)),
    'shape_5_1_q17_b': ((0, 1), (0, 2), (1, 1), (2, 0), (2, 2)),

    # 2 R1
    # 3(0) 2(0) or 2 (0) 3(0)
    'shape_5_2_a1_a': ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1)),
    'shape_5_2_b2_a': ((0, 0), (0, 1), (0, 2), (1, 0), (1, 2)),
    'shape_5_2_b2_b': ((0, 0), (0, 1), (1, 0), (1, 1), (1, 2)),
    'shape_5_2_a1_b': ((0, 0), (0, 2), (1, 0), (1, 1), (1, 2)),

    # 1 3(0) 1(0) or 1(0) 3(0) 1
    'shape_5_2_c3_a': ((0, 1), (1, 0), (1, 1), (1, 2), (2, 0)),
    'shape_5_2_d4_a': ((0, 2), (1, 0), (1, 1), (1, 2), (2, 0)),
    'shape_5_2_d4_b': ((0, 0), (1, 0), (1, 1), (1, 2), (2, 1)),
    'shape_5_2_c3_b': ((0, 0), (1, 0), (1, 1), (1, 2), (2, 2)),

    # 2(0) 2(0) 1 or 1 2(0) 2(0)
    'shape_5_2_e5_a': ((0, 1), (1, 0), (1, 1), (2, 0), (2, 1)),
    'shape_5_2_f6_a': ((0, 1), (1, 0), (1, 2), (2, 0), (2, 1)),
    'shape_5_2_g7_a': ((0, 2), (1, 0), (1, 2), (2, 0), (2, 1)),
    'shape_5_2_h8_a': ((0, 2), (1, 0), (1, 1), (2, 0), (2, 1)),
    'shape_5_2_i9_a': ((0, 1), (1, 0), (1, 1), (2, 0), (2, 2)),
    'shape_5_2_j10_a': ((0, 1), (1, 0), (1, 2), (2, 0), (2, 2)),
    'shape_5_2_k11_a': ((0, 2), (1, 0), (1, 1), (2, 0), (2, 2)),
    'shape_5_2_l12_a': ((0, 2), (1, 0), (1, 2), (2, 0), (2, 2)),
    'shape_5_2_l12_b': ((0, 0), (0, 1), (1, 0), (1, 1), (2, 1)),
    'shape_5_2_k11_b': ((0, 0), (0, 1), (1, 0), (1, 2), (2, 1)),
    'shape_5_2_j10_b': ((0, 0), (0, 1), (1, 0), (1, 1), (2, 2)),
    'shape_5_2_i9_b': ((0, 0), (0, 1), (1, 0), (1, 2), (2, 2)),
    'shape_5_2_h8_b': ((0, 0), (0, 2), (1, 0), (1, 2), (2, 1)),
    'shape_5_2_g7_b': ((0, 0), (0, 2), (1, 0), (1, 1), (2, 1)),
    'shape_5_2_f6_b': ((0, 0), (0, 2), (1, 0), (1, 1), (2, 2)),
    'shape_5_2_e5_b': ((0, 0), (0, 2), (1, 0), (1, 2), (2, 2)),

    # 2 2(0) 1(0) or 1(0) 2(0) 2
    'shape_5_2_m13_a': ((0, 1), (0, 2), (1, 0), (1, 1), (2, 0)),
    'shape_5_2_n14_a': ((0, 1), (0, 2), (1, 0), (1, 2), (2, 0)),
    'shape_5_2_n14_b': ((0, 0), (1, 0), (1, 1), (2, 1), (2, 2)),
    'shape_5_2_m13_b': ((0, 0), (1, 0), (1, 2), (2, 1), (2, 2)),

    # 2(0) 2 1(0) or 1(0) 2 2(0)
    'shape_5_2_o15_a': ((0, 0), (0, 1), (1, 1), (1, 2), (2, 0)),
    'shape_5_2_p16_a': ((0, 0), (0, 2), (1, 1), (1, 2), (2, 0)),
    'shape_5_2_p16_b': ((0, 0), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_5_2_o15_b': ((0, 0), (1, 1), (1, 2), (2, 0), (2, 2)),

    # 2(0) 1(0) 2 or 2 1(0) 2(0)
    'shape_5_2_q17_a': ((0, 0), (0, 1), (1, 0), (2, 1), (2, 2)),
    'shape_5_2_r18_a': ((0, 0), (0, 2), (1, 0), (2, 1), (2, 2)),
    'shape_5_2_r18_b': ((0, 1), (0, 2), (1, 0), (2, 0), (2, 1)),
    'shape_5_2_q17_b': ((0, 1), (0, 2), (1, 0), (2, 0), (2, 2)),

    # 2(0) 1 2(0)
    'shape_5_2_s19_a': ((0, 0), (0, 1), (1, 1), (2, 0), (2, 1)),
    'shape_5_2_t20_a': ((0, 0), (0, 2), (1, 1), (2, 0), (2, 1)),
    'shape_5_2_u21_a': ((0, 0), (0, 1), (1, 2), (2, 0), (2, 1)),
    'shape_5_2_v22_a': ((0, 0), (0, 1), (1, 1), (2, 0), (2, 2)),
    'shape_5_2_v22_b': ((0, 0), (0, 1), (1, 2), (2, 0), (2, 2)),
    'shape_5_2_u21_b': ((0, 0), (0, 2), (1, 1), (2, 0), (2, 2)),
    'shape_5_2_t20_b': ((0, 0), (0, 2), (1, 2), (2, 0), (2, 1)),
    'shape_5_2_s19_b': ((0, 0), (0, 2), (1, 2), (2, 0), (2, 2)),
    # 3 R1
    # 1(0) 3(0) 1(0)
    'shape_5_3_a1': ((0, 0), (1, 0), (1, 1), (1, 2), (2, 0)),

    # 2(0) 2(0) 1(0) or 1(0) 2(0) 2(0)
    'shape_5_3_b2_a': ((0, 0), (0, 1), (1, 0), (1, 1), (2, 0)),
    'shape_5_3_c3_a': ((0, 0), (0, 1), (1, 0), (1, 2), (2, 0)),
    'shape_5_3_d4_a': ((0, 0), (0, 2), (1, 0), (1, 1), (2, 0)),
    'shape_5_3_e5_a': ((0, 0), (0, 2), (1, 0), (1, 2), (2, 0)),
    'shape_5_3_e5_b': ((0, 0), (1, 0), (1, 1), (2, 0), (2, 1)),
    'shape_5_3_d4_b': ((0, 0), (1, 0), (1, 2), (2, 0), (2, 1)),
    'shape_5_3_c3_b': ((0, 0), (1, 0), (1, 1), (2, 0), (2, 2)),
    'shape_5_3_b2_b': ((0, 0), (1, 0), (1, 2), (2, 0), (2, 2)),

    # 2(0) 1(0) 2(0)
    'shape_5_3_f6_a': ((0, 0), (0, 1), (1, 0), (2, 0), (2, 1)),
    'shape_5_3_g7': ((0, 0), (0, 2), (1, 0), (2, 0), (2, 1)),
    'shape_5_3_h8': ((0, 0), (0, 1), (1, 0), (2, 0), (2, 2)),
    'shape_5_3_f6_b': ((0, 0), (0, 2), (1, 0), (2, 0), (2, 2)),

    # 6 positions
    # 1 R1
    # 1 2 3(0) or 3(0) 2 1
    'shape_6_1_a1_a': ((0, 1), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_6_1_b2_a': ((0, 2), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_6_1_b2_b': ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 1)),
    'shape_6_1_a1_b': ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)),

    # 2 3(0) 1 or 1 3(0) 2
    'shape_6_1_c3_a': ((0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 1)),
    'shape_6_1_d4_a': ((0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 2)),
    'shape_6_1_d4_b': ((0, 1), (1, 0), (1, 1), (1, 2), (2, 1), (2, 2)),
    'shape_6_1_c3_b': ((0, 2), (1, 0), (1, 1), (1, 2), (2, 1), (2, 2)),

    # 2 1 3(0) or 3(0) 1 2
    'shape_6_1_e5_a': ((0, 1), (0, 2), (1, 1), (2, 0), (2, 1), (2, 2)),
    'shape_6_1_f6_a': ((0, 1), (0, 2), (1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_6_1_f6_b': ((0, 0), (0, 1), (0, 2), (1, 1), (2, 1), (2, 2)),
    'shape_6_1_e5_b': ((0, 0), (0, 1), (0, 2), (1, 2), (2, 1), (2, 2)),

    # 2(0) 2 2 or 2 2 2(0)
    'shape_6_1_g7_a': ((0, 0), (0, 1), (1, 1), (1, 2), (2, 1), (2, 2)),
    'shape_6_1_h8_a': ((0, 0), (0, 2), (1, 1), (1, 2), (2, 1), (2, 2)),
    'shape_6_1_h8_b': ((0, 1), (0, 2), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_6_1_g7_b': ((0, 1), (0, 2), (1, 1), (1, 2), (2, 0), (2, 2)),
    # 2 2(0) 2
    'shape_6_1_i9_a': ((0, 1), (0, 2), (1, 0), (1, 1), (2, 1), (2, 2)),
    'shape_6_1_i9_b': ((0, 1), (0, 2), (1, 0), (1, 2), (2, 1), (2, 2)),

    # 2 R1
    # 3(0) 3(0)
    'shape_6_2_a1': ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)),

    # 1 2(0) 3(0) or 3(0) 2(0) 1
    'shape_6_2_b2_a': ((0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)),
    'shape_6_2_c3_a': ((0, 1), (1, 0), (1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_6_2_d4_a': ((0, 2), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)),
    'shape_6_2_e5_a': ((0, 2), (1, 0), (1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_6_2_e5_b': ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 1)),
    'shape_6_2_d4_b': ((0, 0), (0, 1), (0, 2), (1, 0), (1, 2), (2, 1)),
    'shape_6_2_c3_b': ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 2)),
    'shape_6_2_b2_b': ((0, 0), (0, 1), (0, 2), (1, 0), (1, 2), (2, 2)),

    # 1(0) 2 3(0)
    'shape_6_2_f6_a': ((0, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_6_2_f6_b': ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 0)),

    # 2(0) 3(0) 1 or 1 3(0) 2(0)
    'shape_6_2_g7_a': ((0, 0), (0, 1), (1, 0), (1, 1), (1, 2), (2, 1)),
    'shape_6_2_h8_a': ((0, 0), (0, 1), (1, 0), (1, 1), (1, 2), (2, 2)),
    'shape_6_2_i9_a': ((0, 0), (0, 2), (1, 0), (1, 1), (1, 2), (2, 1)),
    'shape_6_2_j10_a': ((0, 0), (0, 2), (1, 0), (1, 1), (1, 2), (2, 2)),
    'shape_6_2_j10_b': ((0, 1), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_6_2_i9_b': ((0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_6_2_h8_b': ((0, 1), (1, 0), (1, 1), (1, 2), (2, 0), (2, 2)),
    'shape_6_2_g7_b': ((0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 2)),

    # 2 3(0) 1(0)
    'shape_6_2_k11_a': ((0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0)),
    'shape_6_2_k11_b': ((0, 0), (1, 0), (1, 1), (1, 2), (2, 1), (2, 2)),

    # 2(0) 1 3(0) or 3(0) 1 2(0)
    'shape_6_2_l12_a': ((0, 0), (0, 1), (1, 1), (2, 0), (2, 1), (2, 2)),
    'shape_6_2_m13_a': ((0, 0), (0, 2), (1, 1), (2, 0), (2, 1), (2, 2)),
    'shape_6_2_n14_a': ((0, 0), (0, 1), (1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_6_2_o15_a': ((0, 0), (0, 2), (1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_6_2_o15_b': ((0, 0), (0, 1), (0, 2), (1, 1), (2, 0), (2, 1)),
    'shape_6_2_n14_b': ((0, 0), (0, 1), (0, 2), (1, 1), (2, 0), (2, 2)),
    'shape_6_2_m13_b': ((0, 0), (0, 1), (0, 2), (1, 2), (2, 0), (2, 1)),
    'shape_6_2_l12_b': ((0, 0), (0, 1), (0, 2), (1, 2), (2, 0), (2, 2)),

    # 2 1(0) 3(0)
    'shape_6_2_p16_a': ((0, 1), (0, 2), (1, 0), (2, 0), (2, 1), (2, 2)),
    'shape_6_2_p16_b': ((0, 0), (0, 1), (0, 2), (1, 0), (2, 1), (2, 2)),

    # 2(0) 2(0) 2 or 2 2(0) 2(0)
    'shape_6_2_q17_a': ((0, 0), (0, 1), (1, 0), (1, 1), (2, 1), (2, 2)),
    'shape_6_2_r18_a': ((0, 0), (0, 2), (1, 0), (1, 1), (2, 1), (2, 2)),
    'shape_6_2_s19_a': ((0, 0), (0, 1), (1, 0), (1, 2), (2, 1), (2, 2)),
    'shape_6_2_t20_a': ((0, 0), (0, 2), (1, 0), (1, 2), (2, 1), (2, 2)),
    'shape_6_2_t20_b': ((0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 1)),
    'shape_6_2_r18_b': ((0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)),
    'shape_6_2_s19_b': ((0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 2)),
    'shape_6_2_q17_b': ((0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 2)),

    # 2(0) 2 2(0)
    'shape_6_2_u21_a': ((0, 0), (0, 1), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_6_2_v22': ((0, 0), (0, 2), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_6_2_w23': ((0, 0), (0, 1), (1, 1), (1, 2), (2, 0), (2, 2)),
    'shape_6_2_u21_b': ((0, 0), (0, 2), (1, 1), (1, 2), (2, 0), (2, 2)),

    # 3 R1
    # 1(0) 2(0) 3(0) or 3(0) 2(0) 1(0)
    'shape_6_3_a1_a': ((0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)),
    'shape_6_3_b2_a': ((0, 0), (1, 0), (1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_6_3_b2_b': ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0)),
    'shape_6_3_a1_b': ((0, 0), (0, 1), (0, 2), (1, 0), (1, 2), (2, 0)),

    # 2(0) 3(0) 1(0) or 1(0) 3(0) 2(0)
    'shape_6_3_c3_a': ((0, 0), (0, 1), (1, 0), (1, 1), (1, 2), (2, 0)),
    'shape_6_3_d4_a': ((0, 0), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0)),
    'shape_6_3_d4_b': ((0, 0), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_6_3_c3_b': ((0, 0), (1, 0), (1, 1), (1, 2), (2, 0), (2, 2)),

    # 2(0) 1(0) 3 (0) or 3(0) 1(0) 2(0)
    'shape_6_3_e5_a': ((0, 0), (0, 1), (1, 0), (2, 0), (2, 1), (2, 2)),
    'shape_6_3_f6_a': ((0, 0), (0, 2), (1, 0), (2, 0), (2, 1), (2, 2)),
    'shape_6_3_f6_b': ((0, 0), (0, 1), (0, 2), (1, 0), (2, 0), (2, 1)),
    'shape_6_3_e5_b': ((0, 0), (0, 1), (0, 2), (1, 0), (2, 0), (2, 2)),

    # 2(0) 2(0) 2(0)
    'shape_6_3_g7_a': ((0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)),
    'shape_6_3_h8_a': ((0, 0), (0, 2), (1, 0), (1, 1), (2, 0), (2, 1)),
    'shape_6_3_i9_a': ((0, 0), (0, 1), (1, 0), (1, 2), (2, 0), (2, 1)),
    'shape_6_3_j10_a': ((0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 2)),
    'shape_6_3_j10_b': ((0, 0), (0, 1), (1, 0), (1, 2), (2, 0), (2, 2)),
    'shape_6_3_h8_b': ((0, 0), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)),
    'shape_6_3_i9_b': ((0, 0), (0, 2), (1, 0), (1, 1), (2, 0), (2, 2)),
    'shape_6_3_g7_b': ((0, 0), (0, 2), (1, 0), (1, 2), (2, 0), (2, 2)),
    'shape_6_3_bb': ((0, 0), (0, 2), (0, 0), (0, 2), (2, 0), (2, 2)),
}

SEQ = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
       'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X', '*']

PHOS_RES_EXTEND = ['R', 'K', 'Y', 'H', 'S', 'T', 'A', 'N', 'D', 'C', 'Q', 'E',
                   'G', 'I', 'L', 'M', 'F', 'P', 'W', 'V']  # residues than can interact with phos
PHOS_RES = ['R', 'K', 'Y', 'H', 'S', 'T']

# check if two blades in a shape are same
# if so, delete this shape and print it


def check_blade_identity(shapes):

    bad_shapes = {}
    for k, v in shapes.iteritems():
        for vi in v:
            if v.count(vi) >= 2:
                print 'blades are same:', '\t', k, '\t', v
                bad_shapes[k] = v
    for k in bad_shapes.keys():
        del shapes[k]

# calculate distance between two positions


def position_distance(p1, p2):
    if p1[0] == p2[0]:
        if p1[1] == p2[1]:
            return 0
        else:
            return 1
    elif p2[0] - p1[0] == 1:
        if p1[1] == p2[1]:
            if p1[1] == 0:
                return 1
            else:
                return 2
        elif p2[1] == 2 and p1[1] == 0:
            return 2
        elif p2[1] == 2 and p1[1] == 1:
            return 3
        elif p2[1] == 1 and p1[1] == 0:
            return 1
        elif p2[1] == 1 and p1[1] == 2:
            return 1
        elif p2[1] == 0 and p1[1] == 1:
            return 2
        elif p2[1] == 0 and p1[1] == 2:
            return 1
    elif p2[0] - p1[0] == 2:
        if p1[1] == p2[1]:
            if p1[1] == 0:
                return 2
            else:
                return 3
        elif p2[1] == 2 and p1[1] == 0:
            return 3
        elif p2[1] == 2 and p1[1] == 1:
            return 4
        elif p2[1] == 1 and p1[1] == 0:
            return 2
        elif p2[1] == 1 and p1[1] == 2:
            return 2
        elif p2[1] == 0 and p1[1] == 1:
            return 3
        elif p2[1] == 0 and p1[1] == 2:
            return 2
    elif p2[0] - p1[0] == -1:
        if p1[1] == p2[1]:
            if p1[1] == 0:
                return 1
            else:
                return 2
        elif p2[1] == 2 and p1[1] == 0:
            return 1
        elif p2[1] == 2 and p1[1] == 1:
            return 1
        elif p2[1] == 1 and p1[1] == 0:
            return 2
        elif p2[1] == 1 and p1[1] == 2:
            return 3
        elif p2[1] == 0 and p1[1] == 1:
            return 1
        elif p2[1] == 0 and p1[1] == 2:
            return 2
    elif p2[0] - p1[0] == -2:
        if p1[1] == p2[1]:
            if p1[1] == 0:
                return 2
            else:
                return 4
        elif p2[1] == 2 and p1[1] == 0:
            return 2
        elif p2[1] == 2 and p1[1] == 1:
            return 2
        elif p2[1] == 1 and p1[1] == 0:
            return 3
        elif p2[1] == 1 and p1[1] == 2:
            return 4
        elif p2[1] == 0 and p1[1] == 1:
            return 2
        elif p2[1] == 0 and p1[1] == 2:
            return 3


def check_position_distance(shapes):
    bad_shapes = {}
    for k, v in shapes.iteritems():
        distance = []
        for vi in v:
            v_d = []
            for vj in v:
                if vi != vj:
                    v_d.append(position_distance(vi, vj))
                    distance.append(position_distance(vi, vj))
            if min(v_d) > 1:
                if not k in bad_shapes.keys():
                    bad_shapes[k] = v
                    print 'poistions not near each other:', '\t', k, '\t', v, '\t', v_d
        if max(distance) >= 3:
            if not k in bad_shapes.keys():
                bad_shapes[k] = v
                print 'poistion distance too big:', '\t', k, '\t', v
    for k in bad_shapes.keys():
        del shapes[k]


# @lt.log('check_shape')
def check_shapes():
    shapes = SHAPES
    check_blade_identity(shapes)
    check_position_distance(shapes)

    ofile = open('good_shape.txt', 'w')
    keys = shapes.keys()
    keys = sorted(keys)
    for k in keys:
        v = shapes[k]
        print >> ofile, "'{0}':\t{1}".format(k, v)

check_shapes()

shapes = {
    'shape_3_1_1_a':	((0, 2), (1, 0), (1, 1)),
    'shape_3_1_1_b':	((0, 0), (0, 2), (1, 1)),
    'shape_3_1_2_a':	((0, 2), (1, 0), (1, 2)),
    'shape_3_1_2_b':	((0, 0), (0, 1), (1, 1)),
    'shape_3_1_3_a':	((0, 0), (1, 1), (1, 2)),
    'shape_3_1_3_b':	((0, 1), (0, 2), (1, 0)),
    'shape_3_2_1_a':	((0, 0), (1, 0), (1, 1)),
    'shape_3_2_1_b':	((0, 0), (0, 2), (1, 0)),
    'shape_3_2_2_a':	((0, 0), (1, 0), (1, 2)),
    'shape_3_2_2_b':	((0, 0), (0, 1), (1, 0)),
    'shape_3_3':	    ((0, 0), (0, 1), (0, 2)),
    'shape_4_1_1_a':	((0, 2), (1, 0), (1, 1), (1, 2)),
    'shape_4_1_1_b':	((0, 0), (0, 1), (0, 2), (1, 1)),
    'shape_4_1_2_a':	((0, 1), (0, 2), (1, 0), (1, 1)),
    'shape_4_1_2_b':	((0, 0), (0, 2), (1, 1), (1, 2)),
    'shape_4_1_3_a':	((0, 2), (1, 0), (1, 1), (2, 1)),
    'shape_4_1_3_b':	((0, 2), (1, 0), (1, 2), (2, 1)),
    'shape_4_2_1_a':	((0, 0), (0, 1), (0, 2), (1, 0)),
    'shape_4_2_1_b':	((0, 0), (1, 0), (1, 1), (1, 2)),
    'shape_4_2_2_a':	((0, 0), (0, 1), (1, 0), (1, 1)),
    'shape_4_2_2_b':	((0, 0), (0, 2), (1, 0), (1, 2)),
    'shape_4_2_4_a':	((0, 0), (1, 0), (1, 1), (2, 1)),
    'shape_4_2_4_b':	((0, 2), (1, 0), (1, 2), (2, 0)),
    'shape_4_2_5_a':	((0, 0), (1, 0), (1, 2), (2, 1)),
    'shape_4_2_5_b':	((0, 2), (1, 0), (1, 1), (2, 0)),
    'shape_4_2_6_a':	((0, 0), (0, 2), (1, 0), (2, 1)),
    'shape_4_2_6_b':	((0, 2), (1, 0), (2, 0), (2, 1)),
    'shape_4_3_1_a':	((0, 0), (0, 2), (1, 0), (2, 0)),
    'shape_4_3_1_b':	((0, 0), (1, 0), (2, 0), (2, 1)),
    'shape_4_3_2_a':	((0, 0), (1, 0), (1, 1), (2, 0)),
    'shape_4_3_2_b':	((0, 0), (1, 0), (1, 2), (2, 0)),
    'shape_4_4':	    ((0, 0), (0, 2), (1, 0), (1, 1)),
    'shape_5_1_1':	    ((0, 2), (1, 0), (1, 1), (1, 2), (2, 1)),
    'shape_5_1_2_a':	((0, 0), (0, 2), (1, 1), (1, 2), (2, 1)),
    'shape_5_1_2_b':	((0, 2), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_5_2_3_a':	((0, 0), (0, 1), (0, 2), (1, 0), (1, 1)),
    'shape_5_2_3_b':	((0, 0), (0, 2), (1, 0), (1, 1), (1, 2)),
    'shape_5_2_4_a':	((0, 2), (1, 0), (1, 1), (1, 2), (2, 0)),
    'shape_5_2_4_b':	((0, 0), (1, 0), (1, 1), (1, 2), (2, 1)),
    'shape_5_2_5_a':	((0, 2), (1, 0), (1, 2), (2, 0), (2, 1)),
    'shape_5_2_5_b':	((0, 0), (0, 2), (1, 0), (1, 1), (2, 1)),
    'shape_5_2_6_a':	((0, 2), (1, 0), (1, 1), (2, 0), (2, 1)),
    'shape_5_2_6_b':	((0, 0), (0, 2), (1, 0), (1, 2), (2, 1)),
    'shape_5_2_7_a':	((0, 0), (0, 2), (1, 1), (1, 2), (2, 0)),
    'shape_5_2_7_b':	((0, 0), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_5_3_2_a':	((0, 0), (0, 2), (1, 0), (1, 1), (2, 0)),
    'shape_5_3_2_b':	((0, 0), (1, 0), (1, 2), (2, 0), (2, 1)),
    'shape_5_3_3_a':	((0, 0), (0, 2), (1, 0), (1, 2), (2, 0)),
    'shape_5_3_3_b':	((0, 0), (1, 0), (1, 1), (2, 0), (2, 1)),
    'shape_5_5_1':	    ((0, 0), (1, 0), (1, 1), (1, 2), (2, 0)),
    'shape_5_5_2':	    ((0, 0), (0, 2), (1, 0), (2, 0), (2, 1)),
    'shape_6_2_1_a':	((0, 0), (0, 2), (1, 0), (1, 1), (1, 2), (2, 1)),
    'shape_6_2_1_b':	((0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_6_3_3_a':	((0, 0), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0)),
    'shape_6_3_3_b':	((0, 0), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_6_3_4_a':	((0, 0), (0, 2), (1, 0), (1, 1), (2, 0), (2, 1)),
    'shape_6_3_4_b':	((0, 0), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)),
    'shape_6_6':	    ((0, 0), (0, 2), (1, 1), (1, 2), (2, 0), (2, 1)),
}


def shape_contain(shapes):
    contain = {}
    for ki, vi in shapes.iteritems():
        contain[ki] = []
        for kj, vj in shapes.iteritems():
            contained = 1
            if ki != kj:
                for vjj in vj:
                    if not vjj in vi:
                        contained = 0
                if contained:
                    contain[ki].append(kj)
    ofile = open('good_shape_contain.txt', 'w')
    keys = contain.keys()
    keys = sorted(keys)
    for k in keys:
        v = contain[k]
        print >> ofile, "'{0}':\t{1}".format(k, v)

shape_contain(shapes)

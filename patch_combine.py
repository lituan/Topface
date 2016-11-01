#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import lt

with open(sys.argv[-2]) as o_f:
    lines = o_f.readlines()
    lines = [line.split() for line in lines if len(line.split()) > 0]
    search_1 = lines
with open(sys.argv[-1]) as o_f:
    lines = o_f.readlines()
    lines = [line.split() for line in lines if len(line.split()) > 0]
    search_2 = lines

search_3 = []

for s1 in search_1:
    s3 = s1[0:3]
    for s2 in search_2:
        if s1[1:3] == s2[1:3]:
            s3.extend([s2[-1]])
    if len(s3) == 4:
        pass
    else:
        s3.extend([' '])
    search_3.append(s3)

with open('combine.txt','w') as w_f:

    for len_pros, shape, patch, pros in search_3:
        pros = pros.split(',')
        pros = [p for p in pros if 'HUMAN' in p]
        human_num = len([p for p in pros])
        print >> w_f, '{0:<10}{1:<20}{2:<10}{3:<10}{4:<}'.format(
            len_pros, shape, patch,human_num,','.join(pros))

pros = [s[-1].split(',') for s in search_3]
pros = [pi for p in pros for pi in p if len(pi.split()) > 0]
pros_new = []
for p in pros:
    if not p in pros_new:
        pros_new.append(p)

with open('combine_pro_order.txt','w') as w_f:
    for p in pros_new:
        print >> w_f,p






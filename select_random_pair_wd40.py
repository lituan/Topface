#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
select wd40s with repeats 6n,7n,8n
"""
from wdsp import Wdsp
import numpy as np
from numpy.random import randint

with open('wd648_uniprot_select_cd-hit_90_7.wdsp') as wdsp_f:
    wdsp = Wdsp(wdsp_f)
    pro_num = len(wdsp.pros)
    for i in range(100):
        with open('random_'+str(i)+'.fa','w') as w_f:
            for j in range(2):
                pro = wdsp.pros[randint(0,pro_num)]
                seq = wdsp.seqs[pro]
                print >> w_f,'> ',pro
                for s in [seq[i:i+80] for i in range(0,len(seq),80)]:
                    print >> w_f,s


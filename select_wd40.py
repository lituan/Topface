#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
select wd648s with repeats 6n,7n,8n
"""
from wdsp import Wdsp

with open('wd648_uniprot.wdsp') as wdsp_f:
    wdsp = Wdsp(wdsp_f)

    with open('wd648_uniprot_select.wdsp','w') as w_f:
        for pro,pro_wdsp in wdsp.wdsps.iteritems():
            r_num = wdsp.repeat_num[pro]
            title = pro_wdsp[0].split()
            if r_num%6 == 0:
                n = r_num/6
                if n == 1:
                    for l in pro_wdsp:
                        print >> w_f,l
                else:
                    splits = [pro_wdsp[1:][i:i+6] for i in range(0,r_num,6)]
                    for i in range(n):
                        print >> w_f,' '.join([title[0],title[1]+str(i+1)]+title[2:])
                        for l in splits[i]:
                            print >> w_f,l
            elif r_num%7 == 0:
                n = r_num/7
                if n == 1:
                    for l in pro_wdsp:
                        print >> w_f,l
                else:
                    splits = [pro_wdsp[1:][i:i+7] for i in range(0,r_num,7)]
                    for i in range(n):
                        print >> w_f,' '.join([title[0],title[1]+str(i+1)]+title[2:])
                        for l in splits[i]:
                            print >> w_f,l

            if r_num%8 == 0:
                n = r_num/8
                if n == 1:
                    for l in pro_wdsp:
                        print >> w_f,l
                else:
                    splits = [pro_wdsp[1:][i:i+8] for i in range(0,r_num,8)]
                    for i in range(n):
                        print >> w_f,' '.join([title[0],title[1]+str(i+1)]+title[2:])
                        for l in splits[i]:
                            print >> w_f,l
            else:
                pass

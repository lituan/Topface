#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
select wd648s with repeats 6n,7n,8n
"""
import os
import sys
from collections import OrderedDict


class Wdsp():

    """
    use wdsp output file as input, gives you pros, pro_hotspots, pro_blades
    pro_scores etc.
    """

    def __init__(self, wdsp_f):
        """TODO: to be defined1. """
        self.pros = []
        self.seqs = OrderedDict()
        self.hotspots = OrderedDict()
        self.blades = OrderedDict()
        self.wdsps = OrderedDict()
        self.scores = OrderedDict()
        self.blade_scores = OrderedDict()
        self.repeat_num = OrderedDict()
        self.tetrad_num = OrderedDict()
        self.repeats = OrderedDict()

        self.lines = wdsp_f.readlines()
        self.lines = [line.strip('\n') for line in self.lines]
        self.lines = filter(lambda x: len(x.split()) > 0, self.lines)

        self.get_wdsps()
        self.get_blades()
        self.get_hotspots()

    def get_wdsps(self):
        for line in self.lines:
            words = line.split()
            if words[0] == '>':
                pro_name = words[1]
                # pro_name = words[1].split('|')[2]
                self.pros.append(pro_name)
                try:
                    self.scores[pro_name] = float(words[2])
                except:
                    print words
                pro_wdsp = [line]
            elif len(words) > 4:
                pro_wdsp.append(line)
                self.wdsps[pro_name] = pro_wdsp

    def get_blades(self):
        for pro, wdsp in self.wdsps.iteritems():
            self.blades[pro] = [line.split()[3:-1] for line in wdsp[1:]]
            self.blade_scores[pro] = [
                float(line.split()[-1]) for line in wdsp]
            self.seqs[pro] = ''.join([''.join(blade)
                                      for blade in self.blades[pro]])
            self.repeat_num[pro] = len(wdsp) - 1
            self.tetrad_num[pro] = len(
                [s for s in self.blade_scores[pro] if s >= 44.0])
            self.repeats[pro] = [''.join(b) for b in self.blades[pro]]

    def get_hotspots(self):
        for pro, blades in self.blades.iteritems():
            hotspot = []
            for blade in blades:
                R1 = blade[2][1]
                R1_2 = blade[1][-1]
                if len(blade[5]) <= 5 and blade[5][1] == 'D':
                    D_1 = blade[5][0]
                elif len(blade[5]) == 3 or len(blade[5]) == 2:
                    D_1 = blade[5][0]
                elif 3 <= len(blade[5]) <= 5 and blade[5][2] == 'D':
                    D_1 = blade[5][1]
                elif 4 <= len(blade[5]) <= 5 and blade[5][3] == 'D':
                    D_1 = blade[5][2]
                elif 5 <= len(blade[5]) <= 5 and blade[5][4] == 'D':
                    Di_1 = blade[5][3]
                elif len(blade[5]) <= 5:
                    D_1 = blade[5][1]
                elif len(blade[5]) <= 7:
                    D_1 = blade[5][0]
                else:
                    D_1 = '*'
                hotspot.append(R1 + R1_2 + D_1)
            self.hotspots[pro] = hotspot


def main():
    fname = os.path.split(sys.argv[-1])[1].split('.')[0]

    with open(sys.argv[-1]) as wdsp_f:
        wdsp = Wdsp(wdsp_f)

        repeat_nums = [6,7,8]
        wdsp_select = [[] for i in range(len(repeat_nums))]

        for pro,pro_wdsp in wdsp.wdsps.iteritems():
            r_num = wdsp.repeat_num[pro]
            title = pro_wdsp[0].split()
            for ri,repeat_num in enumerate(repeat_nums):
                if r_num >= repeat_num:
                    if r_num%repeat_num == 0:
                        n = r_num/repeat_num
                        if n == 1:
                            wdsp_select[ri].append([pro,pro_wdsp])
                            # wdsp_6.append([pro,pro_wdsp])
                        elif n > 1:
                            splits = [pro_wdsp[1:][i:i+repeat_num] for i in range(0,r_num,repeat_num)]
                            for i in range(n):
                                pro_i = title[1]+str(i+1)
                                wdsp_i = [' '.join([title[0],title[1]+str(i+1)]+title[2:])] + splits[i]
                                wdsp_select[ri].append([pro_i,wdsp_i])
                                # wdsp_6.append([pro_i,wdsp_i])

        for wi,selected_wdsp in enumerate(wdsp_select):
            print 'repeat_num',repeat_nums[wi],'has',len(selected_wdsp)
            with open(fname+'_'+str(repeat_nums[wi])+'.wdsp','w') as w_f:
                for pro,wdsp in selected_wdsp:
                    for line in wdsp:
                        print >> w_f,line

        total = []
        for w in wdsp_select:
            total += w
        print 'total',len(total)
        with open(fname+'_'+'_'.join(map(str,repeat_nums))+'.wdsp','w') as w_f:
            for pro,wdsp in total:
                for line in wdsp:
                    print >> w_f,line


if __name__ == "__main__":
    main()

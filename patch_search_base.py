# -*- coding:utf-8 -*-
import itertools
import operator
import os
import sys
from collections import defaultdict

import lt
from wdsp import Wdsp


class PatchSearch:

    def __init__(self, pros, seqs, wdsps, hotspots, cutoff=1):
        '''
        shape_name: (blade_num,positin_num(0:R1, 1:R1_2, 2:D_1))
        '''
        self.cutoff = cutoff
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
        self.surround_shapes = {
            'shape_3_1_1_a':	((1, 0), (1, 1), (2, 2), (3, 0), (3, 1)),
            'shape_3_1_1_b':	((0, 0), (0, 2), (1, 1), (2, 0), (2, 2)),
            'shape_3_1_2_a':	((1, 0), (1, 1), (2, 1), (3, 0), (3, 1)),
            'shape_3_1_2_b':	((0, 0), (0, 2), (1, 2), (2, 0), (2, 2)),
            'shape_3_1_3_a':	((0, 0), (0, 2), (1, 1), (1, 2), (2, 0), (3, 0), (3, 1)),
            'shape_3_1_3_b':	((0, 0), (0, 2), (1, 0), (2, 1), (2, 2), (3, 0), (3, 1)),
            'shape_3_2_1_a':	((0, 0), (0, 2), (1, 1), (1, 2), (2, 2), (3, 0), (3, 1)),
            'shape_3_2_1_b':	((0, 0), (0, 2), (1, 1), (2, 1), (2, 2), (3, 0), (3, 1)),
            'shape_3_2_2_a':	((0, 0), (0, 2), (1, 1), (1, 2), (2, 1), (3, 0), (3, 1)),
            'shape_3_2_2_b':	((0, 0), (0, 2), (1, 2), (2, 1), (2, 2), (3, 0), (3, 1)),
            'shape_3_3':	((0, 0), (0, 2), (2, 0), (2, 1)),

            'shape_4_1_1_a':	((1, 0), (1, 1), (3, 0), (3, 1)),
            'shape_4_1_1_b':	((0, 0), (0, 2), (2, 0), (2, 2)),
            'shape_4_1_2_a':	((0, 0), (0, 2), (1, 0), (2, 2), (3, 0), (3, 1)),
            'shape_4_1_2_b':	((0, 0), (0, 2), (1, 1), (2, 0), (3, 0), (3, 1)),
            'shape_4_2_1_a':	((0, 0), (0, 2), (2, 1), (2, 2), (3, 0), (3, 1)),
            'shape_4_2_1_b':	((0, 0), (0, 2), (1, 1), (1, 2), (3, 0), (3, 1)),
            'shape_4_2_2_a':	((0, 0), (0, 2), (1, 2), (2, 2), (3, 0), (3, 1)),
            'shape_4_2_2_b':	((0, 0), (0, 2), (1, 1), (2, 1), (3, 0), (3, 1)),
            'shape_4_4':	((0, 0), (0, 2), (1, 1), (2, 2), (3, 0), (3, 1)),

            'shape_5_2_1_a':	((0, 0), (0, 2), (2, 2), (3, 0), (3, 1)),
            'shape_5_2_1_b':	((0, 0), (0, 2), (1, 1), (3, 0), (3, 1)),
            'shape_5_2_2_a':	((1, 1), (2, 0), (3, 1), (3, 2)),
            'shape_5_2_2_b':	((1, 1), (1, 2), (2, 0), (3, 2)),
        }

        self.aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
                   'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X', '*']
        #residues than can interact with phos using sidechians
        self.phos_res = ['R', 'K', 'Y', 'H', 'S', 'T', 'N', 'Q', 'W']
        # self.phos_res = ['R', 'K', 'Y', 'H', 'S', 'T']

        self.pros = pros
        self.seqs = seqs
        self.wdsps = wdsps
        self.hotspots = hotspots

        self.shape_pro_patches = {}
        self.shape_patch_pros = {}


    def check_patches(self, shape_k, shape_v, patches, surround_patches):
        #do not distinguish between 'S' and 'T'
        # patches = [''.join([p.replace('S','T') for p in patch]) for patch in patches]
        def validify(patch):
            # SCORE = {'R':0.36,'S':0.23,'K':0.12,'Y':0.07,'T':0.06,'H':0.03}
            # score = sum(SCORE.get(i,0) for i in patch)
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
            return 1
        good_patches = []
        for patch,surround_patch in zip(patches,surround_patches):
            if validify(patch) and validify_surround(surround_patch):
                good_patches.append(patch)

        return good_patches


    def get_patches(self):
        for shape_k, shape_v in self.shapes.iteritems():
            surround_shape = self.surround_shapes[shape_k]
            self.shape_pro_patches[shape_k] = {}
            shape_blade_len = max([s[0] for s in shape_v]) + 2
            shape_residue_len = len(shape_v)
            for hotspot_k, hotspot_v in self.hotspots.iteritems():
                hotspot_v = [(hotspot_v + hotspot_v)[i:i + shape_blade_len]
                             for i in xrange(len(hotspot_v))]
                patches = [''.join([h[s[0]][s[1]] for s in shape_v])
                           for h in hotspot_v]
                surround_patches = [''.join([h[s[0]][s[1]] for s in surround_shape])
                           for h in hotspot_v]
                patches = self.check_patches(shape_k, shape_v, patches, surround_patches)
                if len(patches) > 0:
                    self.shape_pro_patches[shape_k][hotspot_k] = patches

    def classify_patches(self):

        for shape, pro_patches in self.shape_pro_patches.iteritems():
            self.shape_patch_pros[shape] = {}
            for pro, patches in pro_patches.iteritems():
                for patch in patches:
                    if not patch in self.shape_patch_pros[shape].keys():
                        self.shape_patch_pros[shape][patch] = [pro]
                    else:
                        if not pro in self.shape_patch_pros[shape][patch]:
                            self.shape_patch_pros[shape][patch].append(pro)


    def write_results(self):

        sta_shape = []
        for shape,patch_pros in self.shape_patch_pros.iteritems():
            pros = set([pro for patch,pros in patch_pros.iteritems() for pro in pros])
            sta_shape.append((shape,len(pros),[(patch,pros) for patch,pros in patch_pros.iteritems()]))
        sta_shape = sorted(sta_shape,key=operator.itemgetter(1),reverse=True)
        with lt.open_file(file_suffix='merged_shape') as w_f:
            for shape,num,detail in sta_shape:
                print >> w_f,'{0:<20}{1:<10}{2:<}'.format(shape,num,detail)


        pros = ' '.join([' '.join([ki for k in v.values() for ki in k])
                         for k, v in self.shape_patch_pros.iteritems()]).split()
        pros = set(pros)
        with lt.open_file(file_suffix='pros_id') as f:
            for pro in pros:
                print >> f, pro

        with lt.open_file(file_suffix='hotspot') as f:
            for pro in pros:
                hotspot = self.hotspots[pro]
                print >> f, '{0:<20}{1}'.format(
                    pro, ' '.join(hotspot))

        sta = [(len(pros), shape, patch, pros) for shape, patch_pros in self.shape_patch_pros.iteritems()
               for patch, pros in patch_pros.iteritems() if len(pros) >= self.cutoff]
        sta = sorted(sta, reverse=True)
        with lt.open_file(file_suffix='sta') as f:
            for len_pros, shape, patch, pros in sta:
                print >> f, '{0:<10}{1:<20}{2:<10}{3:<}'.format(
                    len_pros, shape, patch,','.join(pros))

        for len_pros, shape, patch, pros in sta:
            # do not output these
            break
            dir_name = str(len_pros) + '_' + shape + '_' + patch

            with lt.open_file(file_name=dir_name,file_suffix='hotspot', inner_dir=dir_name) as f:
                pro_hotspots = sorted(
                    [(len(self.hotspots[pro]), pro, self.hotspots[pro]) for pro in pros])
                for _, pro, hotspot in pro_hotspots:
                    print >> f, '{0:<20}{1}'.format(
                        pro, ' '.join(hotspot))

            with lt.open_file(file_name=dir_name,file_extension='.wdsp', inner_dir=dir_name) as f:
                for pro in pros:
                    for w in self.wdsps[pro]:
                        print >> f, w

            with lt.open_file(file_name=dir_name,file_suffix='pros_id', inner_dir=dir_name) as f:
                for pro in pros:
                    print >> f, pro

            with lt.open_file(file_name=dir_name,file_suffix='seqs',file_extension='.fa', inner_dir=dir_name) as f:
                for pro in pros:
                    print >> f, '>{0}'.format(pro)
                    seq = self.seqs[pro]
                    for s in [seq[i:i+80] for i in range(0,len(seq),80)]:
                        print >> f, s

        sta_dic = {}
        for len_pros,shape,patch,pros in sta:
            patch_m = ''.join(sorted(patch))
            if not patch_m in sta_dic.keys():
                sta_dic[patch_m] = [(len_pros,shape,patch)]
            else:
                sta_dic[patch_m].append((len_pros,shape,patch))
        sta_lis = [(k,sum([vi[0] for vi in v]),v) for k,v in sta_dic.iteritems()]
        sta_lis = sorted(sta_lis,key=operator.itemgetter(1),reverse=True)
        with lt.open_file(file_suffix='merged_patch_sta') as f:
            for patch, num, detail in sta_lis:
                # detail = ' '.join(detail)
                print >> f, '{0:<10}{1:<10}{2:<10}{3:<20}{4:<}'.format(
                     patch, num, detail[0][0],detail[0][1],detail[0][2])
                for d in detail[1:]:
                    print >> f, '{0:<20}{1:<10}{2:<20}{3:<}'.format('',d[0],d[1],d[2])
        with lt.open_file(file_suffix='merged_patch_sta_simple') as f:
            for patch, num, detail in sta_lis:
                # detail = ' '.join(detail)
                print >> f, '{0:<10}{1:<10}'.format(
                     patch, num)


@lt.run_time
def main():
    with open('test.wdsp') as wdsp_f:

        wdsp = Wdsp(wdsp_f)
        a = PatchSearch(wdsp.pros,wdsp.seqs,wdsp.wdsps,wdsp.hotspots,1)
        a.get_patches()
        a.classify_patches()
        a.write_results()


if __name__ == "__main__":
    main()

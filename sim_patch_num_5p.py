# -*- coding:utf-8 -*-
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

        self.shapes = {
            'shape_5_2_1_a':	((1, 0), (1, 1), (1, 2), (2, 0), (2, 1)),
            'shape_5_2_1_b':	((1, 0), (1, 2), (2, 0), (2, 1), (2, 2)),
            'shape_5_2_2_a':	((1, 0), (1, 2), (2, 1), (2, 2), (3, 0)),
            'shape_5_2_2_b':	((1, 0), (2, 1), (2, 2), (3, 0), (3, 1)),
         }
        self.phos_res = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
                         'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    def get_patches(self):
        for shape_k, shape_v in self.shapes.iteritems():
            self.shape_pro_patches[shape_k] = {}
            shape_blade_len = max([s[0] for s in shape_v]) + 2
            shape_residue_len = len(shape_v)
            for hotspot_k, hotspot_v in self.hotspots.iteritems():
                hotspot_v = [(hotspot_v + hotspot_v)[i:i + shape_blade_len]
                             for i in xrange(len(hotspot_v))]
                patches = [''.join([h[s[0]][s[1]] for s in shape_v])
                           for h in hotspot_v]
                if len(patches) > 0:
                    self.shape_pro_patches[shape_k][hotspot_k] = patches

    def write_results(self):
        with lt.open_file(file_suffix='sim_patch_num') as w_f:
            for shape,patch_pros in self.shape_patch_pros.iteritems():
                for patch,pros in patch_pros.iteritems():
                    print >> w_f,'{0:<10}{1:<10}{2:<10}'.format(shape,patch,len(pros))


@lt.run_time
def main():
    with open(sys.argv[1]) as wdsp_f:
        CUTOFF=0
        wdsp = Wdsp(wdsp_f)
        a = PatchSearchSpecific(wdsp.pros,wdsp.seqs,wdsp.wdsps,wdsp.hotspots,CUTOFF)
        a.get_patches()
        a.classify_patches()
        a.write_results()

if __name__ == "__main__":
    main()


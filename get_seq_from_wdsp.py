"""
output sequences from a wdsp file
usage: python get_seq_from_wdsp.py wdsp_f
"""
import sys
from wdsp import Wdsp

fname = os.path.split(sys.argv[-1])[1].split('.')[0]
with open(sys.argv[-1]) as o_f:
    wdsp = Wdsp(o_f)
    seqs = wdsp.seqs

    with open(fname+'_seq.fa', 'w') as w_f:
        for pro, seq in seqs.iteritems():
            print >> w_f, '>', pro
            for s in [seq[i:i + 80] for i in range(0, len(seq), 80)]:
                print >> w_f, s

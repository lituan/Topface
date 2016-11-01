"""
get sequences by id
first file contain id, each id per line
second file contain sequences

usage: python get_seq_by_id.py id_f seq_f
"""

import os
import sys

with open(sys.argv[-2]) as id_f:
    lines = id_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line.split()[0] for line in lines if len(line.split()) > 0]
    ids = lines

with open(sys.argv[-1]) as seq_f:
    lines = seq_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    pro_line_num = [i for i, line in enumerate(
        lines) if '>' in line] + [len(lines)]
    seqs = [lines[n:pro_line_num[i + 1]]
            for i, n in enumerate(pro_line_num[:-1])]
    seqs = [(seq[0][1:], ''.join(seq[1:])) for seq in seqs]

seq_pros = [pro for pro, seq in seqs]
found = [i for i in ids if i in seq_pros]
not_found = set(ids).difference(found)
with open('found.fa', 'w') as w_f:
    found_seq = [(pro, seq) for pro, seq in seqs if pro in found]
    for pro, seq in seqs:
        if pro in found:
            print >> w_f, '>{}'.format(pro)
            for i in [seq[i:i + 80] for i in range(0, len(seq), 80)]:
                print >> w_f, i

with open('not_found.txt', 'w') as w_f:
    for i in not_found:
        print >> w_f, i

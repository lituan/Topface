"""
get wdsp output for a list of proteins

usage: python get_wdsp_by_id.py id_f wdsp_f
"""
import sys
import os
from collections import OrderedDict


def get_wdsp(wdsp_f):
    lines = wdsp_f.readlines()
    lines = [line for line in lines if len(line.split()) > 0]
    lines = [line.strip('\n') for line in lines]
    line_mark = [i for i, line in enumerate(lines) if '>' in line]
    line_mark.append(len(lines))
    begin_end = [(line_mark[i], line_mark[i + 1])
                 for i in range(len(line_mark) - 1)]
    wdsps = [lines[i:j] for i, j in begin_end]
    wdsp_dic = OrderedDict()
    for wdsp in wdsps:
        pro = wdsp[0].split()[1]
        wdsp_dic[pro] = wdsp
    return wdsp_dic


def get_id(id_f):
    lines = id_f.readlines()
    lines = [line for line in lines if len(line.split()) > 0]
    lines = [line.strip('\n') for line in lines]
    return lines


def main():
    with open(sys.argv[-2]) as id_f:
        ids = get_id(id_f)
    with open(sys.argv[-1]) as wdsp_f:
        wdsps = get_wdsp(wdsp_f)
    getted = [wdsps[pro] for pro in ids if pro in wdsps.keys()]
    getted_id = [pro for pro in ids if pro in wdsps.keys()]
    failed_id = set(ids).difference(set(getted_id))

    with open('getted.wdsp', 'w') as w_f:
        for pro in getted_id:
            wdsp = wdsps[pro]
            for line in wdsp:
                print >> w_f, line

    with open('failed.id', 'w') as w_f:
        for pro in failed_id:
            print >> w_f, pro

main()

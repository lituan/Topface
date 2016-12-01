#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
workfolw for assess sequence identity and topface identity for a group of orthologs of WD40 proteins
1. read msa, put sequences with gaps at positions where most sequences do not have gaps (90%) into bad_list
2. get non-redundant list of proteins at specific cutoff
3. read wdsp and get topfaces and pairwise topface identity
4. plot heatmap of sequence identity and topface identity
5. plot sequence logo of topface
6. plot scatter and do t-test

notes:
when calculating identity from msa, gaps are ommited if it is present in both seqs, otherwise retained

usage: python wd-oma.py wd-oma.fasta wd-oma.wdsp

"""

import sys
import os
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as spd
from collections import OrderedDict

from matplotlib.path import Path
from matplotlib import ticker
from matplotlib.patches import PathPatch
from svgpath2mpl import parse_path
from xml.dom import minidom

from wdsp import Wdsp


def read_msa(msa_f):
    with open(msa_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines if line]
        pro_line_num = [i for i, line in enumerate(
            lines) if '>' in line] + [len(lines)]
        seqs = [lines[n:pro_line_num[i + 1]]
                for i, n in enumerate(pro_line_num[:-1])]
        seqs = [(seq[0].split()[0][1:], ''.join(seq[1:])) for seq in seqs]

        return seqs


def get_badlist(seqs,cutoff=0.90):
    seq_number = len(seqs)
    seq_len = len(seqs[0][1])

    positions = [[seqs[i][1][j]
                  for i in range(seq_number)] for j in range(seq_len)]
    pfm = [len([a for a in pos if a != '-'])*1.0/seq_number for pos in positions]
    good_pos = [i for i,pos in enumerate(pfm) if pos >= cutoff]
    bad_list = []
    for i in good_pos:
        for pro,seq in seqs:
            if seq[i] == '-':
                bad_list.append(pro)
    print set(bad_list)
    return list(set(bad_list))


def get_pim(seqs):
    def pim(seq1,seq2):
        identity = len([i for i,s in enumerate(seq1) if s == seq2[i] and s !='-'])
        lens = len([i for i,s in enumerate(seq1) if i != '-' and seq2[i] != '-'])
        return identity*1.0/lens

    scores = []
    seqlen = len(seqs)
    for i in range(seqlen):
        score_i = []
        for j in range(seqlen):
            if j < i:
                score_i.append(scores[j][i])
            elif j > i:
                score_i.append(pim(seqs[i][1],seqs[j][1]))
            else:
                score_i.append(1.0)
        scores.append(score_i)
    return scores


def remove_redundancy(seqs,scores,cutoff,good_list=[],bad_list=[]):

    def get_neighbors():
        neighbors = []
        for i,seq in enumerate(seqs):
            seq_neighbors = [seq[0]]
            for j in range(len(seqs)):
                if j != i:
                    if scores[i][j] > cutoff:
                        seq_neighbors.append(seqs[j][0])
            neighbors.append(seq_neighbors)
        return neighbors


    def sort_good_bad(neighbor_list):
        new_good_list = [g for g in good_list if not g in bad_list]
        new_bad_list = [g for g in good_list if g in bad_list] + bad_list
        good = [n for n in neighbor_list if n in new_good_list]
        bad = [n for n in neighbor_list if n in new_bad_list]
        middle = [n for n in neighbor_list if (not n in new_good_list) and (not n in new_bad_list)]
        return good+middle+bad

    neighbors = get_neighbors()

    while max([len(n) for n in neighbors]) > 1:
        neighbors = sorted(neighbors,key=lambda x: len(x),reverse=True)
        neighbors = map(sort_good_bad,neighbors)
        for seq in neighbors[0][1:]:
            if not seq in good_list:
                for n in neighbors:
                    if seq in n:
                        n.pop(n.index(seq))
                neighbors = [n for n in neighbors if len(n) > 0]
        if len(neighbors[0]) > 1:
            for seq in neighbors[0][:-1]:
                for n in neighbors:
                    if seq in n:
                        n.pop(n.index(seq))
                neighbors = [n for n in neighbors if len(n) > 0]


    nr_list = set([n[0] for n in neighbors])
    nr_seqs = [seq for seq in seqs if seq[0] in nr_list ]
    nr_index = [i for i,seq in enumerate(seqs) if seq[0] in nr_list]
    nr_scores = [[score_i[i] for i in nr_index] for score_i in scores]
    nr_scores = [nr_scores[i] for i in nr_index]
    return nr_seqs,nr_scores

def get_wdsp(wdsp_f,nr_list):

    with open(wdsp_f) as w_f:
        wdsp_lines = w_f.readlines()
        wdsp_lines = [line.rstrip('\r\n') for line in wdsp_lines if line]
        wdsp_hotspot = OrderedDict()
        blade = OrderedDict()
        for line in wdsp_lines:
            words = line.split()
            if len(words) >= 2 and words[0] == '>':
                pro_name = words[1]
                pro_blades = []
                pro_seq = ''
            elif len(words) > 4:
                pro_blades.append(words[3:-1])
                pro_seq += ''.join(words[3:-1])
                if pro_name in nr_list:
                    blade[pro_name] = pro_blades
        for pro_name,pro_blades in blade.iteritems():
            hotspot = []
            for blade in pro_blades:
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
                hotspot.append(R1+R1_2+D_1)
            wdsp_hotspot[pro_name] = hotspot

    def get_hots_pim(hots):

        def pim(hot1,hot2):
            hot1 = ''.join(hot1)
            hot2 = ''.join(hot2)
            identity = len([i for i,s in enumerate(hot1) if s == hot2[i]])
            return identity*1.0/len(hot1)

        scores = []
        hotlen = len(hots)
        for i in range(hotlen):
            score_i = []
            for j in range(hotlen):
                if j < i:
                    score_i.append(scores[j][i])
                elif j > i:
                    score_i.append(pim(hots[i][1],hots[j][1]))
                else:
                    score_i.append(1.0)
            scores.append(score_i)
        return scores

    hots = [(k,wdsp_hotspot[k]) for k in nr_list]
    scores = get_hots_pim(hots)
    return hots,scores

def plot_heatmap(labels,scores,file_name,method='complete'):

    if not len(scores) > 4:
        print 'plot_heatmap: number of seqs is too small'
        return

    column_labels = labels
    row_labels = column_labels
    scores = [map(lambda x: float(x), row) for row in scores]
    scores = np.array(scores)
    distances = [map(lambda x: 1.0-x,row) for row in scores]
    linkage = sch.linkage(spd.squareform(distances),method=method)
    df = pd.DataFrame(scores,columns=column_labels, index=row_labels)

    if len(df.columns) > 30:
        sns_plot = sns.clustermap(df,row_linkage=linkage,col_linkage=linkage,xticklabels='',yticklabels='')
    else:
        # sns_plot = sns.clustermap(df,row_linkage=linkage,col_linkage=linkage,annot=True,fmt='3.2f')
        sns_plot = sns.clustermap(df,row_linkage=linkage,col_linkage=linkage)
        plt.setp(sns_plot.ax_heatmap.yaxis.get_majorticklabels(), rotation=20)
        plt.setp(sns_plot.ax_heatmap.xaxis.get_majorticklabels(), rotation=70)
    # plt.yticks(rotation=90)
    sns_plot.savefig(file_name+'.png',dpi=300)
    plt.close('all')

def plotlogo(seqs,filename):
    # chemical colors
    COLORS = {'A': 'black', 'V': 'black', 'L': 'black', 'I': 'black', 'P': 'black', 'W': 'black', 'F': 'black', 'M': 'black', 'G': 'green', 'S': 'green',
              'T': 'green', 'Y': 'green', 'C': 'green', 'Q': 'purple', 'N': 'purple', 'K': 'blue', 'R': 'blue', 'H': 'blue', 'D': 'red', 'E': 'red'}
    AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
          'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    # default font is Roboto Mono from Google
    GLYPHS={}
    GLYPHS["A"]="  M 869,-377  L 383,-377  L 266,0  L 81,0  L 551,-1456  L 706,-1456  L 1168,0  L 984,0  L 869,-377  Z   M 433,-538  L 820,-538  L 628,-1170  L 433,-538  Z    "
    GLYPHS["C"]="  M 1117,-438  Q 1104,-337 1065,-252  Q 1026,-168 963,-107  Q 901,-47 815,-13  Q 730,20 625,20  Q 535,20 461,-5  Q 388,-31 330,-76  Q 273,-121 231,-182  Q 190,-243 162,-315  Q 135,-387 121,-466  Q 108,-545 107,-626  L 107,-829  Q 108,-910 121,-989  Q 135,-1068 162,-1140  Q 190,-1212 231,-1273  Q 273,-1335 330,-1380  Q 388,-1425 461,-1450  Q 534,-1476 625,-1476  Q 734,-1476 820,-1442  Q 906,-1409 968,-1347  Q 1030,-1286 1067,-1200  Q 1105,-1114 1117,-1010  L 932,-1010  Q 923,-1076 901,-1133  Q 879,-1190 842,-1233  Q 805,-1276 751,-1300  Q 698,-1325 625,-1325  Q 559,-1325 508,-1303  Q 458,-1282 421,-1245  Q 385,-1208 360,-1159  Q 336,-1110 321,-1055  Q 306,-1000 299,-942  Q 293,-885 293,-831  L 293,-626  Q 293,-572 299,-514  Q 306,-457 321,-401  Q 336,-346 360,-297  Q 384,-248 421,-210  Q 458,-173 508,-151  Q 558,-130 625,-130  Q 698,-130 751,-153  Q 805,-176 842,-217  Q 879,-259 901,-315  Q 923,-372 932,-438  L 1117,-438  Z    "
    GLYPHS["B"]="  M 172,0  L 172,-1456  L 605,-1456  Q 695,-1455 781,-1432  Q 867,-1410 934,-1363  Q 1001,-1316 1041,-1243  Q 1081,-1170 1080,-1068  Q 1079,-1011 1061,-964  Q 1044,-917 1013,-879  Q 983,-842 942,-814  Q 901,-786 855,-766  Q 913,-749 961,-717  Q 1010,-686 1045,-641  Q 1080,-597 1099,-541  Q 1119,-485 1119,-420  Q 1120,-318 1080,-240  Q 1040,-162 973,-109  Q 906,-57 818,-29  Q 731,-1 638,0  L 172,0  Z   M 358,-681  L 358,-157  L 643,-157  Q 701,-158 753,-176  Q 805,-195 845,-228  Q 885,-261 908,-309  Q 932,-357 932,-418  Q 933,-480 911,-528  Q 890,-576 852,-609  Q 814,-643 763,-661  Q 712,-679 653,-681  L 358,-681  Z   M 358,-835  L 616,-835  Q 668,-836 718,-851  Q 768,-867 807,-896  Q 846,-926 870,-969  Q 894,-1012 894,-1069  Q 894,-1130 870,-1173  Q 847,-1216 807,-1243  Q 768,-1271 716,-1284  Q 665,-1297 611,-1298  L 358,-1298  L 358,-835  Z    "
    GLYPHS["E"]="  M 975,-673  L 367,-673  L 367,-157  L 1076,-157  L 1076,0  L 182,0  L 182,-1456  L 1067,-1456  L 1067,-1298  L 367,-1298  L 367,-830  L 975,-830  L 975,-673  Z    "
    GLYPHS["D"]="  M 155,0  L 155,-1456  L 492,-1456  Q 644,-1454 763,-1404  Q 883,-1355 965,-1266  Q 1048,-1178 1091,-1054  Q 1135,-931 1136,-781  L 1136,-674  Q 1135,-524 1091,-400  Q 1048,-277 965,-188  Q 883,-100 763,-50  Q 644,-1 492,0  L 155,0  Z   M 343,-1304  L 343,-151  L 492,-151  Q 610,-152 696,-192  Q 782,-233 838,-303  Q 895,-373 922,-468  Q 950,-563 951,-674  L 951,-783  Q 950,-894 922,-988  Q 894,-1083 838,-1152  Q 782,-1222 696,-1262  Q 610,-1302 492,-1304  L 343,-1304  Z    "
    GLYPHS["G"]="  M 1116,-191  Q 1024,-83 905,-31  Q 786,21 644,20  Q 554,19 478,-7  Q 403,-34 343,-80  Q 283,-127 238,-190  Q 193,-253 162,-326  Q 132,-400 116,-480  Q 101,-561 100,-643  L 100,-812  Q 101,-893 114,-973  Q 128,-1054 156,-1128  Q 184,-1202 226,-1265  Q 268,-1329 326,-1375  Q 385,-1422 459,-1449  Q 533,-1476 625,-1476  Q 727,-1476 813,-1444  Q 899,-1413 963,-1354  Q 1027,-1296 1066,-1213  Q 1105,-1131 1114,-1029  L 931,-1029  Q 920,-1092 897,-1145  Q 874,-1199 836,-1237  Q 799,-1275 747,-1296  Q 695,-1318 626,-1318  Q 560,-1318 509,-1295  Q 458,-1273 420,-1235  Q 382,-1197 357,-1146  Q 332,-1096 316,-1040  Q 300,-984 293,-926  Q 286,-868 286,-814  L 286,-643  Q 287,-588 295,-529  Q 304,-471 322,-415  Q 340,-359 368,-309  Q 396,-259 435,-221  Q 475,-183 527,-160  Q 579,-138 645,-137  Q 683,-136 724,-140  Q 766,-144 805,-156  Q 844,-168 878,-188  Q 912,-209 935,-242  L 937,-569  L 641,-569  L 641,-725  L 1113,-725  L 1116,-191  Z    "
    GLYPHS["F"]="  M 984,-643  L 378,-643  L 378,0  L 191,0  L 191,-1456  L 1085,-1456  L 1085,-1298  L 378,-1298  L 378,-800  L 984,-800  L 984,-643  Z    "
    GLYPHS["I"]="  M 174,-1456  L 1054,-1456  L 1054,-1295  L 705,-1295  L 705,-160  L 1054,-160  L 1054,0  L 174,0  L 174,-160  L 515,-160  L 515,-1295  L 174,-1295  L 174,-1456  Z    "
    GLYPHS["H"]="  M 1087,0  L 912,0  L 912,-673  L 315,-673  L 315,0  L 141,0  L 141,-1456  L 315,-1456  L 315,-830  L 912,-830  L 912,-1456  L 1087,-1456  L 1087,0  Z    "
    GLYPHS["K"]="  M 523,-676  L 361,-492  L 361,0  L 172,0  L 172,-1456  L 361,-1456  L 361,-745  L 502,-921  L 929,-1456  L 1154,-1456  L 645,-819  L 1188,0  L 963,0  L 523,-676  Z    "
    GLYPHS["J"]="  M 857,-1456  L 1046,-1456  L 1046,-443  Q 1044,-342 1007,-257  Q 971,-172 908,-110  Q 845,-49 759,-14  Q 674,20 573,20  Q 471,20 387,-11  Q 304,-42 242,-99  Q 181,-157 144,-238  Q 107,-320 98,-421  L 286,-421  Q 289,-360 310,-308  Q 332,-256 369,-217  Q 406,-179 457,-158  Q 509,-137 573,-137  Q 639,-137 691,-161  Q 743,-186 779,-228  Q 816,-270 835,-325  Q 855,-381 857,-443  L 857,-1456  Z    "
    GLYPHS["M"]="  M 377,-1456  L 614,-728  L 870,-1456  L 1100,-1456  L 1100,0  L 920,0  L 920,-581  L 935,-1189  L 666,-405  L 560,-405  L 313,-1168  L 328,-581  L 328,0  L 148,0  L 148,-1456  L 377,-1456  Z    "
    GLYPHS["L"]="  M 383,-157  L 1095,-157  L 1095,0  L 198,0  L 198,-1456  L 383,-1456  L 383,-157  Z    "
    GLYPHS["O"]="  M 1121,-644  Q 1120,-566 1107,-486  Q 1095,-407 1069,-333  Q 1043,-259 1002,-195  Q 962,-131 906,-83  Q 850,-35 777,-7  Q 705,20 615,20  Q 525,20 452,-7  Q 380,-35 324,-83  Q 268,-131 227,-195  Q 186,-260 159,-334  Q 133,-408 120,-487  Q 107,-566 106,-644  L 106,-810  Q 107,-888 119,-967  Q 132,-1047 158,-1121  Q 185,-1195 225,-1259  Q 266,-1324 322,-1372  Q 378,-1421 450,-1448  Q 523,-1476 613,-1476  Q 703,-1476 776,-1448  Q 849,-1421 905,-1373  Q 961,-1325 1001,-1260  Q 1042,-1196 1068,-1122  Q 1095,-1048 1107,-968  Q 1120,-888 1121,-810  L 1121,-644  Z   M 938,-812  Q 937,-864 931,-920  Q 925,-977 910,-1032  Q 896,-1088 872,-1138  Q 848,-1189 812,-1227  Q 776,-1266 727,-1288  Q 678,-1311 613,-1311  Q 549,-1311 500,-1288  Q 451,-1265 415,-1226  Q 379,-1188 355,-1137  Q 331,-1087 316,-1031  Q 302,-976 295,-919  Q 289,-863 288,-812  L 288,-644  Q 289,-593 295,-536  Q 302,-479 317,-423  Q 332,-368 356,-317  Q 380,-266 416,-227  Q 452,-189 501,-166  Q 550,-143 615,-143  Q 680,-143 729,-166  Q 779,-189 814,-227  Q 850,-266 873,-316  Q 897,-367 911,-422  Q 926,-478 931,-535  Q 937,-592 938,-644  L 938,-812  Z    "
    GLYPHS["N"]="  M 1086,0  L 898,0  L 333,-1088  L 330,0  L 143,0  L 143,-1456  L 331,-1456  L 896,-370  L 899,-1456  L 1086,-1456  L 1086,0  Z    "
    GLYPHS["Q"]="  M 1134,-663  Q 1133,-582 1120,-500  Q 1107,-418 1080,-342  Q 1054,-266 1012,-199  Q 971,-133 913,-84  L 1164,125  L 1037,246  L 749,2  Q 687,20 615,20  Q 522,20 448,-8  Q 374,-37 316,-86  Q 259,-135 217,-201  Q 175,-268 148,-344  Q 121,-420 108,-501  Q 95,-583 94,-663  L 94,-791  Q 95,-871 108,-953  Q 121,-1035 148,-1111  Q 175,-1187 216,-1253  Q 257,-1320 314,-1369  Q 372,-1419 446,-1447  Q 521,-1476 614,-1476  Q 707,-1476 781,-1447  Q 856,-1419 913,-1370  Q 971,-1321 1012,-1254  Q 1054,-1188 1080,-1112  Q 1107,-1036 1120,-953  Q 1133,-871 1134,-791  L 1134,-663  Z   M 950,-793  Q 949,-848 943,-907  Q 938,-967 923,-1024  Q 909,-1082 884,-1134  Q 860,-1186 823,-1225  Q 786,-1265 734,-1288  Q 683,-1311 614,-1311  Q 546,-1311 495,-1287  Q 444,-1264 407,-1224  Q 370,-1185 345,-1133  Q 321,-1081 306,-1023  Q 291,-966 284,-906  Q 278,-847 278,-793  L 278,-663  Q 278,-609 284,-549  Q 291,-490 306,-432  Q 321,-374 345,-322  Q 370,-270 407,-230  Q 444,-190 495,-166  Q 547,-143 615,-143  Q 684,-143 735,-166  Q 787,-190 824,-229  Q 861,-269 885,-321  Q 909,-373 923,-431  Q 938,-489 943,-548  Q 949,-608 950,-663  L 950,-793  Z    "
    GLYPHS["P"]="  M 376,-584  L 376,0  L 191,0  L 191,-1456  L 663,-1456  Q 761,-1454 848,-1425  Q 936,-1396 1002,-1341  Q 1068,-1286 1106,-1205  Q 1145,-1124 1145,-1019  Q 1145,-914 1106,-833  Q 1068,-753 1002,-698  Q 936,-643 848,-614  Q 761,-585 663,-584  L 376,-584  Z   M 376,-736  L 663,-736  Q 727,-737 781,-756  Q 835,-776 875,-812  Q 915,-848 937,-899  Q 960,-951 960,-1017  Q 960,-1083 937,-1136  Q 915,-1189 875,-1226  Q 836,-1263 781,-1283  Q 727,-1303 663,-1304  L 376,-1304  L 376,-736  Z    "
    GLYPHS["S"]="  M 936,-368  Q 936,-435 905,-481  Q 875,-527 827,-558  Q 779,-590 721,-611  Q 664,-632 611,-649  Q 534,-674 454,-709  Q 375,-744 309,-795  Q 244,-846 202,-915  Q 161,-985 161,-1079  Q 161,-1173 202,-1247  Q 244,-1321 311,-1372  Q 378,-1423 463,-1449  Q 548,-1476 634,-1476  Q 729,-1476 817,-1444  Q 905,-1413 973,-1356  Q 1041,-1299 1082,-1218  Q 1123,-1137 1125,-1037  L 935,-1037  Q 927,-1100 904,-1151  Q 881,-1203 843,-1240  Q 805,-1277 752,-1297  Q 700,-1318 634,-1318  Q 581,-1318 530,-1303  Q 480,-1288 440,-1258  Q 401,-1228 377,-1184  Q 354,-1140 354,-1082  Q 355,-1019 386,-975  Q 417,-932 464,-902  Q 512,-872 567,-852  Q 623,-832 672,-817  Q 726,-800 781,-778  Q 836,-757 887,-729  Q 938,-701 982,-666  Q 1026,-631 1059,-586  Q 1092,-542 1110,-488  Q 1129,-435 1129,-370  Q 1129,-272 1085,-199  Q 1042,-126 973,-77  Q 904,-29 817,-4  Q 730,20 643,20  Q 546,20 453,-10  Q 360,-40 286,-96  Q 213,-153 167,-234  Q 121,-316 118,-420  L 307,-420  Q 316,-352 344,-299  Q 372,-247 416,-210  Q 460,-174 517,-155  Q 575,-137 643,-137  Q 697,-137 749,-150  Q 802,-164 843,-192  Q 884,-221 910,-264  Q 936,-308 936,-368  Z    "
    GLYPHS["R"]="  M 656,-594  L 365,-594  L 365,0  L 181,0  L 181,-1456  L 608,-1456  Q 710,-1454 800,-1427  Q 890,-1400 957,-1346  Q 1025,-1292 1063,-1210  Q 1102,-1129 1102,-1019  Q 1102,-948 1081,-889  Q 1061,-830 1025,-782  Q 989,-734 939,-697  Q 889,-660 829,-634  L 1138,-12  L 1137,0  L 942,0  L 656,-594  Z   M 365,-746  L 613,-746  Q 675,-747 730,-765  Q 785,-784 826,-819  Q 868,-854 892,-904  Q 916,-955 916,-1021  Q 916,-1091 893,-1143  Q 870,-1196 829,-1231  Q 788,-1267 731,-1285  Q 675,-1303 608,-1304  L 365,-1304  L 365,-746  Z    "
    GLYPHS["U"]="  M 1088,-1456  L 1090,-470  Q 1088,-368 1053,-279  Q 1018,-190 955,-123  Q 892,-57 805,-18  Q 718,20 614,20  Q 508,20 421,-18  Q 334,-56 272,-122  Q 210,-189 175,-278  Q 140,-368 139,-470  L 141,-1456  L 317,-1456  L 321,-470  Q 322,-405 341,-345  Q 361,-285 397,-239  Q 434,-193 488,-165  Q 543,-137 614,-137  Q 685,-137 739,-164  Q 793,-192 829,-238  Q 866,-285 885,-345  Q 904,-405 906,-470  L 909,-1456  L 1088,-1456  Z    "
    GLYPHS["T"]="  M 1156,-1298  L 706,-1298  L 706,0  L 526,0  L 526,-1298  L 76,-1298  L 76,-1456  L 1156,-1456  L 1156,-1298  Z    "
    GLYPHS["W"]="  M 896,-394  L 1007,-1456  L 1182,-1456  L 1005,0  L 816,0  L 629,-1097  L 440,0  L 250,0  L 73,-1456  L 249,-1456  L 360,-394  L 547,-1456  L 708,-1456  L 896,-394  Z    "
    GLYPHS["V"]="  M 610,-298  L 954,-1456  L 1151,-1456  L 692,0  L 531,0  L 71,-1456  L 269,-1456  L 610,-298  Z    "
    GLYPHS["Y"]="  M 603,-725  L 935,-1456  L 1145,-1456  L 692,-543  L 689,0  L 517,0  L 514,-543  L 61,-1456  L 272,-1456  L 603,-725  Z    "
    GLYPHS["X"]="  M 625,-885  L 939,-1456  L 1157,-1456  L 734,-734  L 1167,0  L 951,0  L 629,-582  L 306,0  L 87,0  L 521,-734  L 98,-1456  L 315,-1456  L 625,-885  Z    "
    GLYPHS["Z"]="  M 325,-157  L 1079,-157  L 1079,0  L 116,0  L 114,-144  L 839,-1298  L 127,-1298  L 127,-1456  L 1050,-1456  L 1052,-1315  L 325,-157  Z    "

    def get_glyph_patch(path_data, color, x, y, dx, dy, **kwargs):
        kwargs.setdefault('facecolor', color)
        kwargs.setdefault('edgecolor', 'none')
        path = parse_path(path_data)
        # normalize and flip upside down
        path.vertices[:, 0] -= path.vertices[:, 0].min()
        path.vertices[:, 1] -= path.vertices[:, 1].min()
        path.vertices[:, 0] /= path.vertices[:, 0].max()
        path.vertices[:, 1] /= path.vertices[:, 1].max()
        path.vertices[:, 1] = 1 - path.vertices[:, 1]
        # scale then translate
        path.vertices *= [dx, dy]
        path.vertices += [x, y]
        return PathPatch(path, **kwargs)


    def draw_logo(ax, matrix, charwidth):
        # glyphs = read_glyph()
        for i, (_, position) in enumerate(matrix.iterrows()):
            letters_sorted = position.sort_values()
            bottom = 0
            for letter, height in letters_sorted.iteritems():
                patch = get_glyph_patch(GLYPHS[letter], COLORS[
                                        letter], i * charwidth, bottom, charwidth, height)
                ax.add_artist(patch)
                bottom += height


    def plot_seqlogo(ax, pfm, charwidth=1.0, xlim='',**kwargs):

        info_content = np.log2(
            20) - pfm.apply(lambda p: (-p * np.log2(p)).sum(), axis=1)
        matrix = pfm.mul(info_content, axis=0)

        seqlen = len(pfm)
        draw_logo(ax, matrix, charwidth, **kwargs)
        # xlim
        if xlim:
            ax.set_xlim([0, xlim])
        else:
            ax.set_xlim([0, seqlen * charwidth])
        # major ticks
        ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0, seqlen*charwidth,charwidth)))
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.tick_params(which='major', direction='out')
        # minor ticks
        ax.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(0, seqlen*charwidth,charwidth) + charwidth*0.5))
        ax.xaxis.set_minor_formatter(
            ticker.FixedFormatter(np.arange(1, seqlen + 1)))
        ax.tick_params(which='minor', length=0)
        # ylim
        ax.set_ylim([0, 4.5])
        ax.yaxis.set_major_locator(ticker.FixedLocator([0., 1., 2., 3., 4.]))
        # show axis an left and bottom
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        # hide top and right spines
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        # set x y label
        ax.set_xlabel('positions')
        ax.set_ylabel('bits')


    seq_number = len(seqs)
    if seq_number < 4:
        print 'plot_seqlogo: number of seqs is too small'
        return
    seq_len = len(seqs[0][1])
    positions = [[seqs[i][1][j]
                  for i in range(seq_number)] for j in range(seq_len)]
    pfm = [[pos.count(a) * 1.0 / seq_number for a in AA] for pos in positions]
    pfm = [pos for pos in pfm if sum(pos) >= 0.95]
    pfm = pd.DataFrame(pfm, columns=AA)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(1)
    plot_seqlogo(ax, pfm)
    plt.savefig(filename+'.png',dpi=300)
    plt.close('all')

def plot_scatter(seq_identity,topface_identity,filename):
    p = {'seq_identity':seq_identity,'topface_identity':topface_identity}
    p = pd.DataFrame(p)
    g = sns.lmplot('seq_identity','topface_identity',p)
    g.set(xlim=(0,1.0),ylim=(0,1.0))
    test = scipy.stats.ttest_ind(topface_identity,seq_identity,equal_var=False)
    x = [0,0.2,0.4,0.6,0.8,1.0]
    y = [0,0.2,0.4,0.6,0.8,1.0]
    g.ax.plot(x,y,c='black',linestyle='--',lw=1,alpha=0.5)

    plt.title('statistic:'+'{0:4.3f}'.format(test.statistic)+'    p-value:'+'{0:<8.7f}'.format(test.pvalue*0.5))
    plt.savefig(filename+'.png',dpi=300)

def matrix_to_pair(scores):
    pairs = []
    len_s = len(scores)
    for i in range(len_s):
        for j in range(len_s):
            if j < len_s and j > i:
                pairs.append(scores[i][j])
    return pairs

def write_result(nr_seqs,filename):
    with open(filename+'.fasta','w') as w_f:
        for pro,seq in nr_seqs:
            print >> w_f, '>'+pro
            for s in [seq[i:i+80] for i in range(0,len(seq),80)]:
                print >> w_f,s

def main():
    seqs = read_msa(sys.argv[-2])
    filename = os.path.splitext(os.path.split(sys.argv[-2])[1])[0]
    cutoff = 0.9
    good_list=[]
    scores = get_pim(seqs)
    bad_list = get_badlist(seqs)
    nr_seqs,nr_seq_scores = remove_redundancy(seqs,scores,cutoff,good_list,bad_list)
    nr_list = [seq[0] for seq in nr_seqs]
    nr_hots,nr_hots_scores = get_wdsp(sys.argv[-1],nr_list)
    plot_heatmap(nr_list,nr_seq_scores,filename+'_seq'+'_'+str(cutoff))
    plot_heatmap(nr_list,nr_hots_scores,filename+'_topface'+'_'+str(cutoff))
    plot_scatter(matrix_to_pair(nr_seq_scores),matrix_to_pair(nr_hots_scores),filename+'_seq_topface_identity'+'_'+str(cutoff))
    nr_hots = [[k,''.join(v)] for k,v in nr_hots]
    plotlogo(nr_hots,filename+'_topface_logo'+'_'+str(cutoff))
    write_result(nr_seqs,filename+'_nr'+'_'+str(cutoff))

main()





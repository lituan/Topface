#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
use UniProt REST service to search sequences annotated as 'WD40'
"""

import os
import sys
import urllib
import urllib2
import cPickle as pickle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from multiprocessing import Pool

def uniprot_wd40(key='pfam',pdb=False):
    """
    use annotations from different database to query WD40 in uniprot
    return a list of uniprot accesions
    """
    if   key == 'pfam':
        query = 'database:(type:pfam id:PF00400) or database:(type:pfam id:PF12894) or database:(type:pfam id:PF16529) or database:(type:pfam id:PF16756) or database:(type:pfam id:PF17005)'
    elif key == 'smart':
        query = 'database:(type:smart id:SM00320)'
    elif key == 'supfam':
        query = 'database:(type:supfam id:SSF50978)'
    elif key == 'prosite':
        query = 'database:(type:prosite id:PS00678) or database:(type:prosite id:PS50082) or database:(type:prosite id:PS50294)'
    elif key == 'uniprot':
        query = 'keyword:"WD repeat" or annotation:(type:repeat wd) or family:"WD repeat"'
    else:
        print 'wrong query key'
        return

    if pdb:
        query = query + ' AND '+ 'database:(type:pdb)'

    url = ' http://www.uniprot.org/uniprot/?'
    data ={
    'query':query,
    'format':'list',
    }
    data = urllib.urlencode(data)
    req = urllib2.Request(url,data)
    response = urllib2.urlopen(req)
    r = response.readlines()
    lines = set([line.rstrip('\r\n') for line in r])

    print 'uniprot search',key,'is finished'
    return key,lines


def write_lis_lis(lis_lis,filename,cols=[]):
    """align nested list to print a table"""
    lis_lis = [[str(l) for l in lis]
               for lis in lis_lis]  # trans every element to str
    #make all inner lists of the same length
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [lis + (inner_lis_max_len - len(lis)) * [''] for lis in lis_lis]
    #make element in the same list have the same length
    aligned = []
    for lis in lis_lis:
        width = max([len(l) for l in lis])
        lis = [l + (width - len(l)) * ' ' for l in lis]
        aligned.append(lis)
    new_lis_lis = [';'.join([aligned[i][j] for i in range(len(aligned))]) for j in range(len(aligned[0]))]
    with open(filename+'.txt','w') as w_f:
        if cols:
            print >> w_f,'\t;'.join(cols)
        for l in new_lis_lis:
            print >> w_f,l

def get_wdsp_acc():
    with open('wdsp_uniprot_id.txt') as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        lines = [line.split()[0] for line in lines]
        return lines

def main():
    keywords = ['pfam','smart','supfam','prosite','uniprot']
    p = Pool(5)
    result = p.map(uniprot_wd40,keywords)
    p.close()
    wd40s = []
    for k in keywords:
        for r,v in result:
            if r == k:
                wd40s.append(v)
    # wd40s = [v for k in keywords for r,v in result if r == k]
    wdsp = get_wdsp_acc()
    wd40s.append(wdsp)
    keywords.append('wdsp')

    pickle.dump(wd40s,open('wd40s.pickle','w'))
    wd40s = pickle.load(open('wd40s.pickle'))

    # barplot of wd40s annotated by different database or method
    total = []
    for w in wd40s:
        total += w
    total = list(set(total))
    sns.set_color_codes('pastel')
    f,ax = plt.subplots(figsize=(6,5))
    keys = ['Pfam','SMART','Superfamily','Prosite','UniProt','WDSP']
    wd = pd.DataFrame({'Methods':keys,'Sequence Num':map(len,[total for i in range(6)])})
    sns.barplot(x='Methods',y='Sequence Num',data=wd,color='b')
    sns.set_color_codes('muted')
    wd = pd.DataFrame({'Methods':keys,'Sequence Num':map(len,wd40s)})
    wd = wd.sort_values('Sequence Num',ascending=False)
    h = sns.barplot(x='Methods',y='Sequence Num',data=wd,color='b')
    h.figure.subplots_adjust(top=0.9,bottom=0.10,left=0.15,right=0.95)
    ax.set(xlabel='Methods',ylabel='Sequence Num',title='WD40 Proteins Annotated by Different Methods')
    plt.savefig('wd40_annotated_by_different_methods_acc',dpi=300)
    plt.close('all')
    write_lis_lis(wd40s,'wd40_annotated_by_different_databases_accs',keys)

    # percent table (heatmap) of wd40s shared by two different database or method
    table = []
    keys = ['Pfam','SMART','Superfamily','Prosite','UniProt','WDSP']
    for w in wd40s:
        row = [len(set(w).intersection(set(wr)))*1.0/len(w) for wr in wd40s]
        table.append(row)
    data = pd.DataFrame(table,columns=keys,index=keys)
    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111)
    h = sns.heatmap(data,annot=True,fmt='.2f',cmap='Blues')
    h.figure.subplots_adjust(top=0.87,bottom=0.07,left=0.16,right=0.96)
    ax.set_xticklabels(keys,rotation=0)
    ax.set_yticklabels(keys[::-1],rotation=0)
    ax.set_title('Comaration of Different Annotation Methods')
    plt.savefig('Comaration_of Different_Annotation_Methods.png',dpi=300)
    plt.close('all')

    total = set.union(*map(set,wd40s))
    wd40s_score = [[] for i in range(6)]
    def acc_score(acc):
        i = 0
        for w in wd40s:
            if acc in w:
                i += 1
        return i
    for acc in total:
        score = acc_score(acc)
        wd40s_score[score-1].append(acc)

    pickle.dump(wd40s_score,open('wd40s_score.pickle','w'))
    wd40s_score = pickle.load(open('wd40s_score.pickle'))
    sns.set_color_codes('pastel')
    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111)
    wd = pd.DataFrame({'Methods':keys,'Sequence Num':map(len,[total for i in range(6)])})
    sns.barplot(x='Methods',y='Sequence Num',data=wd,color='b')
    sns.set_color_codes('muted')
    wd = pd.DataFrame({'Methods Score':range(1,7),'Sequence Num':map(len,wd40s_score)})
    h = sns.barplot(x='Methods Score',y='Sequence Num',data=wd,color='b')
    h.figure.subplots_adjust(top=0.87,bottom=0.11,left=0.16,right=0.96)
    ax.set(xlabel='Annotation Score',ylabel='Sequence Num',title='Annotation Score of WD40 Proteins in UniProt')
    plt.savefig('wd40_annotation_score_in_uniprot_accs',dpi=300)
    plt.close('all')
    write_lis_lis(wd40s_score,'wd40_annotation_score_in_uniprot_accs',[str(i) for i in range(1,7)])


if __name__ == "__main__":
    main()

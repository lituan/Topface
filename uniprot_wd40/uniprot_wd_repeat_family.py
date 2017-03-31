#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
query wd40 using uniprot annotation
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
from matplotlib_venn import venn2,venn3

def uniprot_wd40(key,pdb=False):

    if key == 'uniprot_keyword':
        query = 'keyword:"WD repeat"'
    elif key == 'uniprot_repeat':
        query = 'annotation:(type:repeat wd)'
    elif key == 'uniprot_family':
        query = 'family:"WD repeat"'
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


def main():
    keywords = ['uniprot_keyword','uniprot_repeat','uniprot_family']
    p = Pool(3)
    result = p.map(uniprot_wd40,keywords)
    p.close()
    uniprot_wd40s = []
    for k in keywords:
        for r,v in result:
            if r == k:
                uniprot_wd40s.append(v)
    uniprot_wd40s = [v for k in keywords for r,v in result if r == k]

    pickle.dump(uniprot_wd40s,open('uniprot_wd40s.pickle','w'))
    uniprot_wd40s = pickle.load(open('uniprot_wd40s.pickle'))

    # barplot of uniprot_wd40s annotated by different database or method
    total = set.union(*map(set,uniprot_wd40s))
    sns.set_color_codes('pastel')
    f,ax = plt.subplots(figsize=(6,5))
    keys = ['UniProt_Keyword','UniProt_Repeat','UniProt_Family']
    wd = pd.DataFrame({'Methods':keys,'Sequence Num':map(len,[total for i in range(3)])})
    sns.barplot(x='Methods',y='Sequence Num',data=wd,color='b')
    sns.set_color_codes('muted')
    wd = pd.DataFrame({'Methods':keys,'Sequence Num':map(len,uniprot_wd40s)})
    wd = wd.sort_values('Sequence Num',ascending=False)
    h = sns.barplot(x='Methods',y='Sequence Num',data=wd,color='b')
    h.figure.subplots_adjust(top=0.9,bottom=0.10,left=0.15,right=0.95)
    ax.set(xlabel='Methods',ylabel='Sequence Num',title='WD40 Proteins Annotated by Different Methods')
    plt.savefig('wd40_annotated_by_different_methods_acc',dpi=300)
    plt.close('all')

    # percent table (heatmap) of uniprot_wd40s shared by two different database or method
    table = []
    keys = ['UniProt_Keyword','UniProt_Repeat','UniProt_Family']
    for w in uniprot_wd40s:
        row = [len(set(w).intersection(set(wr)))*1.0/len(w) for wr in uniprot_wd40s]
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


    sns.set_color_codes('bright')
    set1 = set(uniprot_wd40s[0])
    set2 = set(uniprot_wd40s[1])
    set3 = set(uniprot_wd40s[2])
    set100 = len(set1.difference(set2.union(set3)))
    set110 = len(set1.intersection(set2).difference(set3))
    set010 = len(set2.difference(set1.union(set3)))
    set101 = len(set1.intersection(set3).difference(set2))
    set111 = len(set1.intersection(set2).intersection(set3))
    set011 = len(set2.intersection(set3).difference(set1))
    set001 = len(set3.difference(set1.union(set2)))
    v = venn3(subsets={'100':1, '110':1, '010':1, '101':1, '111':1, '011':1, '001':1}, set_labels = ('UniProt_Repeat', 'UniProt_Keyword','UniProt_Family' ))
    v.get_label_by_id('100').set_text(str(set100))
    v.get_label_by_id('110').set_text(str(set110))
    v.get_label_by_id('010').set_text(str(set010))
    v.get_label_by_id('101').set_text(str(set101))
    v.get_label_by_id('111').set_text(str(set111))
    v.get_label_by_id('011').set_text(str(set011))
    v.get_label_by_id('001').set_text(str(set001))
    plt.savefig('UniProt.png',dpi=300)
    plt.close('all')

    print set3.difference(set1)

if __name__ == "__main__":
    main()

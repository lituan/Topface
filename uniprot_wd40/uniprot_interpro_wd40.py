#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
query wd40 in uniprot
"""
import sys
import os
import urllib
import urllib2
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3,venn2
from collections import OrderedDict
from multiprocessing import Pool
import cPickle as pickle


def uniprot_wd40(key='pfam',pdb=False):
    """
    use annotations from different database to query WD40 in uniprot
    return a list of uniprot accesions
    """
    keywords = ['pfam','pfam1','smart','supfam','prosite1','prosite2','prosite3','interpro_repeat','interpro_domain','interpro_family','uniprot_repeat','uniprot_keyword','uniprot_family']
    if   key == 'pfam':
        query = 'database:(type:pfam id:PF00400) or database:(type:pfam id:PF12894) or database:(type:pfam id:PF16529) or database:(type:pfam id:PF16756) or database:(type:pfam id:PF17005)'
        query = 'database:(type:pfam id:PF00400)'
    elif key == 'pfam1':
        query = 'database:(type:pfam id:PF17005)'
    elif key == 'smart':
        query = 'database:(type:smart id:SM00320)'
    elif key == 'supfam':
        query = 'database:(type:supfam id:SSF50978)'
    elif key == 'prosite1':
        query = 'database:(type:prosite id:PS00678)'
    elif key == 'prosite2':
        query = 'database:(type:prosite id:PS50082)'
    elif key == 'prosite3':
        query = 'database:(type:prosite id:PS50294)'
    elif key == 'interpro_repeat':
        query = 'database:(type:interpro id:IPR001680)'
    elif key == 'interpro_domain':
        query = 'database:(type:interpro id:IPR017986)'
    elif key == 'interpro_family':
        query = 'database:(type:interpro id:IPR031544)'
    elif key == 'uniprot_keyword':
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

    print key,'search finished',len(lines)
    return key,lines

def write_lis_lis(lis_lis,filename,cols=[]):
    """align nested list to print a table"""
    lis_lis = [lis if lis else ['    '] for lis in lis_lis]
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

def venn_standard_3(set1,set2,set3,label1,label2,label3,fname):
    sns.set_color_codes('bright')
    set1 = set(set1)
    set2 = set(set2)
    set3 = set(set3)
    set100 = len(set1.difference(set2.union(set3)))
    set110 = len(set1.intersection(set2).difference(set3))
    set010 = len(set2.difference(set1.union(set3)))
    set101 = len(set1.intersection(set3).difference(set2))
    set111 = len(set1.intersection(set2).intersection(set3))
    set011 = len(set2.intersection(set3).difference(set1))
    set001 = len(set3.difference(set1.union(set2)))
    v = venn3(subsets={'100':1, '110':1, '010':1, '101':1, '111':1, '011':1, '001':1}, set_labels = (label1,label2,label3))
    v.get_label_by_id('100').set_text(str(set100))
    v.get_label_by_id('110').set_text(str(set110))
    v.get_label_by_id('010').set_text(str(set010))
    v.get_label_by_id('101').set_text(str(set101))
    v.get_label_by_id('111').set_text(str(set111))
    v.get_label_by_id('011').set_text(str(set011))
    v.get_label_by_id('001').set_text(str(set001))
    plt.savefig(fname+'.png',dpi=300)
    plt.close('all')

def venn_standard_2(set1,set2,label1,label2,fname):
    sns.set_color_codes('bright')
    set1 = set(set1)
    set2 = set(set2)
    set10 = len(set1.difference(set2))
    set12 = len(set1.intersection(set2))
    set02 = len(set2.difference(set1))
    v = venn2(subsets={'10':1,'11':1,'01':1},set_labels=(label1,label2))
    v.get_label_by_id('10').set_text(str(set10))
    v.get_label_by_id('11').set_text(str(set12))
    v.get_label_by_id('01').set_text(str(set02))
    plt.savefig(fname+'.png',dpi=300)
    plt.close('all')


def main():
    keywords = ['pfam','pfam1','smart','supfam','prosite1','prosite2','prosite3','interpro_repeat','interpro_domain','interpro_family','uniprot_repeat','uniprot_keyword','uniprot_family']
    # p = Pool(13)
    # result = p.map(uniprot_wd40,keywords)
    # p.close()
    # wd40s = []
    # for k in keywords:
        # for r,v in result:
            # if r == k:
                # wd40s.append(v)

    # pickle.dump(wd40s,open('uniprot_wd40s.pickle','w'))
    wd40s = pickle.load(open('uniprot_wd40s.pickle'))

    # interpro_domain and supfam,prosite3
    venn_standard_3(wd40s[8],wd40s[3],wd40s[6],'IPR017986','SSF50978','PS50294','interpro_domain_supfam_prosite3')
    venn_standard_2(wd40s[8],set.union(*[wd40s[3],wd40s[6]]),'IPR017986','SSF50978_PS50294','interpro_domain-supfam_prosite3')

    # interpro_repat and smart,pfam,prosite2
    venn_standard_2(wd40s[7],set.union(*[wd40s[2],wd40s[0],wd40s[5]]),'IPR001680','PS50082_SM00320_PF00400','interpro_repeat_smart_pfam_prosite2')

    # interpro family and pfam1
    venn_standard_2(wd40s[9],wd40s[1],'IPR031544','PF17005','interpro_family_pfam1')

    # interpro
    venn_standard_3(wd40s[7],wd40s[8],wd40s[9],'IPR001680','IPR017986','IPR031544','interpro')

    # interpro uniprot
    set1 = set.union(*[wd40s[7],wd40s[8],wd40s[9]])
    set2 = set.union(*[wd40s[10],wd40s[11],wd40s[12]])
    venn_standard_2(set1,set2,'InterPro','UniProt','interpro_uniprot')

    # interpro_repeat uniprot_repeat
    venn_standard_2(wd40s[7],wd40s[10],'IPR001680','UniProt_Repeat','interpro_repeat_uniprot_repeat')

    # interpro_domain uniprot_keyword
    venn_standard_2(wd40s[8],wd40s[11],'IPR017986','UniProt_Keyword','interpro_domain_uniprot_keyword')

    # uniprot
    venn_standard_3(wd40s[10],wd40s[11],wd40s[12],'UniProt_Repeat','UniProt_Keyword','UniProt_family','uniprot')


if __name__ == "__main__":
    main()


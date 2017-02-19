#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
use UniProt REST service to search sequences annotated as 'WD40'
"""

import sys
import os
import urllib
import urllib2
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3,venn2
from multiprocessing import Pool

import lt
@lt.run_time
def uniprot_wd40(key='pfam',pdb=False):
    """
    use annotations from different database to query WD40 in uniprot
    return a list of uniprot accesions
    """
    if   key == 'pfam':
        query = 'database:(type:pfam id:PF00400) or database:(type:pfam id:PF12894) or database:(type:pfam id:PF16529) or database:(type:pfam id:PF16756)'
    elif key == 'smart':
        query = 'database:(type:smart id:SM00320)'
    elif key == 'supfam':
        query = 'database:(type:supfam id:SSF50978)'
    elif key == 'interpro_repeat':
        query = 'database:(type:interpro id:IPR001680)'
    elif key == 'interpro_domain':
        query = 'database:(type:interpro id:IPR017986)'
    elif key == 'uniprot_keyword':
        query = 'keyword:"WD repeat"'
    elif key == 'uniprot_repeat':
        query = 'annotation:(type:repeat wd)'
    elif key == 'prosite1':
        query = 'database:(type:prosite id:PS00678)'
    elif key == 'prosite2':
        query = 'database:(type:prosite id:PS50082)'
    elif key == 'prosite3':
        query = 'database:(type:prosite id:PS50294)'
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

    return key,lines


@lt.run_time

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

@lt.run_time
def main():
    keywords = ['interpro_repeat','interpro_domain','pfam','smart','supfam','uniprot_repeat','uniprot_keyword','prosite1','prosite2','prosite3']
    p = Pool(10)
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


    # barplot of wd40s annotated by different database or method
    lt.pickle_dump(wd40s,'wd40s')
    f,ax = plt.subplots()
    keys = ['Interpro_repeat','Interpro_domain','Pfam','SMART','Superfamily','UniProt_repeat','UniProt_keyword','Prosite1','Prosite2','Prosite3','WDSP']
    wd = pd.DataFrame({'Database':keys,'Num':map(len,wd40s)})
    wd = wd.sort_values('Num',ascending=False)
    sns.set_color_codes('pastel')
    sns.barplot(y='Database',x='Num',data=wd,color='b')
    ax.set(xlabel='Database',ylabel='Num',title='WD40 Annotated by Different Database')
    # plt.xticks(roation=90)
    plt.savefig('wd40_annotated_by_different_databases_accs',dpi=300)
    plt.close('all')
    write_lis_lis(wd40s,'wd40_annotated_by_different_databases_accs',keys)

    # percent table (heatmap) of wd40s shared by two different database or method
    table = []
    keys = ['Interpro_repeat','Interpro_domain','Pfam','SMART','Superfamily','UniProt_repeat','UniProt_keyword','Prosite1','Prosite2','Prosite3','WDSP']
    for w in wd40s:
        row = [len(set(w).intersection(set(wr)))*1.0/len(w) for wr in wd40s]
        table.append(row)
    data = pd.DataFrame(table,columns=keys,index=keys)
    fig = plt.figure(figsize=(9,10))
    ax = fig.add_subplot(111)
    h = sns.heatmap(data,annot=True,fmt='.2f',cmap='Blues')
    # h.figure.tight_layout()
    h.figure.subplots_adjust(top=0.9,bottom=0.13,left=0.13,right=0.9)
    ax.set_xticklabels(keys,rotation=90)
    ax.set_yticklabels(keys[::-1],rotation=0)
    ax.set_title('Comaration of Different Annotation Methods')
    plt.savefig('Comaration_of Different_Annotation_Methods.png',dpi=300)
    plt.close('all')

    # interpro_repat and smart,pfam,prosite2
    sns.set_color_codes('bright')
    venn2([wd40s[0],wd40s[3]],['Interpro_repeat','SMART'],set_colors=['r','g'])
    plt.savefig('interpro_repat_smart.png',dpi=300)
    plt.close('all')
    venn2([wd40s[0],wd40s[2]],['Interpro_repeat','Pfam'],set_colors=['r','g'])
    plt.savefig('interpro_repat_pfam.png',dpi=300)
    plt.close('all')
    venn2([wd40s[0],wd40s[8]],['Interpro_repeat','Prosite2'],set_colors=['r','g'])
    plt.savefig('interpro_repat_prosite2.png',dpi=300)
    plt.close('all')
    pfam_smart_prosite2 = set.union(*[wd40s[2],wd40s[3],wd40s[8]])
    venn2([wd40s[0],pfam_smart_prosite2],['Interpro_repeat','Pfam_SMART_Prosite2'],set_colors=['r','g'])
    plt.savefig('interpro_repat_pfam_smart_prosite2.png',dpi=300)
    plt.close('all')

    # interpro_domain and supfam,prosite2
    venn2([wd40s[1],wd40s[9]],['Interpro_domain','Prosite3'],set_colors=['r','g'])
    plt.savefig('interpro_domain_prosite3.png',dpi=300)
    plt.close('all')
    venn2([wd40s[1],wd40s[4]],['Interpro_domain','Supfamily'],set_colors=['r','g'])
    plt.savefig('interpro_domain_supfamily.png',dpi=300)
    plt.close('all')
    supfamily_prosite3 = set.union(*[wd40s[4],wd40s[9]])
    venn2([wd40s[1],supfamily_prosite3],['Interpro_domain','Supfamily_Prosite3'],set_colors=['r','g'])
    plt.savefig('interpro_domain_supfamily_prosite3.png',dpi=300)
    plt.close('all')

    # uniprot_repat and uniprot_keyword
    venn2([wd40s[6],wd40s[5]],['Uniprot_keyword','Uniprot_repeat'],set_colors=['r','g'])
    plt.savefig('uniprot_keyword_repeat.png',dpi=300)
    plt.close('all')

    # uniprot and interpro
    interpro = set.union(*[wd40s[0],wd40s[1]])
    venn2([interpro,wd40s[5]],['Interpro','Uniprot_repeat'],set_colors=['r','g'])
    plt.savefig('interpro_uniprot_repeat.png',dpi=300)
    plt.close('all')
    venn2([interpro,wd40s[6]],['Interpro','Uniprot_keyword'],set_colors=['r','g'])
    plt.savefig('interpro_uniprot_keyword.png',dpi=300)
    plt.close('all')

    # prosite
    venn3([wd40s[7],wd40s[8],wd40s[9]],['Prosite1','Prosite2','Prosite3'],set_colors=['r','g','b'])
    plt.savefig('prosite.png',dpi=300)
    plt.close('all')

    keys = keys[2:]
    wd40s = wd40s[2:]

    total = set.union(*map(set,wd40s))
    wd40s_score = [[] for i in range(9)]
    def acc_score(acc):
        i = 0
        for w in wd40s:
            if acc in w:
                i += 1
        return i
    for acc in total:
        num = acc_score(acc)
        wd40s_score[num-1].append(acc)

    lt.pickle_dump(wd40s_score,'wd40s_score')
    f,ax = plt.subplots()
    wd = pd.DataFrame({'Database Score':range(1,10),'Num':map(len,wd40s_score)})
    sns.set_color_codes('pastel')
    sns.barplot(x='Database Score',y='Num',data=wd,color='b')
    ax.set(xlabel='Database Score',ylabel='Num',title='Annotation Score of WD40 in UniProt')
    plt.savefig('wd40_annotation_score_in_uniprot_accs',dpi=300)
    plt.close('all')
    write_lis_lis(wd40s_score,'wd40_annotation_score_in_uniprot_accs',[str(i) for i in range(1,10)])

    # write_lis_lis(wd40s,'uniprot_wd40',keywords)


if __name__ == "__main__":
    main()

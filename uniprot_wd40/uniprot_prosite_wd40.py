#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
query wd40 using prosite annotation
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
import cPickle as pickle

wd40_prosites = {
'PS50082' : 'PS50082',
'PS50294' : 'PS50292',
'PS00678' : 'PS00678'
}

def uniprot_query(key,prosite_id,pdb=False):
    """
    use annotations from different database to query WD40 in uniprot
    return a list of uniprot accesions
    """
    query = 'database:(type:prosite id:'+prosite_id+')'
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
    print 'uniprot search ',key,' is finished'
    return lines

def query(prosite):
    key,prosite_id = prosite[0],prosite[1]
    for i in range(10):
        try:
            ids = uniprot_query(key,prosite_id)
            results = [key,ids]
            print 'query ',key,'is successful'
            print len(ids),' entries got'
            break
        except:
            continue
    return results

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

def main():
    # prosites = [[key,prosite_id] for key,prosite_id in wd40_prosites.iteritems() ]
    # from multiprocessing import Pool
    # p = Pool(3)
    # results = p.map(query,prosites)
    # p.close()
    # keys = [k for k,v in results]
    # values = [v for k,v in results]
    # write_lis_lis(values,'pfam_propellers',cols=keys)

    # pickle.dump(results,open('wd4o_prosites.pickle','w'))
    results = pickle.load(open('wd4o_prosites.pickle'))

    # pfams = pickle.load(open('pfam_propellers.pickle','r'))
    results = OrderedDict([[k,v] for k,v in results if len(v) > 0]) # filter empty entries
    print 'types of propellers: ',len(results)


    #plot all propellers
    f,ax = plt.subplots(figsize=(6,6))
    sns.set_context(rc={'patch.linewidth':0.0})
    sns.set_color_codes('pastel')
    wd = pd.DataFrame({'Prosite Entries':results.keys(),'UniProt Sequence Num':map(len,results.values())})
    wd = wd.sort_values('UniProt Sequence Num',ascending=False)
    h = sns.barplot(y='UniProt Sequence Num',x='Prosite Entries',data=wd,color='b')
    h.figure.subplots_adjust(top=0.88,bottom=0.10,left=0.15,right=0.95)
    ax.set(ylabel='UniProt Sequence Num',xlabel='Prosite Entries',title='Num of WD40s in Uniprot according to Different Prosite Annotation')
    # plt.xticks(rotation='vertical')
    plt.savefig('num_of_different_wd40s_by_prosites',dpi=300)
    plt.close('all')

    #plot wd40s
    f,ax = plt.subplots(figsize=(7,6))
    sns.set_context(rc={'patch.linewidth':0.0})
    wd40_names = ['PS50082','PS50294','PS00678']
    wd40_names = [k for k in wd40_names if k in results.keys()]
    wd40s_large = [results[k] for k in wd40_names]
    wd40s_total = set.union(*map(set,wd40s_large))
    wd = pd.DataFrame({'Prosite entries':wd40_names,'UniProt Sequence Num':map(len,[wd40s_total for i in range(len(wd40_names))])})
    sns.barplot(x='Prosite entries',y='UniProt Sequence Num',data=wd,color='lightblue')
    sns.set_color_codes('muted')
    wd = pd.DataFrame({'Prosite entries':wd40_names,'UniProt Sequence Num':map(len,wd40s_large)})
    wd = wd.sort_values('UniProt Sequence Num',ascending=False)
    sns.barplot(x='Prosite entries',y='UniProt Sequence Num',data=wd,color='b')
    ax.set(xlabel='Prosite entries',ylabel='UniProt Sequence Num',title='Num of WD40s in Uniprot according to Different Prosite Annotation')
    plt.savefig('wd40_in_uniprot_by_prosite',dpi=300)
    plt.close('all')

    with open('wd40s_table.txt','w') as w_f:
        for name,num in zip(wd40_names,map(len,wd40s_large)):
            print >> w_f, '{0:<30}{1:<}'.format(name,num)

    # plot wd40s venns
    sns.set_color_codes('bright')
    set1 = set(wd40s_large[0])
    set2 = set(wd40s_large[1])
    set10 = len(set1.difference(set2))
    set12 = len(set1.intersection(set2))
    set02 = len(set2.difference(set1))
    v = venn2(subsets={'10':4,'11':1,'01':4},set_labels=(wd40_names[0],wd40_names[1]))
    v.get_label_by_id('10').set_text(str(set10))
    v.get_label_by_id('11').set_text(str(set12))
    v.get_label_by_id('01').set_text(str(set02))
    plt.savefig(wd40_names[0]+'_'+wd40_names[1]+'.png',dpi=300)
    plt.close('all')

    set1 = set(wd40s_large[0])
    set2 = set(wd40s_large[2])
    set10 = len(set1.difference(set2))
    set12 = len(set1.intersection(set2))
    set02 = len(set2.difference(set1))
    v = venn2(subsets={'10':1,'11':4,'01':1},set_labels=(wd40_names[0],wd40_names[2]))
    v.get_label_by_id('10').set_text(str(set10))
    v.get_label_by_id('11').set_text(str(set12))
    v.get_label_by_id('01').set_text(str(set02))
    plt.savefig(wd40_names[0]+'_'+wd40_names[2]+'.png',dpi=300)
    plt.close('all')

    set1 = set(wd40s_large[1])
    set2 = set(wd40s_large[2])
    set10 = len(set1.difference(set2))
    set12 = len(set1.intersection(set2))
    set02 = len(set2.difference(set1))
    v = venn2(subsets={'10':4,'11':1,'01':4},set_labels=(wd40_names[1],wd40_names[2]))
    v.get_label_by_id('10').set_text(str(set10))
    v.get_label_by_id('11').set_text(str(set12))
    v.get_label_by_id('01').set_text(str(set02))
    plt.savefig(wd40_names[1]+'_'+wd40_names[2]+'.png',dpi=300)
    plt.close('all')

    set1 = set(wd40s_large[0])
    set2 = set(wd40s_large[1])
    set3 = set(wd40s_large[2])
    set100 = len(set1.difference(set2.union(set3)))
    set110 = len(set1.intersection(set2).difference(set3))
    set010 = len(set2.difference(set1.union(set3)))
    set101 = len(set1.intersection(set3).difference(set2))
    set111 = len(set1.intersection(set2).intersection(set3))
    set011 = len(set2.intersection(set3).difference(set1))
    set001 = len(set3.difference(set1.union(set2)))
    v = venn3(subsets={'100':1, '110':1, '010':1, '101':1, '111':1, '011':1, '001':1}, set_labels = (wd40_names[0], wd40_names[1], wd40_names[2] ))
    v.get_label_by_id('100').set_text(str(set100))
    v.get_label_by_id('110').set_text(str(set110))
    v.get_label_by_id('010').set_text(str(set010))
    v.get_label_by_id('101').set_text(str(set101))
    v.get_label_by_id('111').set_text(str(set111))
    v.get_label_by_id('011').set_text(str(set011))
    v.get_label_by_id('001').set_text(str(set001))
    plt.savefig(wd40_names[0]+'_'+wd40_names[1]+'_'+wd40_names[2]+'.png',dpi=300)


    # plot heatmap of wd40s shared by different annotation
    sns.set_color_codes('pastel')
    table = []
    keys = wd40_names
    for w in wd40s_large:
        row = [len(set(w).intersection(set(wr)))*1.0/len(w) for wr in wd40s_large]
        table.append(row)
    data = pd.DataFrame(table,columns=keys,index=keys)
    fig = plt.figure(figsize=(7,8))
    ax = fig.add_subplot(111)
    h = sns.heatmap(data,annot=True,fmt='.2f',cmap='Blues')
    # h.figure.tight_layout()
    h.figure.subplots_adjust(top=0.9,bottom=0.05,left=0.18,right=0.98)
    ax.set_xticklabels(keys,rotation=0)
    ax.set_yticklabels(keys[::-1],rotation=0)
    ax.set_title('Comaration of Different Annotation Methods')
    plt.savefig('Comaration_of Different_Annotation_Methods.png',dpi=300)
    plt.close('all')


if __name__ == "__main__":
    main()


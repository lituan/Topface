#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""
use different strategies to search WD40 structures in rcsb

"""

import sys
import os
import urllib
import urllib2
import httplib
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import cPickle as pickle
from matplotlib_venn import venn3,venn2
from multiprocessing import Pool


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
    for i in range(20):
        try:
            print 'try time',i+1
            req = urllib2.Request(url,data)
            response = urllib2.urlopen(req)
            r = response.readlines()
            accs = set([line.rstrip('\r\n') for line in r])
            print 'uniprot search ',key,' is finished'
            break
        except httplib.IncompleteRead,e:
            print e
            continue
    return key,accs

def rcsb_acc(accs):
    """
    use beta_sheet, chain length to make sure we get WD40s other than others
    accs be a list, ['Q969H0,P07834']
    return string format: '1A0R:1,1B9X:1,1B9Y:1,1C15:1'
    """
    url = 'http://www.rcsb.org/pdb/rest/search'
    accs = ','.join(accs)

    query ="""
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>
    <description>Simple query for a list of Uniprot Accession IDs: """+accs+"""</description>
    <accessionIdList>"""+accs+"""</accessionIdList>
    </orgPdbQuery>
    """
    req = urllib2.Request(url,data=query)
    response = urllib2.urlopen(req)
    pdbids = response.read()
    pdbids = pdbids.replace('\n',',')
    return pdbids


def get_wdsp_acc():
    with open('wdsp_uniprot_id.txt') as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        lines = [line.split()[0] for line in lines]
        return lines

def rcsb_uniprot():
    keywords = ['pfam','smart','supfam','prosite','uniprot']
    p = Pool(5)
    result = p.map(uniprot_wd40,keywords)
    p.close()
    wd40s = []
    for k in keywords:
        for r,v in result:
            if r == k:
                wd40s.append(v)
    wdsp = get_wdsp_acc()
    wd40s.append(wdsp)
    keywords.append('wdsp')

    total = set.union(*map(set,wd40s))
    uniprot_pdbids = rcsb_acc(total)

    with open('uniprot_wd40.txt','w') as w_f:
        for a in total:
            print >> w_f, a

    with open('uniprot_wd40_pdb.txt','w') as w_f:
        print >> w_f,uniprot_pdbids

    return uniprot_pdbids


def rcsb_pfam():
    url = 'http://www.rcsb.org/pdb/rest/search'

    begin_query="""
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.PfamIdQuery</queryType>"""

    pfam1 = """
    <description>Pfam Accession Number PF00400 WD domain, G-beta repeat</description>
    <pfamID>PF00400</pfamID>
    """
    pfam2 = """
    <description>Pfam Accession Number PF12894 Anaphase-promoting complex subunit 4 WD40 domain</description>
    <pfamID>PF12894</pfamID>
    """
    pfam3 = """
    <description>Pfam Accession Number PF16756 </description>
    <pfamID>PF16756</pfamID>
    """
    pfam4 = """
    <description>Pfam Accession Number PF16529 </description>
    <pfamID>PF12894</pfamID>
    """
    pfam5 = """
    <description>Pfam Accession Number PF17005 </description>
    <pfamID>PF17005</pfamID>
    """

    end_query = """
    </orgPdbQuery>
    """
    pfams = [pfam1,pfam2,pfam3,pfam4,pfam5]
    pdbids = ''
    for pfam in pfams:
        req = urllib2.Request(url,data=begin_query+pfam+end_query)
        response = urllib2.urlopen(req)
        result_pdb = response.read()
        pdbids += result_pdb.replace('\n',',')
    # get customed report
    return pdbids

def rcsb_scop():
    url = 'http://www.rcsb.org/pdb/rest/search'
    scop_query = """
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.TreeQuery</queryType>
    <description>ScopTree Search for WD40 repeat-like ( also contains 8-bladed propellers )</description>
    <t>11</t>
    <n>50978</n>
    <nodeDesc>WD40 repeat-like ( also contains 8-bladed propellers )</nodeDesc>
    </orgPdbQuery>
    """
    req = urllib2.Request(url,data=scop_query)
    response = urllib2.urlopen(req)
    result_pdb = response.read()
    pdbids = result_pdb.replace('\n',',')
    return pdbids
    # get customed report

def rcsb_txt():

    url = 'http://www.rcsb.org/pdb/rest/search'
    txt_query = """
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>
    <description><![CDATA[Text Search for: wd40 or "wd repeat" or "wd 40"]]></description>
    <keywords><![CDATA[wd40 OR "wd repeat" OR "wd 40"]]></keywords>
    </orgPdbQuery>
    """
    req = urllib2.Request(url,data=txt_query)
    response = urllib2.urlopen(req)
    result_pdb = response.read()
    pdbids = result_pdb.replace('\n',',')
    # get customed report
    return pdbids


def main():
    uniprot_pdbids = rcsb_uniprot()
    scop_pdbids = rcsb_scop()
    pfam_pdbids = rcsb_pfam()
    txt_pdbids = rcsb_txt()

    uniprot = set([u.split(':')[0] for u in uniprot_pdbids.split(',') if u])
    scop = set([u.split(':')[0] for u in scop_pdbids.split(',') if u])
    pfam = set([u.split(':')[0] for u in pfam_pdbids.split(',') if u])
    txt = set([u.split(':')[0] for u in txt_pdbids.split(',') if u])

    pickle.dump([uniprot,scop,pfam,txt],open('rcsb_wd40.pickle','w'))
    uniprot,scop,pfam,txt = pickle.load(open('rcsb_wd40.pickle'))
    write_lis_lis([uniprot,scop,pfam,txt],'rcsb_wd40',['uniprot','scop','pfam','txt'])

    # plot heatmap of WD40s shared by different methods
    sns.set_color_codes('pastel')
    table = []
    keys = ['UniProt','SCOP','Pfam','Text']
    total = [uniprot,scop,pfam,txt]
    for w in total:
        row = [len(w.intersection(wr))*1.0/len(w) for wr in total]
        table.append(row)
    data = pd.DataFrame(table,columns=keys,index=keys)
    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111)
    h = sns.heatmap(data,annot=True,fmt='.2f',cmap='Blues')
    h.figure.subplots_adjust(top=0.9,bottom=0.13,left=0.10,right=0.96)
    ax.set_xticklabels(keys,rotation=0)
    ax.set_yticklabels(keys[::-1],rotation=0)
    ax.set_title('Comaration of Different Annotation Methods')
    plt.savefig('Comaration_of Different_Annotation_Methods.png',dpi=300)
    plt.close('all')


    def venn_2(set1,set2,label1,label2):
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
        plt.savefig(label1+'_'+label2+'.png',dpi=300)
        plt.close('all')

    venn_2(uniprot,scop,'UniProt','SCOP')
    venn_2(uniprot,pfam,'UniProt','Pfam')
    venn_2(uniprot,txt,'UniProt','Text')
    spt = set.union(*[scop,pfam,txt])
    venn_2(uniprot,spt,'UniProt','Pfam_SCOP_Text')
    venn_2(pfam,txt,'Pfam','Text')
    venn_2(pfam,scop,'Pfam','SCOP')

    set1 = set(uniprot)
    set2 = set(pfam)
    set3 = set(scop)
    set100 = len(set1.difference(set2.union(set3)))
    set110 = len(set1.intersection(set2).difference(set3))
    set010 = len(set2.difference(set1.union(set3)))
    set101 = len(set1.intersection(set3).difference(set2))
    set111 = len(set1.intersection(set2).intersection(set3))
    set011 = len(set2.intersection(set3).difference(set1))
    set001 = len(set3.difference(set1.union(set2)))
    v = venn3(subsets={'100':1, '110':1, '010':1, '101':1, '111':1, '011':1, '001':1}, set_labels = ('UniProt', 'Pfam', 'SCOP' ))
    v.get_label_by_id('100').set_text(str(set100))
    v.get_label_by_id('110').set_text(str(set110))
    v.get_label_by_id('010').set_text(str(set010))
    v.get_label_by_id('101').set_text(str(set101))
    v.get_label_by_id('111').set_text(str(set111))
    v.get_label_by_id('011').set_text(str(set011))
    v.get_label_by_id('001').set_text(str(set001))
    plt.savefig('UniProt_Pfam_SCOP.png',dpi=300)
    plt.close('all')

    with open('not_uniprot.txt','w') as w_f:
        print >> w_f, 'in txt out of uniprot'
        for p in set(txt).difference(set(uniprot)):
            print >> w_f,p
        print >> w_f, 'in pfam out of uniprot'
        for p in set(pfam).difference(set(uniprot)):
            print >> w_f,p
        print >> w_f, 'in scop out of uniprot'
        for p in set(scop).difference(set(uniprot)):
            print >> w_f,p

    # plot barplot
    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111)
    total = set.union(*map(set,[uniprot,pfam,txt,scop]))
    sns.set_color_codes('pastel')
    methods = ['UniProt','Pfam','Text','SCOP']
    wd = pd.DataFrame({'Methods':methods,'PDB Num':map(len,[total,total,total,total])})
    sns.barplot(x='Methods',y='PDB Num',data=wd,color='b')
    sns.set_color_codes('muted')
    wd = pd.DataFrame({'Methods':methods,'PDB Num':map(len,[uniprot,pfam,txt,scop])})
    sns.barplot(x='Methods',y='PDB Num',data=wd,color='b')
    ax.set(xlabel='Methods',ylabel='PDB Num',title='WD40 Structures in RCSB')
    plt.savefig('wd40_in_RCSB_pdbs',dpi=300)
    plt.close('all')

if __name__ == "__main__":
    main()







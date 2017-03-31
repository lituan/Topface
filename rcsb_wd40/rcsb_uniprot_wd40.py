#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""
use different strategies to search WD40 structures in rcsb
all search has additional restrictions: resolution below 3.5, chain length longer than 150, beta strands more than 14

caution: use a low resolution as restriction will eliminate EM structures and NMR structures
"""

import os
import sys
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


def uniprot_wd40(key,pdb=False):
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
    for i in range(10):
        try:
            print 'try time',i+1
            req = urllib2.Request(url,data)
            response = urllib2.urlopen(req)
            r = response.readlines()
            accs = set([line.rstrip('\r\n') for line in r])
            print 'uniprot search ',key,' is finished'
            return key,list(accs)
        except e:
            print e
            continue


def rcsb_acc(accs,beta,chain_len,resolution):
    """
    use beta_sheet, chain length to make sure we get WD40s other than others
    accs be a list, ['Q969H0,P07834']
    return string format: '1A0R:1,1B9X:1,1B9Y:1,1C15:1'
    """
    url = 'http://www.rcsb.org/pdb/rest/search'
    accs = ','.join(accs)

    query ="""
    <orgPdbCompositeQuery version="1.0">
    <queryRefinement>
    <queryRefinementLevel>0</queryRefinementLevel>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>
    <description>Simple query for a list of Uniprot Accession IDs: """+accs+"""</description>
    <accessionIdList>"""+accs+"""</accessionIdList>
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>1</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.SecondaryStructureQuery</queryType>
    <description>Secondary structure has:  """+ str(beta) +""" or more Beta Sheets</description>
    <polyStats.helixPercent.comparator>between</polyStats.helixPercent.comparator>
    <polyStats.helixCount.comparator>between</polyStats.helixCount.comparator>
    <polyStats.sheetPercent.comparator>between</polyStats.sheetPercent.comparator>
    <polyStats.sheetCount.comparator>between</polyStats.sheetCount.comparator>
    <polyStats.sheetCount.min>""" + str(beta) + """</polyStats.sheetCount.min>
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>2</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.SequenceLengthQuery</queryType>
    <description>Sequence Length is between """ + str(chain_len) + """ and more </description>
    <v_sequence.chainLength.min>""" + str(chain_len) + """</v_sequence.chainLength.min>
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>3</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ResolutionQuery</queryType>
    <description>Resolution is """ + str(resolution) +""" or less </description>
    <refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>
    <refine.ls_d_res_high.max>""" + str(resolution) + """</refine.ls_d_res_high.max>
    </orgPdbQuery>
    </queryRefinement>
    </orgPdbCompositeQuery>
    """
    req = urllib2.Request(url,data=query)
    response = urllib2.urlopen(req)
    pdbids = response.read()
    pdbids = pdbids.replace('\n',',')
    return pdbids

def rcsb_customreport(pdbids):
    """
    pdbids be a string '1A0R:1,1B9X:1,1B9Y:1,1C15:1,1CWW:1,1CY5:1'
    return format:
    [['structureId','chainId','uniprotAcc','resolution','chainLength','releaseDate'],
     ['"1A0R"', '"B"', '"P62871"', '"2.8"', '"340"', '"1998-12-30"'],
     ['"1B9X"', '"A"', '"P62871"', '"3.0"', '"340"', '"1999-02-23"']]
    """

    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbids,
        'customReportColumns':'structureId,uniprotAcc,entityId,resolution,chainLength,releaseDate',
        'service':'wsfile',
        'format':'csv',
    }
    data = urllib.urlencode(data)
    req = urllib2.Request(url,data)
    response = urllib2.urlopen(req)
    lines = response.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line for line in lines if line]
    lines = [line.split(',') for line in lines]
    lines = [[w.strip('"') for w in line] for line in lines]
    lines = [[line[0],line[1],line[2],line[3],float(line[4]),int(line[5]),line[6]] for line in lines[1:]]
    # line format: '3V7D', 'D', 'P07834', '2', 2.3', 464, '2012-05-02'
    return lines


def rcsb_acc_customreport(accs,beta,chain_len,resolution):
    """
    when use uniprot accessions to query pdbs and use returned pdbs to get custom report, some chimera pdb will return unrelated uniprot accessions
    """
    pdbids = rcsb_acc(accs,beta,chain_len,resolution)
    lines = rcsb_customreport(pdbids)
    # remove accs associated with chimera pdbs
    new_lines = []
    for line in lines:
        if '#' in line[2]:
            line2s = line[2].split('#')
            for l in line2s:
                if l in accs:
                    new_lines.append(line[0:2]+[l]+line[3:])
        else:
            new_lines.append(line)
    # line format: '3V7D', 'D', 'P07834', '2', 2.3', 464, '2012-05-02'
    return new_lines

def get_wdsp_acc():
    with open('wdsp_uniprot_id.txt') as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        lines = [line.split()[0] for line in lines]
        return lines

def rcsb_uniprot(beta,chain_len,resolution):
    # keywords = ['pfam','smart','supfam','prosite','uniprot']
    # p = Pool(5)
    # result = p.map(uniprot_wd40,keywords)
    # p.close()
    # wd40s = []
    # for k in keywords:
        # for r,v in result:
            # if r == k:
                # wd40s.append(v)
    # wdsp = get_wdsp_acc()
    # wd40s.append(wdsp)
    # keywords.append('wdsp')
    # pickle.dump(wd40s,open('wd40s.pickle','w'))
    wd40s = pickle.load(open('wd40s.pickle'))

    total = set.union(*map(set,wd40s))
    def acc_score(acc):
        i = 0
        for w in wd40s:
            if acc in w:
                i += 1
        return i

    # report = rcsb_acc_customreport(total,beta,chain_len,resolution)
    # pdb_scores = []
    # for p in report:
        # pdb_scores.append(p+[acc_score(p[2])])
    # pdb_scores = sorted(pdb_scores,key=lambda x:x[-1],reverse=True)
    # pickle.dump(pdb_scores,open('pdb_scores.pickle','w'))
    pdb_scores = pickle.load(open('pdb_scores.pickle'))
    with open('uniprot_pdb_scores.txt','w') as w_f:
        print >> w_f,'{0:<15}{1:<10}{2:<8}{3:<8}{4:<15}{5:<12}{6:<18}{7:<8}'.format('acc','pdb','chain','entity','resolution','chain_len','release','score')
        for p in pdb_scores:
            print >> w_f,'{0:<15}{1:<10}{2:<8}{3:<8}{4:<15}{5:<12}{6:<18}{7:<8}'.format(p[2],p[0],p[1],p[3],p[4],p[5],p[6],p[7])


    # plot wd40 structures annotated by different database
    total_pdb_accs = set([p[2] for p in pdb_scores])
    wd40s_pdb_accs = [[a for a in total_pdb_accs if a in w] for w in wd40s]
    pickle.dump(wd40s_pdb_accs,open('wd40_pdb_acc_databases.pickle','w'))
    sns.set_color_codes('pastel')
    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111)
    keys = ['Pfam','SMART','Superfamily','Prosite','UniProt','WDSP']
    wd = pd.DataFrame({'Methods':keys,'Protein Num':map(len,[total_pdb_accs for i in range(6)])})
    sns.barplot(y='Methods',x='Protein Num',data=wd,color='b')
    sns.set_color_codes('muted')
    wd = pd.DataFrame({'Methods':keys,'Protein Num':map(len,wd40s_pdb_accs)})
    wd = wd.sort_values('Protein Num',ascending=False)
    h = sns.barplot(y='Methods',x='Protein Num',data=wd,color='b')
    h.figure.subplots_adjust(top=0.9,bottom=0.10,left=0.21,right=0.95)
    ax.set(xlabel='Protein Num',title='WD40 Proteins (with structure) Annotated by Different Methods')
    # plt.xticks(roation=90)
    plt.savefig('wd40_structures_accs_annotated_by_different_database',dpi=300)
    plt.close('all')
    write_lis_lis(wd40s_pdb_accs,'wd40_structures_accs_annotated_by_different_database',keys)

    # plot annotation score of wd40 structures
    pdb_acc_scores = [[] for i in range(6)]
    for p in pdb_scores:
        pdb_acc_scores[p[-1]-1].append(p[2])
    pdb_acc_scores = map(set,pdb_acc_scores)
    pickle.dump(pdb_acc_scores,open('pdb_acc_scores','w'))


    sns.set_color_codes('pastel')
    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111)
    wd = pd.DataFrame({'Methods':keys,'Protein Num':map(len,[total_pdb_accs for i in range(6)])})
    sns.barplot(x='Methods',y='Protein Num',data=wd,color='b')
    sns.set_color_codes('muted')
    wd = pd.DataFrame({'Annotation Score':range(1,7),'Protein Num':map(len,pdb_acc_scores)})
    sns.barplot(x='Annotation Score',y='Protein Num',data=wd,color='b')
    ax.set(xlabel='Annotation Score',ylabel='Protein Num',title='Annotation Score of WD40 Proteins (with structure)')
    plt.savefig('wd40_structures_annotation_score_accs',dpi=300)
    plt.close('all')
    write_lis_lis(pdb_acc_scores,'wd40_structures_annotation_score_accs',[str(i) for i in range(1,7)])


def main():
    beta,chain_len,resolution = 20,200,3.5
    rcsb_uniprot(beta,chain_len,resolution)


if __name__ == "__main__":
    main()







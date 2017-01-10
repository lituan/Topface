#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""
use different strategies to search WD40 structures in rcsb
1. use pfam annotaion PF00400
2. use scop annotation
3. use txt search, 'wd40' or 'wd 40' or 'wd repeat'
4. use uniprot accs, which are annotated as wd40 by different protein family database
all search has additional restrictions: resolution below 3.5, chain length longer than 150, beta strands more than 14

for uniprot accs, there are 10 different annotations, so each acc can be given a score according to the number of annotations
caution: use resolution as restriction will eliminate EM structures and NMR structures
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
    return lines

@lt.run_time
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


@lt.run_time
def uniprot_wd40(key='pfam',pdb=False):
    """
    use annotations from different database to query WD40 in uniprot
    return a list of uniprot accesions
    """
    if   key == 'pfam':
        query = 'database:(type:pfam id:PF00400) OR database:(type:pfam id:PF12894) OR database:(type:pfam id:PF16529) OR database:(type:pfam id:PF16756)'
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
    for i in range(10):
        try:
            req = urllib2.Request(url,data)
            response = urllib2.urlopen(req)
            r = response.readlines()
            lines = set([line.rstrip('\r\n') for line in r])
            print 'uniprot search ',key,' is finished'
            break
        except IncompleteRead:
            continue
    return key,lines

@lt.run_time
def rcsb_acc(accs,beta=15,chain_len=150,resolution=3.5):
    """
    accs be a list, ['Q969H0,P07834']
    return string format: '1A0R:1,1B9X:1,1B9Y:1,1C15:1'
    """
    url = 'http://www.rcsb.org/pdb/rest/search'
    accs = ','.join(accs)

    query ="""
    <orgPdbCompositeQuery version="1.0">
    <resultCount>11</resultCount>
    <queryId>B0EC4FC2</queryId>
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

def rcsb_acc_customreport(accs,beta=15,chain_len=150,resolution=3.5):
    new_lines = [['structureId','chainId','uniprotAcc','entityId','resolution','chainLength','releaseDate']]
    for acc in accs:
        pdbids = rcsb_acc([acc],beta=15,chain_len=150,resolution=3.5)
        if pdbids:
            lines = rcsb_customreport(pdbids)
            for line in lines[1:]:
                if '#' in line[2]:
                    line2s = line[2].split('#')
                    for l in line2s:
                        if l == acc:
                            new_lines.append(line[0:2]+[l]+line[3:])
                else:
                    new_lines.append(line)
    return new_lines

def rcsb_acc_customreport(accs,beta=15,chain_len=150,resolution=3.5):
    """
    when use uniprot accessions to query pdbs and use returned pdbs to get custom report, some chimera pdb will return unrelated uniprot accessions
    """
    pdbids = rcsb_acc(accs,beta=15,chain_len=150,resolution=3.5)
    lines = rcsb_customreport(pdbids)
    # remove accs associated with chimera pdbs
    new_lines = [lines[0]]
    for line in lines[1:]:
        if '#' in line[2]:
            line2s = line[2].split('#')
            for l in line2s:
                new_lines.append(line[0:2]+[l]+line[3:])
        else:
            new_lines.append(line)
    real_lines = [line for line in new_lines[1:] if line[2] in accs]
    real_lines = [line for line in real_lines if int(line[5]) >= chain_len]
    real_lines = [line for line in real_lines if float(line[4]) <= resolution]
    real_lines = [lines[0]] + real_lines
    return pdbids,real_lines

def get_wdsp_acc():
    with open('wdsp_uniprot_id.txt') as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        lines = [line.split()[0] for line in lines]
        return lines

@lt.run_time
def rcsb_uniprot(beta=15,chain_len=150,resolution=3.5):
    keywords = ['pfam','smart','supfam','uniprot_repeat','uniprot_keyword','prosite1','prosite2','prosite3']
    p = Pool(8)
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

    lt.pickle_dump(wd40s,'wd40s')
    f,ax = plt.subplots()
    keys = ['Pfam','SMART','Superfamily','UniProt1','UniProt2','Prosite1','Prosite2','Prosite3','WDSP']
    wd = pd.DataFrame({'Database':keys,'Num':map(len,wd40s)})
    wd = wd.sort_values('Num',ascending=True)
    sns.set_color_codes('pastel')
    sns.barplot(x='Database',y='Num',data=wd,color='b')
    ax.set(xlabel='Database',ylabel='Num',title='WD40 Annotated by Different Database')
    # plt.xticks(roation=90)
    plt.savefig('wd40_annotated_by_different_databases_accs',dpi=300)
    plt.close('all')
    write_lis_lis(wd40s,'wd40_annotated_by_different_databases_accs',keys)

    total = set.union(*map(set,wd40s))
    # if an entry apears in n different querys, its score is n
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

    # uniprot_pdbids = rcsb_acc(total,beta,chain_len,resolution)
    # report = rcsb_customreport(uniprot_pdbids)
    uniprot_pdbids,report = rcsb_acc_customreport(total,beta,chain_len,resolution)
    pdb_scores = []
    for p in report[1:]:
        pdb_scores.append(p+[acc_score(p[2])])
    pdb_scores = sorted(pdb_scores,key=lambda x:x[-1],reverse=True)
    with open('uniprot_pdb_scores.txt','w') as w_f:
        print >> w_f,'{0:<15}{1:<10}{2:<8}{3:<8}{4:<15}{5:<8}{6:<18}{7:<8}'.format('acc','pdb','chain','entity','resolution','chain_len','release','score')
        for p in pdb_scores:
            print >> w_f,'{0:<15}{1:<10}{2:<8}{3:<8}{4:<15}{5:<8}{6:<18}{7:<8}'.format(p[2],p[0],p[1],p[3],p[4],p[5],p[6],p[7])


    # plot wd40 structures annotated by different database
    total_pdb = set([p[2] for p in pdb_scores])
    wd40s_pdb = [[a for a in w if a in total_pdb] for w in wd40s]
    lt.pickle_dump(wd40s_pdb,'pdb_acc_databases')
    f,ax = plt.subplots()
    keys = ['Pfam','SMART','Superfamily','UniProt1','UniProt2','Prosite1','Prosite2','Prosite3','WDSP']
    wd = pd.DataFrame({'Database':keys,'Num':map(len,wd40s_pdb)})
    wd = wd.sort_values('Num',ascending=True)
    sns.set_color_codes('pastel')
    sns.barplot(x='Database',y='Num',data=wd,color='b')
    ax.set(xlabel='Database',ylabel='Num',title='WD40 Structures Annotated by Different Database')
    # plt.xticks(roation=90)
    plt.savefig('wd40_structures_annotated_by_different_database',dpi=300)
    plt.close('all')
    write_lis_lis(wd40s_pdb,'wd40_structures_annotated_by_different_database',keys)

    # plot annotation score of wd40 structures
    pdb_acc_scores = [[] for i in range(9)]
    for p in pdb_scores:
        pdb_acc_scores[p[-1]-1].append(p[2])
    pdb_acc_scores = map(set,pdb_acc_scores)

    lt.pickle_dump(pdb_acc_scores,'pdb_acc_scores')
    f,ax = plt.subplots()
    wd = pd.DataFrame({'Database Score':range(1,10),'Num':map(len,pdb_acc_scores)})
    sns.set_color_codes('pastel')
    sns.barplot(x='Database Score',y='Num',data=wd,color='b')
    ax.set(xlabel='Database Score',ylabel='Num',title='Annotation Score of WD40 Structures')
    plt.savefig('wd40_structures_annotation_score_accs',dpi=300)
    plt.close('all')
    write_lis_lis(pdb_acc_scores,'wd40_structures_annotation_score_accs',[str(i) for i in range(1,10)])

    print 'uniprot search is finished'
    return uniprot_pdbids,pdb_scores


@lt.run_time
def rcsb_pfam(beta=15,chain_len=150,resolution=3.5):
    url = 'http://www.rcsb.org/pdb/rest/search'
    pfam_query="""<orgPdbCompositeQuery version="1.0">
    <queryRefinement>
    <queryRefinementLevel>0</queryRefinementLevel>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.PfamIdQuery</queryType>
    <description>Pfam Accession Number PF00400 WD domain, G-beta repeat</description>
    <pfamID>PF00400</pfamID>
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
    req = urllib2.Request(url,data=pfam_query)
    response = urllib2.urlopen(req)
    result_pdb = response.read()
    pdbids = result_pdb.replace('\n',',')
    # get customed report
    return pdbids

@lt.run_time
def rcsb_scop(beta=15,chain_len=150,resolution=3.5):
    url = 'http://www.rcsb.org/pdb/rest/search'
    scop_query = """
    <orgPdbCompositeQuery version="1.0">
    <queryRefinement>
    <queryRefinementLevel>0</queryRefinementLevel>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.TreeQuery</queryType>
    <description>ScopTree Search for WD40 repeat-like ( also contains 8-bladed propellers )</description>
    <nodeDesc>WD40 repeat-like ( also contains 8-bladed propellers )</nodeDesc>
    <t>11</t>
    <n>50978</n>
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
    req = urllib2.Request(url,data=scop_query)
    response = urllib2.urlopen(req)
    result_pdb = response.read()
    pdbids = result_pdb.replace('\n',',')
    return pdbids
    # get customed report

@lt.run_time
def rcsb_txt(beta=15,chain_len=150,resolution=3.5):

    url = 'http://www.rcsb.org/pdb/rest/search'
    txt_query = """<orgPdbCompositeQuery version="1.0">
    <queryRefinement>
    <queryRefinementLevel>0</queryRefinementLevel>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>
    <description><![CDATA[Text Search for: wd40 or "wd repeat" or "wd 40"]]></description>
    <keywords><![CDATA[wd40 OR "wd repeat" OR "wd 40"]]></keywords>
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
    req = urllib2.Request(url,data=txt_query)
    response = urllib2.urlopen(req)
    result_pdb = response.read()
    pdbids = result_pdb.replace('\n',',')
    # get customed report
    return pdbids


@lt.run_time
def main():
    beta,chain_len,resolution = 16,160,3.0
    uniprot_pdbids,pdb_scores = rcsb_uniprot(beta,chain_len,resolution)
    scop_pdbids = rcsb_scop(beta,chain_len,resolution)
    pfam_pdbids = rcsb_pfam(beta,chain_len,resolution)
    txt_pdbids = rcsb_pfam(beta,chain_len,resolution)
    uniprot = set([u.split(':')[0] for u in uniprot_pdbids.split(',')[:-1]])
    scop = set([u.split(':')[0] for u in scop_pdbids.split(',')[:-1]])
    pfam = set([u.split(':')[0] for u in pfam_pdbids.split(',')[:-1]])
    txt = set([u.split(':')[0] for u in txt_pdbids.split(',')[:-1]])


    sns.set_color_codes('colorblind')
    venn2([uniprot,scop],['uniprot','scop'])
    plt.savefig('uniprot_scop.png',dpi=300)
    plt.close('all')
    venn2([uniprot,pfam],['uniprot','pfam'])
    plt.savefig('uniprot_pfam.png',dpi=300)
    plt.close('all')
    venn2([uniprot,txt],['uniprot','txt'])
    plt.savefig('uniprot_pfam.png',dpi=300)
    plt.close('all')
    venn2([uniprot,set.union(*[scop,pfam,txt])],['uniprot','pfam_scop_txt'])
    plt.savefig('uniprot__pfam_scop_txt.png',dpi=300)
    plt.close('all')
    venn2([pfam,txt],['pfam','txt'])
    plt.savefig('pfam_txt.png',dpi=300)
    plt.close('all')
    venn2([pfam,scop],['pfam','scop'])
    plt.savefig('pfam_scop.png',dpi=300)
    plt.close('all')
    venn3([pfam,scop,txt],['pfam','scop','txt'])
    plt.savefig('pfam_scop_txt.png',dpi=300)
    plt.close('all')

    lt.pickle_dump([uniprot,scop,pfam,txt],'search_method')
    write_lis_lis([uniprot,scop,pfam,txt],'rcsb_wd40_pdb',['uniprot','scop','pfam','txt'])

    f,ax = plt.subplots()
    total = set.union(*map(set,[uniprot,scop,pfam,txt]))
    sns.set_color_codes('pastel')
    methods = ['UniProt','Pfam','Txt','SCOP']
    wd = pd.DataFrame({'Search Method':methods,'Num':map(len,[total,total,total,total])})
    sns.barplot(x='Search Method',y='Num',data=wd,color='b')
    sns.set_color_codes('muted')
    wd = pd.DataFrame({'Search Method':methods,'Num':map(len,[uniprot,pfam,txt,scop])})
    sns.barplot(x='Search Method',y='Num',data=wd,color='b')
    ax.set(xlabel='Search Method',ylabel='Num',title='WD40 Structures in RCSB')
    plt.savefig('wd40_in_RCSB_pdbs',dpi=300)
    plt.close('all')

if __name__ == "__main__":
    main()







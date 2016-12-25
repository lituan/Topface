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
"""

import sys
import os
import urllib
import urllib2
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3,venn2
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
        'customReportColumns':'structureId,uniprotAcc,entityId,resolution,chainLength,releaseDate,pfamAccession',
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
        query = 'database:(type:pfam id:PF00400)'
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
    for i in range(10):
        try:
            print 'try uniprot ',key,'for the ',i+1,' time'
            data = urllib.urlencode(data)
            req = urllib2.Request(url,data)
            response = urllib2.urlopen(req)
            r = response.readlines()
            lines = set([line.rstrip('\r\n') for line in r])
        except:
            continue
        break
    print 'uniprot search ',key,' is finished'
    return lines

@lt.run_time
def rcsb_acc(accs):
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
    <description>Secondary structure has:  14 or more Beta Sheets</description>
    <polyStats.helixPercent.comparator>between</polyStats.helixPercent.comparator>
    <polyStats.helixCount.comparator>between</polyStats.helixCount.comparator>
    <polyStats.sheetPercent.comparator>between</polyStats.sheetPercent.comparator>
    <polyStats.sheetCount.comparator>between</polyStats.sheetCount.comparator>
    <polyStats.sheetCount.min>14</polyStats.sheetCount.min>
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>2</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.SequenceLengthQuery</queryType>
    <description>Sequence Length is between 150 and 10000 </description>
    <v_sequence.chainLength.min>150</v_sequence.chainLength.min>
    <v_sequence.chainLength.max>10000</v_sequence.chainLength.max>
    </orgPdbQuery>
    </queryRefinement>
    </orgPdbCompositeQuery>
    """
    for i in range(10):
        try:
            print 'try uniprot acc in rcsb'
            req = urllib2.Request(url,data=query)
            response = urllib2.urlopen(req)
            pdbids = response.read()
            pdbids = pdbids.replace('\n',',')
        except:
            continue
        break
    return pdbids

def get_wdsp_acc():
    with open('wdsp_uniprot_id.txt') as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        lines = [line.split()[0] for line in lines]
        return lines

@lt.run_time
def rcsb_uniprot():
    keywords = ['pfam','smart','supfam','uniprot_repeat','uniprot_keyword','prosite1','prosite2','prosite3']
    wd40s = []
    for key in keywords:
        for i in range(10):
            try:
                wd40s.append(uniprot_wd40(key))
            except:
                continue
            break
    wdsp = get_wdsp_acc()
    wd40s.append(wdsp)
    keywords.append('wdsp')

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

    uniprot_pdbids = rcsb_acc(total)
    report = rcsb_customreport(uniprot_pdbids)
    pdb_scores = [[] for i in range(9)]
    for p in report[1:]:
        for i,w in enumerate(wd40s_score):
            if p[2] in w:
                pdb_scores[i].append('_'.join(p[:3]))
    pdb_acc_scores = map(lambda x:x.split('_')[2],pdb_scores)
    pdb_acc_scores = map(set,pdb_acc_scores)


    write_lis_lis(pdb_scores,'uniprot_wd40_pdb_scores',[str(i) for i in range(1,10)])
    write_lis_lis(pdb_acc_scores,'uniprot_wd40_pdb_acc_scores',[str(i) for i in range(1,10)])

    print 'uniprot search is finished'
    return uniprot_pdbids,pdb_scores


@lt.run_time
def rcsb_pfam():
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
    <description>Secondary structure has:  14 or more Beta Sheets</description>
    <polyStats.helixPercent.comparator>between</polyStats.helixPercent.comparator>
    <polyStats.helixCount.comparator>between</polyStats.helixCount.comparator>
    <polyStats.sheetPercent.comparator>between</polyStats.sheetPercent.comparator>
    <polyStats.sheetCount.comparator>between</polyStats.sheetCount.comparator>
    <polyStats.sheetCount.min>14</polyStats.sheetCount.min>
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>2</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.SequenceLengthQuery</queryType>
    <description>Sequence Length is 150 and more</description>
    <v_sequence.chainLength.min>150</v_sequence.chainLength.min>
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>3</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ResolutionQuery</queryType>
    <description>Resolution is between 0.1 and 3.5 </description>
    <refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>
    <refine.ls_d_res_high.min>0.1</refine.ls_d_res_high.min>
    <refine.ls_d_res_high.max>3.5</refine.ls_d_res_high.max>
    </orgPdbQuery>
    </queryRefinement>
    </orgPdbCompositeQuery>
    """
    for i in range(10):
        try:
            print 'try pfam query in rcsb'
            req = urllib2.Request(url,data=pfam_query)
            response = urllib2.urlopen(req)
            result_pdb = response.read()
            pdbids = result_pdb.replace('\n',',')
        except:
            continue
        break
    # get customed report
    return pdbids

@lt.run_time
def rcsb_scop():
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
    <description>Secondary structure has:  14 or more Beta Sheets</description>
    <polyStats.helixPercent.comparator>between</polyStats.helixPercent.comparator>
    <polyStats.helixCount.comparator>between</polyStats.helixCount.comparator>
    <polyStats.sheetPercent.comparator>between</polyStats.sheetPercent.comparator>
    <polyStats.sheetCount.comparator>between</polyStats.sheetCount.comparator>
    <polyStats.sheetCount.min>14</polyStats.sheetCount.min>
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>2</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.SequenceLengthQuery</queryType>
    <description>Sequence Length is 150 and more</description>
    <v_sequence.chainLength.min>150</v_sequence.chainLength.min>
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>3</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ResolutionQuery</queryType>
    <description>Resolution is between 0.1 and 3.5 </description>
    <refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>
    <refine.ls_d_res_high.min>0.1</refine.ls_d_res_high.min>
    <refine.ls_d_res_high.max>3.5</refine.ls_d_res_high.max>
    </orgPdbQuery>
    </queryRefinement>
    </orgPdbCompositeQuery>
    """
    for i in range(10):
        try:
            print 'try scop query in rcsb'
            req = urllib2.Request(url,data=scop_query)
            response = urllib2.urlopen(req)
            result_pdb = response.read()
            pdbids = result_pdb.replace('\n',',')
        except:
            continue
        break
    return pdbids
    # get customed report

@lt.run_time
def rcsb_txt():

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
    <description>Secondary structure has:  14 or more Beta Sheets</description>
    <polyStats.helixPercent.comparator>between</polyStats.helixPercent.comparator>
    <polyStats.helixCount.comparator>between</polyStats.helixCount.comparator>
    <polyStats.sheetPercent.comparator>between</polyStats.sheetPercent.comparator>
    <polyStats.sheetCount.comparator>between</polyStats.sheetCount.comparator>
    <polyStats.sheetCount.min>14</polyStats.sheetCount.min>
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>2</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.SequenceLengthQuery</queryType>
    <description>Sequence Length is 150 and more</description>
    <v_sequence.chainLength.min>150</v_sequence.chainLength.min>
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>3</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ResolutionQuery</queryType>
    <description>Resolution is between 0.1 and 3.5 </description>
    <refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>
    <refine.ls_d_res_high.min>0.1</refine.ls_d_res_high.min>
    <refine.ls_d_res_high.max>3.5</refine.ls_d_res_high.max>
    </orgPdbQuery>
    </queryRefinement>
    </orgPdbCompositeQuery>
    """
    for i in range(10):
        try:
            print 'try txt query in rcsb'
            req = urllib2.Request(url,data=txt_query)
            response = urllib2.urlopen(req)
            result_pdb = response.read()
            pdbids = result_pdb.replace('\n',',')
        except:
            continue
        break
    # get customed report
    return pdbids


@lt.run_time
def main():
    uniprot_pdbids,pdb_scores = rcsb_uniprot()
    scop_pdbids = rcsb_scop()
    pfam_pdbids = rcsb_pfam()
    txt_pdbids = rcsb_pfam()
    uniprot = set([u.split(':')[0] for u in uniprot_pdbids.split(',')[0]])
    scop = set([u.split(':')[0] for u in scop_pdbids.split(',')[0]])
    pfam = set([u.split(':')[0] for u in pfam_pdbids.split(',')[0]])
    txt = set([u.split(':')[0] for u in txt_pdbids.split(',')[0]])


    venn2([uniprot,scop],['uniprot','scop'])
    plt.savefig('uniprot_scop.png','w')
    venn2([uniprot,pfam],['uniprot','pfam'])
    plt.savefig('uniprot_pfam.png','w')
    venn2([uniprot,txt],['uniprot','txt'])
    plt.savefig('uniprot_pfam.png','w')
    venn2([pfam,txt],['pfam','txt'])
    plt.savefig('pfam_txt.png','w')
    venn2([pfam,scop],['pfam','scop'])
    plt.savefig('pfam_scop.png','w')
    plt.close('all')

    write_lis_lis([uniprot,scop,pfam,txt],'rcsb_wd40',['uniprot','scop','pfam','txt'])

if __name__ == "__main__":
    main()







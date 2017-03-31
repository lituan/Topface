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
import cPickle as pickle
import lxml.etree as et
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
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

    report = rcsb_acc_customreport(total,beta,chain_len,resolution)
    pdb_scores = []
    for p in report:
        pdb_scores.append(p+[acc_score(p[2])])
    pdb_scores = sorted(pdb_scores,key=lambda x:x[-1],reverse=True)
    pickle.dump(pdb_scores,open('pdb_scores.pickle','w'))
    pdb_scores = pickle.load(open('pdb_scores.pickle'))
    with open('uniprot_pdb_scores.txt','w') as w_f:
        print >> w_f,'{0:<15}{1:<10}{2:<8}{3:<8}{4:<15}{5:<15}{6:<18}{7:<8}'.format('acc','pdb','chain','entity','resolution','chain_len','release','score')
        for p in pdb_scores:
            print >> w_f,'{0:<15}{1:<10}{2:<8}{3:<8}{4:<15}{5:<15}{6:<18}{7:<8}'.format(p[2],p[0],p[1],p[3],p[4],p[5],p[6],p[7])

    return pdb_scores


def rcsb_txt(beta,chain_len,resolution):

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
    pdbids = [p.split(':')[0] for p in pdbids.split(',') if p]
    return pdbids


def rcsb_scop(beta,chain_len,resolution):
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
    pdbids = [p.split(':')[0] for p in pdbids.split(',') if p]
    return pdbids


def rcsb_pfam(beta,chain_len,resolution):
    url = 'http://www.rcsb.org/pdb/rest/search'

    begin_query="""<orgPdbCompositeQuery version="1.0">
    <queryRefinement>
    <queryRefinementLevel>0</queryRefinementLevel>
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
    pfams = [pfam1,pfam2,pfam3,pfam4,pfam5]
    pdbids = ''
    for pfam in pfams:
        req = urllib2.Request(url,data=begin_query+pfam+end_query)
        response = urllib2.urlopen(req)
        result_pdb = response.read()
        pdbids += result_pdb.replace('\n',',')
    pdbids = [p.split(':')[0] for p in pdbids.split(',') if p]
    return pdbids


def rcsb_phos(pdbids):

    pdb_num = len(pdbids)
    pdbids = ', '.join(pdbids)
    query_begin = """
    <orgPdbCompositeQuery version="1.0">
    <queryRefinement>
    <queryRefinementLevel>0</queryRefinementLevel>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.StructureIdQuery</queryType>
    """
    query_pdb ="""
    <description>Simple query for a list of PDB IDs """ + '(' + str(pdb_num) + " IDs) : " + pdbids+"""</description>
    <structureIdList>""" + pdbids + """</structureIdList>
    """
    query_middle = """
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>1</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ChemCompIdQuery</queryType>
    """
    query_tpo = """
    <description>Chemical ID(s):  TPO and Polymeric type is Polymeric</description>
    <chemCompId>TPO</chemCompId>
    """
    query_sep = """
    <description>Chemical ID(s):  SEP and Polymeric type is Polymeric</description>
    <chemCompId>SEP</chemCompId>
    """
    query_ptr = """
    <description>Chemical ID(s):  PTR and Polymeric type is Polymeric</description>
    <chemCompId>PTR</chemCompId>
    """

    query_end = """
    <polymericType>Polymeric</polymericType>
    </orgPdbQuery>
    </queryRefinement>
    </orgPdbCompositeQuery>
    """

    query1 = query_begin + query_pdb + query_middle + query_tpo + query_end
    query2 = query_begin + query_pdb + query_middle + query_sep + query_end
    query3 = query_begin + query_pdb + query_middle + query_ptr + query_end

    phos_pdbs = []

    for query in [query1,query2,query3]:
        for i in range(20):
            try:
                url = 'http://www.rcsb.org/pdb/rest/search'
                req = urllib2.Request(url,data=query)
                response = urllib2.urlopen(req)
                pdbids = response.read()
                pdbids = pdbids.replace('\n',',')
                pdbids = pdbids.split(',')
                pdbids = [p for p in pdbids if p]
                phos_pdbs += pdbids
                break
            except e:
                print e
                continue
    return phos_pdbs

def check_pdb_status(pdbid):
    """Returns the status and up-to-date entry in the PDB for a given PDB ID"""
    url = 'http://www.rcsb.org/pdb/rest/idStatus?structureId=%s' % pdbid
    xmlf = urllib2.urlopen(url)
    xml = et.parse(xmlf)
    xmlf.close()
    status = None
    current_pdbid = pdbid
    for df in xml.xpath('//record'):
        status = df.attrib['status']  # Status of an entry can be either 'UNKWOWN', 'OBSOLETE', or 'CURRENT'
        if status == 'OBSOLETE':
            current_pdbid = df.attrib['replacedBy']  # Contains the up-to-date PDB ID for obsolete entries
    return [status, current_pdbid.lower()]

def check_current_directory(pdbid,directory):
    if not os.path.exists(directory):
        return False
    for root,dirs,files in os.walk(directory):
        for f in files:
            if f[-4:] == '.pdb' and len(f) == 8:
                if f[:4].lower() == pdbid.lower():
                    return True
    return False

def fetch_pdb(p):
    pdbid,pdbpath = p
    if check_current_directory(pdbid,pdbpath):
        return

    pdbid = pdbid.lower()
    print 'Checking status of PDB ID ' + pdbid
    state,current_entry = check_pdb_status(pdbid)
    if state == 'OBSOLETE':
        print 'entry is obsolete, getting ' + current_entry + ' instead'
    elif state == 'CURRENT':
        print 'entry is up to date'
    elif state == 'UNKNOWN':
        print 'Invalid PDB ID ' + pdbid
    print 'Downloading file from PDB...'
    pdburl = 'http://www.rcsb.org/pdb/files/'+current_entry+'.pdb'
    try:
        pdbfile = urllib2.urlopen(pdburl).read()
        if 'sorry' in pdbfile:
            print 'No file in PDB format available from wwwPDB for the given PDB ID ' + current_entry
        else:
            if not os.path.exists(pdbpath):
                os.makedirs(pdbpath)
            pdbpath = os.path.join(pdbpath,pdbid+'.pdb')
            with open(pdbpath,'w') as w_f:
                w_f.write(pdbfile)
                print pdbid+' is downloaded'
                return pdbpath
    except urllib2.HTTPError:
        print 'No file in PDB format available from wwwPDB for the given PDB ID ' + current_entry


def main():
    beta,chain_len,resolution = 20,200,3.5
    pdb_scores = rcsb_uniprot(beta,chain_len,resolution)
    pdbids = [p[0] for p in pdb_scores]
    phos_pdbs = rcsb_phos(pdbids)
    print 'num of uniprot search',len(phos_pdbs)

    phos_pdb_scores = [p for p in pdb_scores if p[0] in phos_pdbs]
    pickle.dump(phos_pdb_scores,open('wd40_phos_pdb_scores.pickle','w'))
    with open('wd40_phos_pdb_scores.txt','w') as w_f:
        print >> w_f,'{0:<15}{1:<10}{2:<8}{3:<8}{4:<15}{5:<15}{6:<18}{7:<8}'.format('acc','pdb','chain','entity','resolution','chain_len','release','score')
        for p in phos_pdb_scores:
            print >> w_f,'{0:<15}{1:<10}{2:<8}{3:<8}{4:<15}{5:<15}{6:<18}{7:<8}'.format(p[2],p[0],p[1],p[3],p[4],p[5],p[6],p[7])

    parameters = [(p,'uniprot_wd40_phos_pdb') for p in phos_pdbs]
    p = Pool(len(parameters))
    p.map(fetch_pdb,parameters)
    p.close()


    def other_special(other_phos_pdbs,other_name):
        other_phos_pdbs_special = set(other_phos_pdbs).difference(set(phos_pdbs))
        if other_phos_pdbs_special:
            with open(other_name+'_special_wd40_phos_pdbs.txt','w') as w_f:
                for p in other_phos_pdbs_special:
                    print >> w_f, p

            parameters = [(p,other_name+'_special_wd40_phos_pdb') for p in other_phos_pdbs_special]
            p = Pool(len(parameters))
            p.map(fetch_pdb,parameters)
            p.close()

    txt_pdbids = rcsb_txt(beta,chain_len,resolution)
    print 'txt',txt_pdbids
    txt_phos_pdbs = rcsb_phos(txt_pdbids)
    print 'txt search phos',txt_phos_pdbs
    other_special(txt_phos_pdbs,'txt')

    pfam_pdbids = rcsb_pfam(beta,chain_len,resolution)
    print 'pfam',pfam_pdbids
    pfam_phos_pdbs = rcsb_phos(pfam_pdbids)
    print 'pfam search phos',pfam_phos_pdbs
    other_special(pfam_phos_pdbs,'pfam')

    scop_pdbids = rcsb_scop(beta,chain_len,resolution)
    print 'scop',scop_pdbids
    scop_phos_pdbs = rcsb_phos(scop_pdbids)
    print 'scop search phos',scop_phos_pdbs
    other_special(scop_phos_pdbs,'scop')


    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111)
    total = set.union(*map(set,[phos_pdbs,scop_phos_pdbs,pfam_phos_pdbs,txt_phos_pdbs]))
    sns.set_color_codes('pastel')
    methods = ['UniProt','Pfam','Text','SCOP']
    wd = pd.DataFrame({'Methods':methods,'PDB Num':map(len,[total,total,total,total])})
    sns.barplot(x='Methods',y='PDB Num',data=wd,color='b')
    sns.set_color_codes('muted')
    wd = pd.DataFrame({'Methods':methods,'PDB Num':map(len,[phos_pdbs,pfam_phos_pdbs,txt_phos_pdbs,scop_phos_pdbs])})
    sns.barplot(x='Methods',y='PDB Num',data=wd,color='b')
    ax.set(xlabel='Methods',ylabel='PDB Num',title='Complex of WD40 Protein and Phosphorylated Protein in RCSB')
    plt.savefig('wd40_phos_in_RCSB_pdbs.png',dpi=300)
    plt.close('all')



if __name__ == "__main__":
    main()







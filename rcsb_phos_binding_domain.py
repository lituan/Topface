#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
use pfam accs and chemical ID to search structures Phos-binding
"""
import os
import sys
import urllib
import urllib2
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3,venn2
from multiprocessing import Pool


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

def rcsb_pfam(p):
    """
    use beta_sheet, chain length to make sure we get WD40s other than others
    accs be a list, ['Q969H0,P07834']
    return string format: '1A0R:1,1B9X:1,1B9Y:1,1C15:1'
    """
    url = 'http://www.rcsb.org/pdb/rest/search'

    query1 = """
    <orgPdbCompositeQuery version="1.0">
    <queryRefinement>
    <queryRefinementLevel>0</queryRefinementLevel>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.PfamIdQuery</queryType>
    """
    pfam1 = """
    <description>Pfam Accession Number PF00017 SH2 domain</description>
    <pfamID>PF00017</pfamID>
    """
    pfam2 = """
    <description>Pfam Accession Number PF00244 14-3-3 protein</description>
    <pfamID>PF00244</pfamID>
    """
    pfam3 = """
    <description>Pfam Accession Number PF00533 BRCA1 C Terminus (BRCT) domain</description>
    <pfamID>PF00533</pfamID>
    """
    pfam4 = """
    <description>Pfam Accession Number PF00498 FHA domain</description>
    <pfamID>PF00498</pfamID>
    """
    pfam5 = """
    <description>Pfam Accession Number PF03166 MH2 domain</description>
    <pfamID>PF03166</pfamID>
    """
    pfam6 = """
    <description>Pfam Accession Number PF00786 P21-Rho-binding domain</description>
    <pfamID>PF00786</pfamID>
    """
    pfam7 = """
    <description>Pfam Accession Number PF08416 Phosphotyrosine-binding domain</description>
    <pfamID>PF08416</pfamID>
    """
    pfam8 = """
    <description>Pfam Accession Number PF00397 WW domain</description>
    <pfamID>PF00397</pfamID>
    """
    pfam9 = """
    <description>Pfam Accession Number PF00168 C2 domain</description>
    <pfamID>PF00168</pfamID>
    """
    pfam10 = """
    <description>Pfam Accession Number PF00400 WD domain, G-beta repeat</description>
    <pfamID>PF00400</pfamID>
    """

    query2 = """
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>1</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ChemCompIdQuery</queryType>
    """
    phos1_any = """
    <description>Chemical ID(s):  PTR and Polymeric type is Any</description>
    <chemCompId>PTR</chemCompId>
    """
    phos2_any = """
    <description>Chemical ID(s):  TPO and Polymeric type is Any</description>
    <chemCompId>TPO</chemCompId>
    """
    phos3_any = """
    <description>Chemical ID(s):  SEP and Polymeric type is Any</description>
    <chemCompId>SEP</chemCompId>
    """
    phos1_poly = """
    <description>Chemical ID(s):  PTR and Polymeric type is Any</description>
    <chemCompId>PTR</chemCompId>
    """
    phos2_poly = """
    <description>Chemical ID(s):  TPO and Polymeric type is Any</description>
    <chemCompId>TPO</chemCompId>
    """
    phos3_poly = """
    <description>Chemical ID(s):  SEP and Polymeric type is Any</description>
    <chemCompId>SEP</chemCompId>
    """
    query3 = """
    <polymericType>Any</polymericType>
    </orgPdbQuery>
    </queryRefinement>
    </orgPdbCompositeQuery>
    """
    pfam_labels = ['SH2','14333','BRCT','FHA','MH2','PBD','PTB','WW','C2','WD40']
    phos_labels = ['PTR','TPO','SEP']
    pfams = [pfam1,pfam2,pfam3,pfam4,pfam5,pfam6,pfam7,pfam8,pfam9,pfam10]
    phoses_any = [phos1_any,phos2_any,phos3_any]
    phoses_poly = [phos1_poly,phos2_poly,phos3_poly]

    pfami,phosi = p
    query = query1 + pfams[pfami] + query2 + phoses_poly[phosi] + query3
    req = urllib2.Request(url,data=query)
    response = urllib2.urlopen(req)
    pdbids = response.read()
    pdbids = pdbids.replace('\n',',')
    return (pfam_labels[pfami],phos_labels[phosi],pdbids)


def main():

    parameters = [(pfi,phi) for pfi in range(10) for phi in range(3)]
    p = Pool(6)
    r = p.map(rcsb_pfam,parameters)
    p.close()

    for pf,ph,pd in r:
        if len(pd) > 0:
            print pf,ph,len(pd)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""
use different strategies to search WD40 structures in rcsb
all search has additional restrictions: resolution below 3.5, chain length longer than 150, beta strands more than 14

caution: use a low resolution as restriction will eliminate EM structures and NMR structures
"""

import sys
import os
import urllib
import urllib2
import httplib
import cPickle as pickle
from multiprocessing import Pool


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
    for i in range(20):
        try:
            print 'try time',i+1
            req = urllib2.Request(url,data)
            response = urllib2.urlopen(req)
            r = response.readlines()
            accs = set([line.rstrip('\r\n') for line in r])
            print 'uniprot search ',key,' is finished'
            break
        except httplib.IncompleteRead, e:
            continue
    return key,accs

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
    lines = [[line[0],line[1],line[2],line[3],line[3],float(line[4]),int(line[5]),line[6]] for line in lines[1:]]
    # line format: '3V7D', 'D', 'P07834', '2', 2.3', 464, '2012-05-02'
    return lines


def rcsb_uniprot(beta,chain_len,resolution):
    keywords = ['pfam','smart','supfam','prosite']
    p = Pool(4)
    result = p.map(uniprot_wd40,keywords)
    p.close()
    wd40s = []
    for k in keywords:
        for r,v in result:
            if r == k:
                wd40s.append(v)

    total = set.union(*map(set,wd40s))

    pdbids = rcsb_acc(total,beta,chain_len,resolution)
    report = rcsb_customreport(pdbids)
    pickle.dump(report,open('wd40s_rcsb_report.pickle','w'))

    chimera = [r for r in report if '#' in r[2]]

    fname = '_'.join(map(str,[beta,chain_len,resolution]))
    with open('chimera_'+fname+'.txt','w') as w_f:
        for c in chimera:
            print >> w_f,'{0:<8}{1:<6}{2:<20}{3:<6}{4:<8}{5:<10}{6:<}'.format(c[0],c[1],c[2],c[3],c[4],c[5],c[6])


def main():
    # beta,chain_len,resolution = 16,160,3.0
    # rcsb_uniprot(beta,chain_len,resolution)
    beta,chain_len,resolution = 20,200,3.0
    rcsb_uniprot(beta,chain_len,resolution)


if __name__ == "__main__":
    main()








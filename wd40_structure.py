#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
query WD40 structures using three methods
1. use pfam accession: pf00400 to search rcsb
2. use keyword: WD repeat AND database(type:pdb) to search uniprot
3. use wdsp_ids to query rcsb
then plot a venn chart
"""

import sys
import os
import urllib
import urllib2
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from bioservices.uniprot import UniProt


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
        'customReportColumns':'structureId,uniprotAcc,resolution,chainLength,releaseDate',
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

def rcsb_acc(accs):
    """
    accs be a list, ['Q969H0,P07834']
    return string format: '1A0R:1,1B9X:1,1B9Y:1,1C15:1'
    """
    url = 'http://www.rcsb.org/pdb/rest/search'
    accs = ','.join(accs)

    query_head = '<orgPdbQuery><queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>'
    query_body = '<description>Simple query for a list of Uniprot Accession IDs: '+accs+'</description><accessionIdList>'+accs+'</accessionIdList></orgPdbQuery>'
    req = urllib2.Request(url,data=query_head+query_body)
    response = urllib2.urlopen(req)
    pdbids = response.read()
    pdbids = pdbids.replace('\n',',')
    return pdbids

def rcsb_check_acc(accs,resolution_cutoff=5.0,seq_len_cutoff=200):
    """
    accs be a list,
    check if related pdb chains satisfies certain cutoff, like resoution, seq_len_cutoff
    """
    pdbids = rcsb_acc(accs)
    results = rcsb_customreport(pdbids)
    lines = [line for line in results[1:]]
    lines = [line for line in lines if line[3]] # select entries with resolution infomation
    lines = [line for line in lines if float(line[3]) < resolution_cutoff]  # select entries with resolution better than resolution_cutoff
    lines = [line for line in lines if float(line[4]) > seq_len_cutoff] # select entries with sequence length longer than seq_len_cutoff
    good_accs = set([line[2] for line in lines])
    return good_accs


def rcsb_pfam():
    """
    query rcsb using pfam: PF00400
    return a set
    """
    url = 'http://www.rcsb.org/pdb/rest/search'
    pfam_query = """
    <orgPdbQuery>
      <queryType>org.pdb.query.simple.PfamIdQuery</queryType>
      <description>Pfam Accession Number:  PF00400/description>
        <pfamID>PF00400</pfamID>
    </orgPdbQuery>
    """
    req = urllib2.Request(url,data=pfam_query)
    response = urllib2.urlopen(req)
    result_pdb = response.read()
    pdbids = result_pdb.replace('\n',',')
    # get customed report
    lines = rcsb_customreport(pdbids)
    rcsb_pfam_acc = set([line[2].strip('"') for line in lines])
    return rcsb_pfam_acc


def uniprot_wd40():
    """
    use keyword:WD repeat AND database(type:pdb) to query uniprot
    return a list of uniprot accesions
    """
    url = ' http://www.uniprot.org/uniprot/?'
    data ={
    'query':'keyword:"WD repeat" database:(type:pdb)',
    'columns':'id,entry name, length,organism,database(PDB)',
    'format':'tab',
    'compress':'no',
    'inclue':'no',
    }
    data = urllib.urlencode(data)
    req = urllib2.Request(url,data)
    response = urllib2.urlopen(req)
    r = response.readlines()
    lines = [line.rstrip('\r\n') for line in r[1:]]
    lines = [line for line in lines if line]
    lines = [line.split('\t') for line in lines]
    uniprot_wd40_acc = set([line[0] for line in lines])

    return uniprot_wd40_acc

def uniprot_map(query, f, t, format='tab'):
    """
    map one kind of id to another
    possible mapping, please check http://www.uniprot.org/help/programmatic_access
    return format: [['FBW1A_HUMAN', 'Q9Y297'], ['K3WSI0', 'K3WSI0']]
    """
    url = 'http://www.uniprot.org/mapping/'
    data = {
            'from':f,
            'to':t,
            'format':format,
            'query':query
            }
    data = urllib.urlencode(data)
    req = urllib2.Request(url,data)
    response = urllib2.urlopen(req)
    lines = response.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line for line in lines if line]
    lines = [line.split('\t') for line in lines[1:]]
    return lines

def rcsb_wdsp(wd_id_f):
    """
    map ids from wdsp_output to uniprot, then query rcsb, return a set of uniprot accesions
    """
    def read_wdsp_id(wd_id_f):
        with open(wd_id_f) as o_f:
            lines = o_f.readlines()
            lines = [line.rstrip('\r\n') for line in lines]
            lines = [line for line in lines if line]
            ids = [line.split()[0] for line in lines]
            return ids
    wdsp_ids = read_wdsp_id(wd_id_f)
    wdsp_ids = ','.join(wdsp_ids)
    mapped_ids = uniprot_map(wdsp_ids,'ACC+ID','ACC')
    wdsp_acc = [mi[1] for mi in mapped_ids]
    pdbids = rcsb_acc(wdsp_acc)

    lines = rcsb_customreport(pdbids)
    rcsb_wdsp_acc = set([line[2].strip('"') for line in lines])
    return rcsb_wdsp_acc

def rcsb_wdsp(wdsp_acc_f):
    """
    read a file contian mapped uniprot accesions, mapping may take a long time
    """
    def read_wdsp_acc(wdsp_acc_f):
        with open(wdsp_acc_f) as o_f:
            lines = o_f.readlines()
            lines = [line.rstrip('\r\n') for line in lines]
            lines = [line for line in lines if line]
            lines = [line.split()[0] for line in lines]
            return lines
    wdsp_acc = read_wdsp_acc(wdsp_acc_f)
    pdbids = rcsb_acc(wdsp_acc)

    lines = rcsb_customreport(pdbids)
    rcsb_wdsp_acc = set([line[2].strip('"') for line in lines])
    return rcsb_wdsp_acc

def plot_venn(set1,set2,set3,labels,filename):
    set1 = set(set1)
    set2 = set(set2)
    set3 = set(set3)
    venn3([set1,set2,set3],labels)
    plt.savefig('test.png',dpi=300)


def main():
    rcsb_pfam_acc = rcsb_pfam()
    uniprot_wd_acc = uniprot_wd40()
    rcsb_wdsp_acc = rcsb_wdsp('wd648_id.txt')

    total_accs = set.union(rcsb_pfam_acc,uniprot_wd_acc,rcsb_wdsp_acc)
    good_accs = rcsb_check_acc(total_accs,resolution_cutoff=5.0,seq_len_cutoff=200)

    rcsb_pfam_acc = rcsb_pfam_acc.intersection(good_accs)
    uniprot_wd_acc = uniprot_wd_acc.intersection(good_accs)
    rcsb_wdsp_acc = rcsb_wdsp_acc.intersection(good_accs)

    plot_wenn(rcsb_pfam_acc,uniprot_acc,rcsb_wdsp_acc)







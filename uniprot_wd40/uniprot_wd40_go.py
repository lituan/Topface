#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
use UniProt REST service to search sequences annotated as 'WD40'
"""

import os
import sys
import urllib
import urllib2
import cPickle as pickle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from multiprocessing import Pool
import time

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
    req = urllib2.Request(url,data)
    response = urllib2.urlopen(req)
    r = response.readlines()
    lines = set([line.rstrip('\r\n') for line in r])

    print 'uniprot search',key,'is finished'
    return key,lines


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


def uniprot_go():
    query0 = 'goa:("RNA polymerase II C-terminal domain phosphoserine binding [1990269]")' # RNA polymerase domain phosphoserine binding
    query1 = 'goa:("phosphoserine binding [50815]")' # phosphoserine binding
    query2 = 'goa:("phosphothreonine binding [50816]")' # phosphothreonine binding
    query3 = 'goa:("phosphotyrosine binding [1784]")' # phosphotyrosine binding
    query4 = 'goa:("protein phosphorylated amino acid binding [45309]")' # protein phosphorylated amino acid binding
    query5 = 'goa:("phosphoprotein binding [51219]")' # phosphoprotein binding

    accs = []
    for query in [query0,query1,query2,query3,query4,query5]:

        url = ' http://www.uniprot.org/uniprot/?'
        data ={
        'query':query,
        'format':'list',
        }

        for i in range(10):
            try:
                data = urllib.urlencode(data)
                req = urllib2.Request(url,data)
                response = urllib2.urlopen(req)
                result = response.readlines()
                lines = set([line.rstrip('\r\n') for line in result])
                accs += lines
                break
            except Exception,e:
                continue
    return accs


def uniprot_acc(acc):
    query = 'accession:'+acc
    url = ' http://www.uniprot.org/uniprot/?'
    data ={
    'query':query,
    'columns':'entry name,go',
    'format':'tab',
    }
    for i in range(10):
        try:
            data = urllib.urlencode(data)
            req = urllib2.Request(url,data)
            response = urllib2.urlopen(req)
            result = response.readlines()

            return result[1].split('\t')
        except Exception,e:
            continue


def main():
    begin = time.time()

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
    # wd40s = pickle.load(open('wd40s.pickle'))

    # total = set.union(*map(set,wd40s))


    # phos_binding_accs = uniprot_go()

    # phos_binding_wd40s = set(total).intersection(set(phos_binding_accs))

    # pickle.dump(phos_binding_wd40s,open('phos_binding_wd40s.pickle','w'))
    phos_binding_wd40s = pickle.load(open('phos_binding_wd40s.pickle'))

    p = Pool(4)
    result = p.map(uniprot_acc,phos_binding_wd40s)
    p.close()

    phos_go = ['0050815','0050816','0001784','0045309','0051219','1990269']
    with open('phos_binding_wd40s.txt','w') as w_f:
        for r in result:
            entry_name = r[0]
            go_terms = r[1].split(';')
            phos_go_terms = [g for g in go_terms if g.split('GO:')[-1].split(']')[0] in phos_go]
            print >> w_f,entry_name,phos_go_terms

    end = time.time()
    print end-begin


if __name__ == "__main__":
    main()

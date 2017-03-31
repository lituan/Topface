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


def uniprot_ec_acc(acc):

    query = 'accession:'+acc
    url = ' http://www.uniprot.org/uniprot/?'
    data ={
    'query':query,
    'columns':'entry name,ec',
    'format':'tab',
    }
    for i in range(10):
        try:
            data = urllib.urlencode(data)
            req = urllib2.Request(url,data)
            response = urllib2.urlopen(req)
            result = response.readlines()

            print acc,result
            return result[1].split('\t')
        except Exception,e:
            continue


def uniprot_ec():
    query = 'database:(type:pfam id:PF00400) or database:(type:pfam id:PF12894) or database:(type:pfam id:PF16529) or database:(type:pfam id:PF16756) or database:(type:pfam id:PF17005) \
             or database:(type:smart id:SM00320) \
             or database:(type:supfam id:SSF50978) \
             or database:(type:prosite id:PS00678) or database:(type:prosite id:PS50082) or database:(type:prosite id:PS50294) \
             or keyword:"WD repeat" or annotation:(type:repeat wd) or family:"WD repeat"'

    url = ' http://www.uniprot.org/uniprot/?'
    data ={
    'query':query,
    'columns':'entry name,ec',
    'format':'tab',
    }
    for i in range(10):
        try:
            data = urllib.urlencode(data)
            req = urllib2.Request(url,data)
            response = urllib2.urlopen(req)
            result = response.readlines()

            result = [line.split('\t') for line in result[1:]]
            return result
        except Exception,e:
            continue

def main():
    begin = time.time()

    # keywords = ['pfam','smart','supfam','prosite','uniprot']
    # p = Pool(5)
    # result = p.map(uniprot_wd40,keywords)
    # p.close()

    # total_uniprot_wd40s = set.union(*[r[1] for r in result])
    # pickle.dump(total_uniprot_wd40s,open('uniprot_wd40.pickle','w'))
    # total_uniprot_wd40s = pickle.load(open('uniprot_wd40.pickle'))

    # wdsp = get_wdsp_acc()
    # wdsp_special = set(wdsp).difference(set(total_uniprot_wd40s))
    # pickle.dump(wdsp_special,open('wdsp_special.pickle','w'))
   # wdsp_special = pickle.load(open('wdsp_special.pickle'))

    # p = Pool(4)
    # wdsp_special_ecs = p.map(uniprot_ec_acc,list(wdsp_special))
    # p.close()

    # pickle.dump(wdsp_special_ecs,open('wdsp_special_ecs.pickle','w'))

    # uniprot_ecs = uniprot_ec()
    # pickle.dump(uniprot_ecs,open('uniprot_ecs.pickle','w'))

    wdsp_special_ecs = pickle.load(open('wdsp_special_ecs.pickle'))
    uniprot_ecs = pickle.load(open('uniprot_ecs.pickle'))

    with open('wdsp_special_ecs.txt','w') as w_f:
        for w in wdsp_special_ecs:
            print >> w_f,w

    with open('uniprot_ecs.txt','w') as w_f:
        for w in uniprot_ecs:
            print >> w_f,w

    wd40_ec = [w for w in uniprot_ecs if len(w[1]) > 4] + [w for w in wdsp_special_ecs if len(w[1]) > 4]

    with open('wd40_ec.txt','w') as w_f:
        for w in wd40_ec:
            print >> w_f, '{0:<30}{1:<}'.format(w[0],w[1])

    end = time.time()
    print end-begin


if __name__ == "__main__":
    main()

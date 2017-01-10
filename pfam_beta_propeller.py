#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
query propeller annotation in uniprot accordint to pfam
"""
import sys
import os
import urllib
import urllib2
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3,venn2
from collections import OrderedDict
import cPickle as pickle

pfam_propellers = {
'WD40'            : 'PF00400',
'WD40_3'          : 'PF15911',
'WD40_4'          : 'PF16300',
'ANAPC4_WD40'     : 'PF12894',
'Ge1_WD40'        : 'PF16529',
'PALB2_WD40'      : 'PF16756',
'WD40_like'       : 'PF17005',
'Arylesterase'    : 'PF01731',
'Arylsulfotran_2' : 'PF14269',
'Arylsulfotrans'  : 'PF05935',
'BBS2_Mid'        : 'PF14783',
'Beta_propel'     : 'PF09826',
'Coatomer_WDAD'   : 'PF04053',
'CPSF_A'          : 'PF03178',
'Cytochrom_D1'    : 'PF02239',
'DPPIV_N'         : 'PF00930',
'DUF1513'         : 'PF07433',
'DUF1668'         : 'PF07893',
'DUF2415'         : 'PF10313',
'DUF4221'         : 'PF13970',
'DUF4934'         : 'PF16288',
'DUF5046'         : 'PF16465',
'DUF5050'         : 'PF16472',
'DUF5122'         : 'PF17164',
'DUF5128'         : 'PF17170',
'DUF839'          : 'PF05787',
'eIF2A'           : 'PF08662',
'FG-GAP'          : 'PF01839',
'FG-GAP_2'        : 'PF14312',
'Frtz'            : 'PF11768',
'Glu_cyclase_2'   : 'PF05096',
'Gmad1'           : 'PF10647',
'GSDH'            : 'PF07995',
'IKI3'            : 'PF04762',
'Itfg2'           : 'PF15907',
'Kelch_1'         : 'PF01344',
'Kelch_2'         : 'PF07646',
'Kelch_3'         : 'PF13415',
'Kelch_4'         : 'PF13418',
'Kelch_5'         : 'PF13854',
'Kelch_6'         : 'PF13964',
'Lactonase'       : 'PF10282',
'Ldl_recept_b'    : 'PF00058',
'Lgl_C'           : 'PF08596',
'LVIVD'           : 'PF08309',
'Me-amine-dh_H'   : 'PF06433',
'MRJP'            : 'PF03022',
'Nbas_N'          : 'PF15492',
'Neisseria_PilC'  : 'PF05567',
'NHL'             : 'PF01436',
'Nucleoporin_N'   : 'PF08801',
'Nup160'          : 'PF11715',
'PD40'            : 'PF07676',
'Pectate_lyase22' : 'PF14583',
'Peptidase_S9_N'  : 'PF02897',
'Phytase-like'    : 'PF13449',
'PQQ'             : 'PF01011',
'PQQ_2'           : 'PF13360',
'PQQ_3'           : 'PF13570',
'RAG2'            : 'PF03089',
'RCC1'            : 'PF00415',
'RCC1_2'          : 'PF13540',
'Reg_prop'        : 'PF07494',
'SBBP'            : 'PF06739',
'SBP56'           : 'PF05694',
'SdiA-regulated'  : 'PF06977',
'SGL'             : 'PF08450',
'Str_synth'       : 'PF03088',
'TcdB_toxin_midN' : 'PF12256',
'TolB_like'       : 'PF15869',
'VCBS'            : 'PF13517'
}

def uniprot_query(key='wd40',pfam_id='PF00400',pdb=False):
    """
    use annotations from different database to query WD40 in uniprot
    return a list of uniprot accesions
    """
    query = 'database:(type:pfam id:'+pfam_id+')'
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
    print 'uniprot search ',key,' is finished'
    return lines

def query(pfam):
    key,pfam_id = pfam[0],pfam[1]
    for i in range(10):
        try:
            ids = uniprot_query(key,pfam_id)
            propellers = [key,ids]
            print 'query ',key,'is successful'
            print len(ids),' entries got'
            break
        except:
            continue
    return propellers

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

def main():
    pfams = [[key,pfam_id] for key,pfam_id in pfam_propellers.iteritems() ]
    from multiprocessing import Pool
    p = Pool(71)
    propellers = p.map(query,pfams)
    p.close()
    keys = [k for k,v in propellers]
    values = [v for k,v in propellers]
    write_lis_lis(values,'pfam_propellers',cols=keys)

    pickle.dump(propellers,open('pfam_propellers.pickle','w'))

    # pfams = pickle.load(open('pfam_propellers.pickle','r'))
    propellers = OrderedDict([[k,v] for k,v in propellers if len(v) > 0]) # filter empty entries
    print 'types of propellers: ',len(propellers)


    #plot all propellers
    f,ax = plt.subplots(figsize=(10,8))
    sns.set_color_codes('pastel')
    wd = pd.DataFrame({'Propellers':propellers.keys(),'Num':map(len,propellers.values())})
    wd = wd.sort_values('Num',ascending=True)
    sns.barplot(x='Propellers',y='Num',data=wd,color='b')
    ax.set(xlabel='Propellers',ylabel='Num',title='Num of Different Propellers in Uniprot')
    plt.xticks(rotation='vertical')
    plt.savefig('num_of_different_propellers_in_uniprot',dpi=300)
    plt.close('all')

    #plot wd40s
    f,ax = plt.subplots()
    wd40_names = ['WD40','WD40_3','WD40_4','ANAPC4_WD40','Ge1_WD40','PALB2_WD40']
    wd40_names = [k for k in wd40_names if k in propellers.keys()]
    wd40s_large = [propellers[k] for k in wd40_names]
    wd40s_total = set.union(*map(set,wd40s_large))
    wd = pd.DataFrame({'WD40_entries':wd40_names,'Num':map(len,[wd40s_large for i in range(len(wd40_names))])})
    sns.barplot(x='WD40_entries',y='Num',data=wd,color='b')
    sns.set_color_codes('muted')
    wd = pd.DataFrame({'WD40_entries':wd40_names,'Num':map(len,wd40s_large)})
    sns.barplot(x='WD40_entries',y='Num',data=wd,color='b')
    ax.set(xlabel='WD40 entries in Pfam',ylabel='Num',title='WD40 in UniProt (Pfam annnotation)')
    plt.savefig('wd40_in_uniprot_by_pfam',dpi=300)
    plt.close('all')

    # plot wd40s and non-wd40s

    f,ax = plt.subplots()
    all_propellers = set.union(*map(set,[v for k,v in propellers.iteritems()]))
    wd40s_small =  propellers['WD40']
    non_wd40s = all_propellers.difference(wd40s_total)
    wd = pd.DataFrame({'Propellers':['WD40','non-WD40'],'Num':map(len,[wd40s_total,non_wd40s])})
    sns.barplot(x='Propellers',y='Num',data=wd,color='b')
    ax.set(xlabel='Propellers',ylabel='Num',title='WD40 vs non-WD40')
    plt.savefig('WD40_vs_nonWD40',dpi=300)
    plt.close('all')

if __name__ == "__main__":
    main()


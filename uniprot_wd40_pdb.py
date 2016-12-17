#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import urllib
import urllib2

import lt

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
    data = urllib.urlencode(data)
    req = urllib2.Request(url,data)
    response = urllib2.urlopen(req)
    r = response.readlines()
    lines = set([line.rstrip('\r\n') for line in r])

    return lines

def uniprot_txt_parser(uniprot_txt_lines):
    """
    parser protein family annotation and pdb annotation
    return data structure:
    uniprot['P06738'] = {
    'interpro': [('IPR001680',1)]
    'pfam': [('PF00400',7)]
    'smart':[('SM00230',3)]
    'pdb':[('1nex','a',263-744),('2p63',222-273)]
    ...
    }
    """
    uniprot = {}
    entry_line = [i for i,l in enumerate(uniprot_txt_lines) if l[:2]=='ID']
    entry_line.append(len(uniprot_txt_lines))
    begin_end = [(begin,entry_line[i+1]) for i,begin in enumerate(entry_line[:-1])]
    for begin,end in begin_end:
        acc = uniprot_txt_lines[begin+1].replace(';',' ').split()[1]
        uniprot[acc] = {}
        for line in uniprot_txt_lines[begin:end]:
            line = line.rstrip('\r\n')
            line = line.rstrip('.')
            line = line.replace(';',' ')
            words = line.split()
            if words[0] == 'DR' and words[1] =='InterPro':
                if uniprot[acc].has_key('interpro'):
                    uniprot[acc]['interpro'].append((words[2],1))
                else:
                    uniprot[acc]['interpro'] = [(words[2],1)]
            elif words[0] == 'DR' and words[1] == 'Pfam':
                if uniprot[acc].has_key('pfam'):
                    uniprot[acc]['pfam'].append((words[2],int(words[-1])))
                else:
                    uniprot[acc]['pfam'] = [(words[2],int(words[-1]))]
            elif words[0] == 'DR' and words[1] == 'SMART':
                if uniprot[acc].has_key('smart'):
                    uniprot[acc]['smart'].append((words[2],words[-1]))
                else:
                    uniprot[acc]['smart'] = [(words[2],words[-1])]
            elif words[0] == 'DR' and words[1] == 'SUPFAM':
                if uniprot[acc].has_key('supfam'):
                    uniprot[acc]['supfam'].append((words[2],words[-1]))
                else:
                    uniprot[acc]['supfam'] = [(words[2],words[-1])]
            elif words[0] == 'DR' and words[1] == 'PROSITE':
                if uniprot[acc].has_key('prosite'):
                    uniprot[acc]['prosite'].append((words[2],words[-1]))
                else:
                    uniprot[acc]['prosite'] = [(words[2],words[-1])]
            elif words[0] == 'DR' and words[1] =='PDB':
                w = words[-1].replace('/',' ')
                w = w.replace('=',' ')
                w = w.replace('-',' ')
                w = w.split()
                w = words[2:-1]+w

                if uniprot[acc].has_key('pdb'):
                    uniprot[acc]['pdb'].append(w)
                else:
                    uniprot[acc]['pdb'] = [w]

    return uniprot


@lt.run_time
def uniprot_accs(accs):
    results = []
    for acc in accs:
        url = ' http://www.uniprot.org/uniprot/?'
        data ={
        'query':'id: '+ acc,
        'format':'txt',
        'compress':'no',
        'include':'no',
        }
        data = urllib.urlencode(data)
        req = urllib2.Request(url,data)
        for i in range(10):
            try:
                response = urllib2.urlopen(req)
                r = response.readlines()
            except:
                continue
            results += r
            break

    return results


@lt.run_time
def align_lis_lis(lis_lis):
    """align and trans nested list to print a table"""
    lis_lis = [[str(l) for l in lis]
               for lis in lis_lis]  # trans every element to str
    #make all inner lists of the same length
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [lis + (inner_lis_max_len - len(lis)) * [''] for lis in lis_lis]
    #trans list, so that the elements of the same column are in one list
    lis_lis = [[lis[i] for lis in lis_lis] for i in range(inner_lis_max_len)]
    #make element in the same list have the same length
    aligned = []
    for lis in lis_lis:
        width = max([len(l) for l in lis])
        lis = [l + (width - len(l)) * ' ' for l in lis]
        aligned.append(lis)
    #trans list_list to the original list_list
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [[lis[i] for lis in aligned] for i in range(inner_lis_max_len)]
    return lis_lis

@lt.run_time
def write_lis_lis(lis_lis,filename,cols=[]):
    new_lis_lis = align_lis_lis(lis_lis)
    new_lis_lis = ['\t;'.join([lis_lis[i][j] for i in range(new_lis_lis)]) for j in range(len(new_lis_lis[0]))]
    with open(filename+'.txt','w') as w_f:
        if cols:
            print >> w_f,'\t;'.join(cols)
        for l in new_lis_lis:
            print >> w_f,l

def main():
    keywords = ['pfam','smart','supfam','interpro_repeat','interpro_domain','uniprot_repeat','uniprot_keyword','prosite1','prosite2','prosite3']
    # keywords= ['pfam','smart','supfam','prosite1','prosite2','prosite3']
    wd40s = []
    for key in keywords:
        for i in range(10):
            try:
                wd40s.append(uniprot_wd40(key,pdb=True))
            except:
                continue
            break

    total = set.union(*map(set,wd40s))
    total_txt= uniprot_accs(total)

    # i = 10
    # while i > 0 :
        # try:
            # total_txt= uniprot_accs(total)
        # except:
            # continue
        # i -= 1
        # break

    total_annotation = uniprot_txt_parser(total_txt)
    # select pdbs with length > 200
    wd40_pdbs = []
    for acc,anno in total_annotation.iteritems():
        good_pdbs = []
        for p in anno['pdb']:
            try:
                if int(p[-1]) - int(p[-2]) > 150:
                    good_pdbs.append(p)
            except:
                print p
        if good_pdbs:
            wd40_pdbs.append((acc,good_pdbs))

    with open('uniprot_good_pdbs.txt','w') as w_f:
        for acc,good_pdbs in wd40_pdbs:
            for good_pdb in good_pdbs:
                print >> w_f,'\t;'.join([acc]+good_pdb)

    wd40s_good = [w[0] for w in wd40_pdbs]

    wd40s = map(lambda x: [xi for xi in x if xi in wd40s_good],wd40s)

    total = set.union(*map(set,wd40s))
    total_repeat = []
    for w in wd40s:
        total_repeat += w
    # if an entry apears in n different querys, its score is n
    wd40s_score = [[] for i in range(10)]
    for i in total:
        num = total_repeat.count(i)
        wd40s_score[num-1].append(i)
    write_lis_lis(wd40s_score,'uniprot_wd40_pdb_score',[str(i) for i in range(1,11)])

    write_lis_lis(wd40s,'uniprot_wd40_pdb',keywords)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import urllib
import urllib2

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
    'compress':'no',
    'inclue':'no',
    }
    data = urllib.urlencode(data)
    req = urllib2.Request(url,data)
    response = urllib2.urlopen(req)
    r = response.readlines()
    lines = set([line.rstrip('\r\n') for line in r])
    return lines

def uniprot_accs(accs):
    results = ''
    for acc in accs:
        url = ' http://www.uniprot.org/uniprot/?'
        data ={
        'query':'id: '+ acc,
        'columns':columns,
        'format':'txt',
        'compress':'no',
        'include':'no',
        }
        data = urllib.urlencode(data)
        req = urllib2.Request(url,data)
        response = urllib2.urlopen(req)
        r = response.readlines()
        results += r

    return results

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
        for line in uniprot_txt_lines[begin:end]:
            line = line.rstrip('\r\n')
            line = line.rstrip('.')
            line = line.replace(';',' ')
            words = line.split()
            if words[0] == 'AC':
                acc = words[1]
                uniprot[acc] = {}
            elif words[0] == 'DR' and words[1] =='InterPro':
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
            # elif words[0] == 'DR' and words[1] =='PDB':
                # w = words[-1].replace('/',' ')
                # w = w.replace('=',' ')
                # w = w.replace('-',' ')
                # w = w.split()
                # w = words[2:-1]+w

                # if uniprot[acc].has_key('pdb'):
                    # uniprot[acc]['pdb'].append(w)
                # else:
                    # uniprot[acc]['pdb'] = [w]

    return uniprot


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


def main():
    # keywords = ['pfam','smart','supfam','interpro_repeat','interpro_domain','uniprot_repeat','uniprot_keyword','prosite1','prosite2','prosite3']
    keywords= ['pfam','smart','supfam','prosite1','prosite2','prosite3']
    wd40s = []
    for key in keywords:
        for i in range(10):
            try:
                wd40s.append(uniprot_wd40(key))
            except:
                continue
            break

    total = set.union(*map(set,wd40s))
    while True:
        try:
            total_txt= uniprot_accs(total)
        except:
            continue
        break
    total_annotation = uniprot_txt_parser(total_txt)

    keywords= ['pfam','smart','supfam','prosite1','prosite2','prosite3']
    keyids = ['PF00400','SM00320','SSF50978','PS00678','PS50082','PS50294']
    wd40s = [[] for i i range(len(keywords))]
    # select accs with hits > 6
    for acc,anno in total_annotation.iteritems():
        for i,key in enumerate(kewwords):
            for a in anno[key]:
                if a[0] ==  keyids[i] and  a[1] > 6:


    wd40s = align_lis_lis(wd40s)
    longest = max(map(len,wd40s))

    with open('uniprot_wd40.txt','w') as w_f:
        print >> w_f, '{0:<20}{1:<20}{2:<20}{3:<20}{4:<20}{5:<20}'.format(keywords[0],keywords[1],keywords[2],keywords[3],keywords[4],keywords[5])
        for i in range(longest):
            print >> w_f, '{0:<20}{1:<20}{2:<20}{3:<20}{4:<20}{5:<20}'.format(wd40s[0][i],wd40s[1][i],wd40s[2][i],wd40s[3][i],wd40s[4][i],wd40s[5][i])

    total = set.union(*map(set,wd40s))
    total_repeat = []
    for w in wd40s:
        total_repeat += w
    # if an entry apears in n different querys, its score is n
    wd40s_score = [[] for i in range(6)]
    for i in total:
        num = total_repeat.count(i)
        wd40s_score[num-1].append(i)
    scores = [str(i) for i in range(1,7)]
    with open('uniprot_wd40_score.txt','w') as w_f:
        print >> w_f, '{0:<20}{1:<20}{2:<20}{3:<20}{4:<20}{5:<20}'.format(scores[0],scores[1],scores[2],scores[3],scores[4],scores[5])
        for i in range(longest):
            print >> w_f, '{0:<20}{1:<20}{2:<20}{3:<20}{4:<20}{5:<20}'.format(wd40s_score[0][i],wd40s_score[1][i],wd40s_score[2][i],wd40s_score[3][i],wd40s_score[4][i],wd40s_score[5][i])

if __name__ == "__main__":
    main()

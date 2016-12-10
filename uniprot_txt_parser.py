#!/usr/bin/env python
# -*- coding: utf-8 -*-


def uniprot_txt_parser(uniprot_txt_lines):
    """
    parser protein family annotation and pdb annotation
    return data structure:
    uniprot['P06738'] = {
    'interpro': [('IPR001680',1),]
    'pfam': [('PF00400',7),]
    'smart':[('SM00230',3),]
    'pdb:[['3V7D', 'X-ray', '2.31', 'A', 'B', 'D', '263', '744'],]
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


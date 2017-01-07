#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3
from collections import OrderedDict
import cPickle as pickle


pfams = pickle.load(open('pfam_propellers.pickle','r'))
pfams = OrderedDict([[k,v] for k,v in pfams if len(v) > 0]) # filter empty entries
print 'types of propellers: ',len(pfams)


#plot all propellers
f,ax = plt.subplots(figsize=(10,8))
sns.set_color_codes('pastel')
wd = pd.DataFrame({'Propellers':pfams.keys(),'Num':map(len,pfams.values())})
wd = wd.sort_values('Num',ascending=True)
sns.barplot(x='Propellers',y='Num',data=wd,color='b')
ax.set(xlabel='Propellers',ylabel='Num',title='Num of Different Propellers in Uniprot')
plt.xticks(rotation='vertical')
plt.savefig('num_of_different_propellers_in_uniprot',dpi=300)
plt.close('all')

#plot wd40s
f,ax = plt.subplots()
wd40_names = ['WD40','WD40_3','WD40_4','ANAPC4_WD40','Ge1_WD40','PALB2_WD40']
wd40_names = [k for k in wd40_names if k in pfams.keys()]
wd40s_large = [pfams[k] for k in wd40_names]
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
all_propellers = set.union(*map(set,[v for k,v in pfams.iteritems()]))
wd40s_small =  pfams['WD40']
non_wd40s = all_propellers.difference(wd40s_total)
wd = pd.DataFrame({'Propellers':['WD40','non-WD40'],'Num':map(len,[wd40s_total,non_wd40s])})
sns.barplot(x='Propellers',y='Num',data=wd,color='b')
ax.set(xlabel='Propellers',ylabel='Num',title='WD40 vs non-WD40')
plt.savefig('WD40_vs_nonWD40',dpi=300)
plt.close('all')

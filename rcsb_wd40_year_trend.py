#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""
plot wd40 structures per year
"""

import sys
import os
import urllib
import urllib2
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3,venn2
from multiprocessing import Pool
from datetime import datetime
import lt

with open(sys.argv[-1]) as o_f:
    lines = o_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line for line in lines if line]
    lines = [line.split() for line in lines]
    pdb_scores = lines[1:]

    begin_year = min([int(str(p[6]).split('-')[0]) for p in pdb_scores])
    end_year = datetime.now().year
    num = end_year-begin_year+1
    years = [[] for i in range(num)]
    for p in pdb_scores:
        acc = p[0]
        year = int(str(p[6]).split('-')[0])
        years[year-begin_year].append(acc)
    wd = pd.DataFrame({'Year':range(begin_year,end_year+1),'Num of WD40':map(len,years)})
    f,ax = plt.subplots()
    sns.set_color_codes('pastel')
    sns.barplot(x='Year',y='Num of WD40',data=wd,color='b')
    ax.set(xlabel='Year',ylabel='Num of WD40',title='WD40 Structures per Year')
    plt.xticks(rotation=90)
    plt.savefig('WD40 Structures per Year',dpi=300)
    plt.close('all')


















#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
use some code from https://github.com/nvictus/svgpath2mpl/blob/master/examples/seqlogo.ipynb
plot sequence logo using glyph
usage: python seqlogo.py *.msa (fasta format)
"""
import sys
import os
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib import ticker
import seaborn.apionly as sns
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

from svgpath2mpl import parse_path
from xml.dom import minidom


def read_glyph():
    glyphs = {}
    for i in range(65,91):
        svg_file_name = 'svg/'+str(i)+'.svg'
        svg = minidom.parse(svg_file_name)
        glyphs[chr(i)]= [path.getAttribute('d') for path in svg.getElementsByTagName('path')][-1]
    return glyphs

def get_glyph_patch(path_data, color, x, y, dx, dy, **kwargs):
    kwargs.setdefault('facecolor', color)
    kwargs.setdefault('edgecolor', 'none')
    path = parse_path(path_data)
    # normalize and flip upside down
    path.vertices[:, 0] -= path.vertices[:, 0].min()
    path.vertices[:, 1] -= path.vertices[:, 1].min()
    path.vertices[:, 0] /= path.vertices[:, 0].max()
    path.vertices[:, 1] /= path.vertices[:, 1].max()
    path.vertices[:, 1] = 1 - path.vertices[:, 1]
    # scale then translate
    path.vertices *= [dx, dy]
    path.vertices += [x, y]
    return PathPatch(path, **kwargs)

def draw_logo(ax,matrix,charwidth):

    # R K H is shown blue, Y S T is show red
    colors = {'A':'black','C':'green','D':'red','E':'red','F':'black','G':'green','H':'blue','I':'black','K':'blue','L':'black','M':'black','N':'green','P':'black','Q':'green','R':'blue','S':'green','T':'green','V':'black','W':'black','Y':'green'}
    glyphs = read_glyph()
    for i, (_, position) in enumerate(matrix.iterrows()):
        letters_sorted = position.sort_values()
        bottom = 0
        for letter, height in letters_sorted.iteritems():
            patch = get_glyph_patch(glyphs[letter],colors[letter],i*charwidth,bottom,charwidth,height)
            ax.add_artist(patch)
            bottom += height

def plot_seqlogo(ax,pfm,charwidth=1.0,**kwargs):
    info_content = np.log2(20) - pfm.apply(lambda p: (-p*np.log2(p)).sum(),axis=1)
    matrix = pfm.mul(info_content,axis=0)

    seqlen = len(pfm)
    draw_logo(ax,matrix,charwidth,**kwargs)
    # xlim
    ax.set_xlim([0,seqlen*charwidth])
    # major ticks
    ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0, seqlen)))
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.tick_params(which='major', direction='out')
    # minor ticks
    ax.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(0, seqlen) + 0.5))
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(np.arange(1, seqlen+1)))
    ax.tick_params(which='minor', length=0)
    # ylim
    ax.set_ylim([0, 4.5])
    ax.yaxis.set_major_locator(ticker.FixedLocator([0., 1., 2.,3.,4.]))
    # show axis an left and bottom
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.splines['top'].set_visible(False)
    ax.splines['right'].set_visible(False)

def read_msa(msa_f):
    AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',\
                'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    with open(msa_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines if line]
        pro_line_num = [i for i, line in enumerate(
            lines) if '>' in line] + [len(lines)]
        seqs = [lines[n:pro_line_num[i + 1]]
                for i, n in enumerate(pro_line_num[:-1])]
        seq_number = len(seqs)
        seq_len = len(seqs[0][1])
        seqs = [(seq[0][1:], ''.join(seq[1:])) for seq in seqs]
    positions = [[seqs[i][1][j] for i in range(seq_number)] for j in range(seq_len)]
    pfm = [[pos.count(a)*1.0/seq_number for a in AA] for pos in positions]
    pfm = pd.DataFrame(pfm,columns=AA)
    return pfm


def main():
    pfm = read_msa(sys.argv[-1])
    fig = plt.figure()
    ax = fig.add_subplot(211)
    plot_seqlogo(ax,pfm)
    ax.set_aspect(1)
    ax.set_xlabel('position')
    ax.set_ylabel('bits')
    plt.savefig('test.png')
    plt.close('all')

main()





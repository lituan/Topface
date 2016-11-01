#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
plot patches
"""

import os
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


# R K H is shown blue, Y S T is show red
RES_HASH = {'K': 0, 'R': 0, 'H': 0, 'D': 8, 'E': 8, 'F': 8, 'W': 8, 'Y': 1, 'S': 1, 'T': 1,
            'N': 8, 'Q': 8, 'V': 8, 'L': 8, 'I': 8, 'M': 8, 'A': 8, 'C': 8, 'P': 8, 'G': 8, '*': 8}
RES_COLORS = {0: 'blue', 1: 'red', 2: 'green', 3: 'white',
              4: 'purple', 5: 'brown', 6: 'yellow', 7: 'cyan', 8: 'none'}
P_COLORS = {0: 'red', 1: 'blue', 2: 'purple'}

SHAPES = {
    'shape_3_1_1_a':	((1, 2), (2, 0), (2, 1)),
    'shape_3_1_1_b':	((1, 0), (1, 2), (2, 1)),
    'shape_3_1_2_a':	((1, 2), (2, 0), (2, 2)),
    'shape_3_1_2_b':	((1, 0), (1, 1), (2, 1)),
    'shape_3_1_3_a':	((1, 0), (2, 1), (2, 2)),
    'shape_3_1_3_b':	((1, 1), (1, 2), (2, 0)),
    'shape_3_2_1_a':	((1, 0), (2, 0), (2, 1)),
    'shape_3_2_1_b':	((1, 0), (1, 2), (2, 0)),
    'shape_3_2_2':        ((1, 0), (1, 1), (1, 2)),

    'shape_4_1_1_a':	((1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_4_1_1_b':	((1, 0), (1, 1), (1, 2), (2, 1)),
    'shape_4_1_2_a':	((1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_4_1_2_b':	((1, 0), (1, 2), (2, 1), (2, 2)),
    'shape_4_2_1_a':	((1, 0), (1, 1), (1, 2), (2, 0)),
    'shape_4_2_1_b':	((1, 0), (2, 0), (2, 1), (2, 2)),
    'shape_4_2_2_a':	((1, 0), (1, 1), (2, 0), (2, 1)),
    'shape_4_2_2_b':	((1, 0), (1, 2), (2, 0), (2, 2)),
    'shape_4_2_3':	((1, 0), (1, 2), (2, 0), (2, 1)),

    'shape_5_2_1_a':	((1, 0), (1, 1), (1, 2), (2, 0), (2, 1)),
    'shape_5_2_1_b':	((1, 0), (1, 2), (2, 0), (2, 1), (2, 2)),
    'shape_5_2_2_a':	((1, 0), (1, 2), (2, 1), (2, 2), (3, 0)),
    'shape_5_2_2_b':	((1, 0), (2, 1), (2, 2), (3, 0), (3, 1)),
}


def polar_to_rect(theta, r):
    return (r * np.cos(theta) + 0.5, r * np.sin(theta) + 0.25)


def plot_hotspot(ax, blades, title):

    area_in = 0.064  # area fo inner circle
    area_out = 0.064  # area of outter circle

    font_size = 26
    title_posi = 0.15
    circle_alpha = 0.6
    line_alpha = 0.6
    line_style = 'dotted'
    line_width = 0.4
    line_color = 'none'
    line_fill = False

    # ax.axis('off')
    fontdict = {'fontsize': font_size}
    ax.set_title(title, fontdict, position=(0.5, title_posi))

    for blade in blades:
        # plot circles
        circ_in = patches.Circle(blade[0][0], area_in, alpha=circle_alpha, color=blade[
                                 2][0], transform=ax.transAxes)
        circ_out_1 = patches.Circle(blade[0][1], area_in, alpha=circle_alpha, color=blade[
                                    2][1], transform=ax.transAxes)
        circ_out_2 = patches.Circle(blade[0][2], area_in, alpha=circle_alpha, color=blade[
                                    2][2], transform=ax.transAxes)
        ax.add_patch(circ_in)
        ax.add_patch(circ_out_1)
        ax.add_patch(circ_out_2)
        # add text
        ax.text(blade[0][0][0], blade[0][0][1], blade[1][0], transform=ax.transAxes,
                horizontalalignment='center', verticalalignment='center', **fontdict)
        ax.text(blade[0][1][0], blade[0][1][1], blade[1][1], transform=ax.transAxes,
                horizontalalignment='center', verticalalignment='center', **fontdict)
        ax.text(blade[0][2][0], blade[0][2][1], blade[1][2], transform=ax.transAxes,
                horizontalalignment='center', verticalalignment='center', **fontdict)
        # plot triangles connecting circles
        trip = patches.Polygon(blade[0], alpha=line_alpha, ls=line_style, lw=line_width,
                               fill=line_fill, facecolor=line_color, transform=ax.transAxes)
        ax.add_patch(trip)


def get_blades(shape_name, patch_seq):
    # blade format
    # [
    #    ((positions),...)
    #    (text,...)
    #    (colors,...)
    #        ]

    r_in = 0.2  # position of inner circles, less than r_out
    r_out = 0.4  # position of outter circles, more than r_in

    shape = SHAPES[shape_name]

    blades = []
    blade_num = set([s[0] for s in shape])
    for b in blade_num:
        blade_full = []
        for r in (0, 1, 2):
            if not (b, r) in shape:
                blade_full.append(((b, r), ''))
            else:
                blade_full.append(((b, r), patch_seq[shape.index((b, r))]))
        blades.append(blade_full)

    texts = [[b[1] for b in blade]
             for blade in blades]  # format (('R','R',''),...)

    # color by residue
    colors = [[RES_COLORS.get(RES_HASH.get(t, ''), 'none')
               for t in text] for text in texts]
    # color by position (empty is not shown)
    colors = [[P_COLORS[b[0][1]] if b[1] else 'none' for b in blade]
              for blade in blades]
    # color by position (empty is shown)
    # colors = [[P_COLORS[b[0][1]] for b in blade] for blade in blades]

    blade_positions = []
    if len(blade_num) == 1:
        blade_positions.append([polar_to_rect(
            np.pi * 0.5, r_in), polar_to_rect(np.pi * 0.375, r_out), polar_to_rect(np.pi * 0.625, r_out)])
    elif len(blade_num) == 2:
        blade_positions.append([polar_to_rect(
            np.pi * 0.30, r_in), polar_to_rect(np.pi * 0.2, r_out), polar_to_rect(np.pi * 0.4, r_out)])
        blade_positions.append([polar_to_rect(
            np.pi * 0.70, r_in), polar_to_rect(np.pi * 0.6, r_out), polar_to_rect(np.pi * 0.8, r_out)])
    elif len(blade_num) == 3:
        blade_positions.append([polar_to_rect(np.pi * 0.25, r_in), polar_to_rect(
            np.pi * 1 / 7, r_out), polar_to_rect(np.pi * 2 / 7, r_out)])
        blade_positions.append([polar_to_rect(np.pi * 0.50, r_in), polar_to_rect(
            np.pi * 3 / 7, r_out), polar_to_rect(np.pi * 4 / 7, r_out)])
        blade_positions.append([polar_to_rect(np.pi * 0.75, r_in), polar_to_rect(
            np.pi * 5 / 7, r_out), polar_to_rect(np.pi * 6 / 7, r_out)])

    return zip(blade_positions, texts, colors)


def plot_patch(shape_name, patch_seq, dirsuffix='', file_name=''):

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    # ax = fig.add_subplot(111)
    ax.set_xlim([0.1, 0.9])
    ax.set_ylim([0.25, 0.75])
    title = shape_name + ' ' + patch_seq

    blades = get_blades(shape_name, patch_seq)
    # plot_hotspot(ax, blades, title)
    plot_hotspot(ax, blades, '')

    if file_name != '':
        ofile = os.path.join(dirsuffix, file_name)
    else:
        ofile = os.path.join(dirsuffix, title)

    # fig.set_size_inches(6.0,3.0)
    fig.savefig(ofile, transparent=True, bbox_inches='tight',
                pad_inches=-0.65, dpi=300)
    plt.close('all')

# usage example
plot_patch('shape_4_2_2_a', 'TDWT')

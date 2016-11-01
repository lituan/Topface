#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
plot patch shapes
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def polar_to_rect(theta, r):
    return (r * np.cos(theta) + 0.5, r * np.sin(theta) + 0.1)


def plot_shape(ax, shape, blade_num, title):
    # patch format (((0,0),'R'),...)
    ax.axis('off')
    fontdict = {'fontsize': 12}
    # ax.set_title(title,fontdict,position=(0.5,0))
    ax.set_title(title, fontdict, position=(0.5, -0.1))
    # ax.set_title('good',fontdict,position=(0.5,0.0))

    color_in = {}
    color_out = {}
    for i in range(blade_num):
        color_in[i] = 'white'
    for i in range(blade_num * 2):
        color_out[i] = 'white'
    for i, p in enumerate(shape):
        b = p[0]
        r = p[1]
        if r == 0:
            color_in[b] = 'red'
        elif r == 1:
            color_out[b * 2 + r - 1] = 'blue'
        elif r == 2:
            color_out[b * 2 + r - 1] = 'purple'
    if blade_num == 1:
        num_in = blade_num
        theta_in = [np.pi * 0.5]
        num_out = blade_num * 2
        blade_bet = np.pi / 4
        theta_out = [blade_bet * 1.5, blade_bet * 2.5]
        r_in = 0.1
        area_in = 0.064
        r_out = 0.4
        area_out = 0.064
    elif blade_num == 2:
        num_in = blade_num
        theta_in = [np.pi * 0.25, np.pi * 0.75]
        num_out = blade_num * 2
        blade_bet = np.pi / 5
        theta_out = [blade_bet * 1, blade_bet *
                     2, blade_bet * 3, blade_bet * 4]
        r_in = 0.15
        area_in = 0.064
        r_out = 0.4
        area_out = 0.064
    elif blade_num == 3:
        num_in = blade_num
        blade_bet = np.pi / 4
        theta_in = [blade_bet * 1, blade_bet * 2, blade_bet * 3]
        num_out = blade_num * 2
        blade_bet = np.pi / 7
        theta_out = [blade_bet * (i + 1) for i in range(num_out)]
        r_in = 0.23
        area_in = 0.064
        r_out = 0.4
        area_out = 0.064

    center_in = []
    for i in range(num_in):
        center_in.append(polar_to_rect(theta_in[i], r_in))
        circ = patches.Circle(center_in[i], area_in, alpha=0.6, color=color_in[
                              i], transform=ax.transAxes)
        ax.add_patch(circ)
    center_out = []
    colors = ['blue', 'purple']
    for i in range(num_out):
        center_out.append(polar_to_rect(theta_out[i], r_out))
        circ = patches.Circle(center_out[i], area_out, alpha=0.6, color=color_out[
                              i], transform=ax.transAxes)
        ax.add_patch(circ)

    for i in range(num_in):
        a = center_in[i]
        b = center_out[i * 2]
        c = center_out[i * 2 + 1]
        vx = [(a[0], a[1]), (b[0], b[1]), (c[0], c[1])]
        trip = patches.Polygon(vx, alpha=0.6, ls='dotted', lw=0.2,
                               fill=False, facecolor='none', transform=ax.transAxes)
        ax.add_patch(trip)
        # ax.triplot([a[0],b[0],c[0]],[a[1],b[1],c[1]],transform=ax.transAxes)


def plot_shapes(shapes, filesuffix='', dirsuffix=''):
    # pro_hots format: {pro_name:['RRR','KKK','YYYY',...],...}
    if not type(shapes) is dict:
        return
    sha_names = shapes.keys()
    sha_names = sorted(sha_names, reverse=True)
    fig_num = len(shapes)
    c_num = 4
    r_num = 4
    if fig_num % (c_num * r_num) == 0:
        p_num = fig_num // (c_num * r_num)
    else:
        p_num = fig_num // (c_num * r_num) + 1
    for p in range(p_num):
        fig = plt.figure(figsize=(8, 8))
        for i in range(c_num * r_num):
            try:
                sha_name = sha_names.pop()
                ax = fig.add_subplot(r_num, c_num, i + 1, aspect='equal')
                shape = shapes[sha_name]
                blade_num = shape[-1][0] + 1
                title = str(sha_name)
                plot_shape(ax, shape, blade_num, title)
            except:
                if filesuffix == '':
                    ofile_name = str(p + 1)
                else:
                    ofile_name = filesuffix + '_' + str(p + 1)
                ofile = os.path.join(dirsuffix, ofile_name)
                fig.savefig(ofile, transparent=True,
                            bbox_inches='tight', pad_inches=0, dpi=1000)
                plt.close('all')
                return
        if filesuffix == '':
            ofile_name = str(p + 1)
        else:
            ofile_name = filesuffix + '_' + str(p + 1)
        ofile = os.path.join(dirsuffix, ofile_name)
        fig.savefig(ofile, transparent=True,
                    bbox_inches='tight', pad_inches=0, dpi=1000)
        plt.close('all')

shapes = {

    # 'shape_3_1_1_a':	((0, 2), (1, 0), (1, 1)),
    # 'shape_3_1_1_b':	((0, 0), (0, 2), (1, 1)),
    # 'shape_3_1_2_a':	((0, 2), (1, 0), (1, 2)),
    # 'shape_3_1_2_b':	((0, 0), (0, 1), (1, 1)),
    # 'shape_3_1_3_a':	((0, 0), (1, 1), (1, 2)),
    # 'shape_3_1_3_b':	((0, 1), (0, 2), (1, 0)),
    # 'shape_3_2_1_a':	((0, 0), (1, 0), (1, 1)),
    # 'shape_3_2_1_b':	((0, 0), (0, 2), (1, 0)),
    # 'shape_3_2_2':      ((0, 0), (0, 1), (0, 2)),

    # 'shape_3_1':	((0, 2), (1, 0), (1, 1)),
    # 'shape_3_2':	((0, 0), (0, 2), (1, 1)),
    # 'shape_3_3':	((0, 2), (1, 0), (1, 2)),
    # 'shape_3_4':	((0, 0), (0, 1), (1, 1)),
    # 'shape_3_5':	((0, 0), (1, 1), (1, 2)),
    # 'shape_3_6':	((0, 1), (0, 2), (1, 0)),
    # 'shape_3_7':	((0, 0), (1, 0), (1, 1)),
    # 'shape_3_8':	((0, 0), (0, 2), (1, 0)),
    # 'shape_3_9':        ((0, 0), (0, 1), (0, 2)),
    # 'shape_3_4':        ((0, 0), (1, 1), (2, 2)),

    #'shape_4_1_1_a':	((0, 2), (1, 0), (1, 1), (1, 2)),
    #'shape_4_1_1_b':	((0, 0), (0, 1), (0, 2), (1, 1)),
    #'shape_4_1_2_a':	((0, 1), (0, 2), (1, 0), (1, 1)),
    #'shape_4_1_2_b':	((0, 0), (0, 2), (1, 1), (1, 2)),
    #'shape_4_2_1_a':	((0, 0), (0, 1), (0, 2), (1, 0)),
    #'shape_4_2_1_b':	((0, 0), (1, 0), (1, 1), (1, 2)),
    'shape_4_2_2_a':	((0, 0), (0, 1), (1, 0), (1, 1)),
    #'shape_4_2_2_b':	((0, 0), (0, 2), (1, 0), (1, 2)),
    #'shape_4_2_3':        ((0, 0), (0, 2), (1, 0), (1, 1)),

    # 'shape_4_1':	((0, 2), (1, 0), (1, 1), (1, 2)),
    # 'shape_4_2':	((0, 0), (0, 1), (0, 2), (1, 1)),
    # 'shape_4_3':	((0, 1), (0, 2), (1, 0), (1, 1)),
    # 'shape_4_4':	((0, 0), (0, 2), (1, 1), (1, 2)),
    # 'shape_4_5':	((0, 0), (0, 1), (0, 2), (1, 0)),
    # 'shape_4_6':	((0, 0), (1, 0), (1, 1), (1, 2)),
    # 'shape_4_7':	((0, 0), (0, 1), (1, 0), (1, 1)),
    # 'shape_4_8':	((0, 0), (0, 2), (1, 0), (1, 2)),
    # 'shape_4_9':	        ((0, 0), (0, 2), (1, 0), (1, 1)),

    # 'shape_5_2_1_a':	((0, 0), (0, 1), (0, 2), (1, 0), (1, 1)),
    # 'shape_5_2_1_b':	((0, 0), (0, 2), (1, 0), (1, 1), (1, 2)),
    # 'shape_5_2_2_a':	((0, 0), (0, 2), (1, 1), (1, 2), (2, 0)),
    # 'shape_5_2_2_b':	((0, 0), (1, 1), (1, 2), (2, 0), (2, 1)),
}

plot_shapes(shapes)

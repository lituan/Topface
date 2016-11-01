import sys
import os
import itertools
import operator
import numpy as np
import lt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import defaultdict

def polar_to_rect(theta,r):
    return (r*np.cos(theta)+0.5,r*np.sin(theta)+0.5)

def plot_hotspot(ax,blade_num,patch,title,ft=6,title_posi=-0.06,line_width=1.0):
#patch format (((0,0),'R'),...)
    ax.axis('off')
    fontdict = {'fontsize':ft}
    ax.set_title(title,fontdict,position=(0.5,title_posi))

    basic = ['K','R','H']
    acid = ['D','E']
    aromatic = ['F','W','Y']
    polar =['S','T','N','Q']
    branch_phobic = ['V','L','I','M','A']
    special_branch = ['P','G']
    sulf = ['C']
    res_hash = {'K':0,'R':0,'H':0,'D':1,'E':1,'F':2,'W':2,'Y':2,'S':3,'T':3,\
                'N':3,'Q':3,'V':4,'L':4,'I':4,'M':4,'A':4,'C':5,'P':6,'G':7,'*':8}
    res_hash = {'K':0,'R':0,'H':0,'D':8,'E':8,'F':8,'W':8,'Y':0,'S':0,'T':0,\
                'N':8,'Q':8,'V':8,'L':8,'I':8,'M':8,'A':8,'C':8,'P':8,'G':8,'*':8}
    res_hash = {'K':0,'R':0,'H':0,'D':8,'E':8,'F':1,'W':1,'Y':1,'S':2,'T':2,\
                'N':8,'Q':8,'V':8,'L':8,'I':8,'M':8,'A':8,'C':8,'P':8,'G':8,'*':8}
    colors = {0:'blue',1:'red',2:'green',3:'white',4:'purple',5:'brown',6:'yellow',7:'cyan',8:'none'}

    color_in = {}
    color_out = {}
    for i in range(blade_num):
        color_in[i] = 'none'
    for i in range(blade_num*2):
        color_out[i] = 'none'

    text_in_num = []
    text_out_num = []
    text_in = []
    text_out = []
    for i,p in enumerate(patch):
        b = p[0][0]
        r = p[0][1]
        if r == 0:
            text_in_num.append(b)
            text_in.append(p[1])
            color_in[b] = colors.get(res_hash.get(p[1],8),'none')
        else:
            text_out.append(p[1])
            text_out_num.append(b*2+r-1)
            color_out[b*2+r-1] = colors.get(res_hash.get(p[1],8),'none')

    num_in = blade_num
    blade_bet = np.pi*2/num_in
    theta_in = [blade_bet*i for i in range(num_in)]
    r_in = 0.2
    area_in = 0.064
    center_in = []
    for i in range(num_in):
        center_in.append(polar_to_rect(theta_in[i],r_in))
        circ = patches.Circle(center_in[i],area_in,alpha=0.6,color=color_in[i],transform=ax.transAxes)
        ax.add_patch(circ)

    num_out = blade_num*2
    blade_bet = np.pi*2/num_out
    theta_out = [blade_bet*(i-0.50) for i in range(num_out)]
    r_out = 0.4
    area_out = 0.064
    center_out = []
    colors = ['blue','purple']
    for i in range(num_out):
        center_out.append(polar_to_rect(theta_out[i],r_out))
        circ   = patches.Circle(center_out[i],area_out,alpha=0.6,color=color_out[i],transform=ax.transAxes)
        ax.add_patch(circ)

    for i,n in enumerate(text_in_num):
        ax.text(center_in[n][0],center_in[n][1],text_in[i],transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',**fontdict)
    for i,n in enumerate(text_out_num):
        ax.text(center_out[n][0],center_out[n][1],text_out[i],transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',**fontdict)

    for i in range(num_in):
        a = center_in[i]
        b = center_out[i*2]
        c = center_out[i*2+1]
        vx=[(a[0],a[1]),(b[0],b[1]),(c[0],c[1])]
        trip=patches.Polygon(vx,alpha=0.9,ls='dotted',lw=line_width,fill=False,facecolor='none',transform=ax.transAxes)
        ax.add_patch(trip)
        #ax.triplot([a[0],b[0],c[0]],[a[1],b[1],c[1]],transform=ax.transAxes)

def plot_top_face(pro_hots,dirsuffix=''):
    #pro_hots format: {pro_name:['RRR','KKK','YYYY',...],...}
    for pro_name,pro_blade in pro_hots.iteritems():
        fig = plt.figure()
        ax = fig.add_subplot(111,aspect='equal')
        blade_num = len(pro_blade)
        title = str(pro_name) + ' ' + 'bn:' + str(blade_num)
        patch = []
        for i,vi in enumerate(pro_blade):
            patch.append(((i,0),vi[0]))
            patch.append(((i,1),vi[1]))
            patch.append(((i,2),vi[2]))
        plot_hotspot(ax,blade_num,patch,title,ft=12)
        ofile = os.path.join(dirsuffix,str(pro_name))
        fig.savefig(ofile,transparent=True,bbox_inches='tight',dpi=1000)
        plt.close('all')

def plot_top_faces(pro_hots,dirsuffix=''):
    #pro_hots format: {pro_name:['RRR','KKK','YYYY',...],...}
    pro_names = pro_hots.keys()
    fig_num = len(pro_hots)
    c_num = 3
    r_num = 3
    if fig_num%(c_num*r_num) == 0:
        p_num = fig_num//(c_num*r_num)
    else:
        p_num = fig_num//(c_num*r_num) + 1
    for p in range(p_num):
        fig = plt.figure()
        for i in range(c_num*r_num):
            try:
                pro_name = pro_names.pop()
                ax = fig.add_subplot(r_num,c_num,i+1,aspect='equal')
                pro_blade = pro_hots[pro_name]
                blade_num = len(pro_blade)
                title = str(pro_name) + ' ' + 'bn:' + str(blade_num)
                patch = []
                for i,vi in enumerate(pro_blade):
                    patch.append(((i,0),vi[0]))
                    patch.append(((i,1),vi[1]))
                    patch.append(((i,2),vi[2]))
                plot_hotspot(ax,blade_num,patch,title,title_posi=-0.10,line_width=0.6)
            except:
                ofile_name = str(p+1)
                ofile = os.path.join(dirsuffix,ofile_name)
                fig.savefig(ofile,transparent=True,bbox_inches='tight',dpi=1000)
                plt.close('all')
                return
        ofile_name = str(p+1)
        ofile = os.path.join(dirsuffix,ofile_name)
        fig.savefig(ofile,transparent=True,bbox_inches='tight',dpi=1000)
        plt.close('all')


def get_hotspot(wdsp_f):
    wdsp_lines = wdsp_f.readlines()
    wdsp_hotspot = {}
    blade = {}
    for line in wdsp_lines:
        words = line.split()
        if len(words) >= 2 and words[0] == '>':
            pro_name = words[1]
            pro_blades = []
            pro_seq = ''
        elif len(words) > 4:
            pro_blades.append(words[3:-1])
            pro_seq += ''.join(words[3:-1])
            blade[pro_name] = pro_blades
    for pro_name,pro_blades in blade.iteritems():
        hotspot = []
        for blade in pro_blades:
            R1 = blade[2][1]
            R1_2 = blade[1][-1]
            if len(blade[5]) <= 5 and blade[5][1] == 'D':
                D_1 = blade[5][0]
            elif len(blade[5]) == 3 or len(blade[5]) == 2:
                D_1 = blade[5][0]
            elif 3 <= len(blade[5]) <= 5 and blade[5][2] == 'D':
                D_1 = blade[5][1]
            elif 4 <= len(blade[5]) <= 5 and blade[5][3] == 'D':
                D_1 = blade[5][2]
            elif 5 <= len(blade[5]) <= 5 and blade[5][4] == 'D':
                Di_1 = blade[5][3]
            elif len(blade[5]) <= 5:
                D_1 = blade[5][1]
            elif len(blade[5]) <= 7:
                D_1 = blade[5][0]
            else:
                D_1 = '*'
            hotspot.append(R1+R1_2+D_1)
        wdsp_hotspot[pro_name] = hotspot
    return wdsp_hotspot

@lt.run_time
def main():
    wdsp_f = open(sys.argv[-1])
    hotspot_d = get_hotspot(wdsp_f)

    file_path,file_name = os.path.split(sys.argv[-1])
    script_short_name, script_extension = os.path.splitext(sys.argv[0])
    file_short_name, file_extension = os.path.splitext(file_name)
    result_path = os.path.join(file_path,file_short_name + '_' +script_short_name+'_'+'result')
    if not os.path.exists(result_path):
        os.makedirs(result_path)
    plot_top_faces(hotspot_d,result_path)

main()






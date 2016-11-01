import sys
import os
import itertools
import operator
import numpy as np
import lt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import defaultdict

basic = ['K','R','H']
acid = ['D','E']
aromatic = ['F','W','Y']
polar =['S','T','N','Q']
branch_phobic = ['V','L','I','M','A']
special_branch = ['P','G']
sulf = ['C']
RES_HASH = {'K':0,'R':0,'H':0,'D':1,'E':1,'F':2,'W':2,'Y':2,'S':3,'T':3,\
            'N':3,'Q':3,'V':4,'L':4,'I':4,'M':4,'A':4,'C':5,'P':6,'G':7,'*':8}
RES_HASH = {'K':0,'R':0,'H':0,'D':8,'E':8,'F':8,'W':8,'Y':0,'S':0,'T':0,\
            'N':8,'Q':8,'V':8,'L':8,'I':8,'M':8,'A':8,'C':8,'P':8,'G':8,'*':8}
# R K H is shown blue, Y S T is show green, D E is shown red
RES_HASH = {'K':0,'R':0,'H':0,'D':1,'E':1,'F':8,'W':8,'Y':2,'S':2,'T':2,\
            'N':8,'Q':8,'V':8,'L':8,'I':8,'M':8,'A':8,'C':8,'P':8,'G':8,'*':8}
# R K H is shown blue, Y S T is show red
RES_HASH = {'K':0,'R':0,'H':0,'D':8,'E':8,'F':8,'W':8,'Y':1,'S':1,'T':1,\
            'N':8,'Q':8,'V':8,'L':8,'I':8,'M':8,'A':8,'C':8,'P':8,'G':8,'*':8}
RES_COLORS = {0:'blue',1:'red',2:'green',3:'white',4:'purple',5:'brown',6:'yellow',7:'cyan',8:'none'}
def polar_to_rect(theta,r):
    return (r*np.cos(theta)+0.5,r*np.sin(theta)+0.5)

def plot_hotspot(ax,blade_num,patch,title,font_size=6,title_posi=-0.2,circle_alpha=0.6,line_alpha=0.6,line_style='dotted',line_width=0.4,line_color='none',line_fill=False,color_res=1):
#patch format (((0,0),'R'),...)
# plot hotspots on topface, things to change are in args
    ax.axis('off')
    fontdict = {'fontsize':font_size}
    ax.set_title(title,fontdict,position=(0.5,title_posi))

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
            if color_res:
                color_in[b] = RES_COLORS.get(RES_HASH.get(p[1],8),'none')
            else:
                color_in[b]= 'red'
        else:
            text_out.append(p[1])
            text_out_num.append(b*2+r-1)
            if color_res:
                color_out[b*2+r-1] = RES_COLORS.get(RES_HASH.get(p[1],8),'none')
            else:
                if r == 1:
                    color_out[b*2+r-1] = 'blue'
                elif r == 2:
                    color_out[b*2+r-1] = 'purple'
    # fbxw7 color
    # color_in = {0:'none',1:'none',2:'blue',3:'blue',4:'none',5:'none',6:'none',7:'none'}
    # color_out = {0:'none',1:'none',2:'none',3:'none',4:'none',5:'blue',6:'none',7:'red',8:'none',9:'none',10:'none',11:'none',12:'none',13:'none',14:'none',15:'none'}
    # fbw1a color
    # color_in = {0:'red',1:'none',2:'none',3:'none',4:'none',5:'none',6:'none',7:'none'}
    # color_out = {0:'none',1:'blue',2:'red',3:'red',4:'none',5:'none',6:'none',7:'none',8:'none',9:'none',10:'none',11:'none',12:'none',13:'none',14:'none',15:'none'}

    num_in = blade_num
    blade_bet = np.pi*2/num_in
    theta_in = [blade_bet*i for i in range(num_in)]
    r_in = 0.2 # position of inner circles, less than r_out
    area_in = 0.064 # area fo inner circle
    center_in = []
    for i in range(num_in):
        center_in.append(polar_to_rect(theta_in[i],r_in))
        circ = patches.Circle(center_in[i],area_in,alpha=circle_alpha,color=color_in[i],transform=ax.transAxes)
        ax.add_patch(circ)

    num_out = blade_num*2
    blade_bet = np.pi*2/num_out
    theta_out = [blade_bet*(i-0.50) for i in range(num_out)]
    r_out = 0.4 # position of outter circles, more than r_in
    area_out = 0.064 # area of outter circle
    center_out = []
    for i in range(num_out):
        center_out.append(polar_to_rect(theta_out[i],r_out))
        circ   = patches.Circle(center_out[i],area_out,alpha=circle_alpha,color=color_out[i],transform=ax.transAxes)
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
        trip=patches.Polygon(vx,alpha=line_alpha,ls=line_style,lw=line_width,fill=line_fill,facecolor=line_color,transform=ax.transAxes)
        ax.add_patch(trip)
        #ax.triplot([a[0],b[0],c[0]],[a[1],b[1],c[1]],transform=ax.transAxes)

def plot_topface(pro_hots,dirsuffix='',font_size=20,title_posi=0.99,circle_alpha=0.6,line_alpha=0.6,line_style='dotted',line_width=0.4,line_color='none',line_fill=False,color_res=1):
    #pro_hots format: {pro_name:['RRR','KKK','YYYY',...],...}
    # plot each topface on one figure
    for pro_name,pro_blade in pro_hots:
        fig = plt.figure()
        ax = fig.add_subplot(111,aspect='equal')
        blade_num = len(pro_blade)
        title = str(pro_name) + ' ' + 'bn:' + str(blade_num)
        title = str(pro_name)
        patch = []
        for i,vi in enumerate(pro_blade):
            patch.append(((i,0),vi[0]))
            patch.append(((i,1),vi[1]))
            patch.append(((i,2),vi[2]))
        plot_hotspot(ax,blade_num,patch,title,font_size=font_size,title_posi=title_posi,circle_alpha=circle_alpha,line_alpha=line_alpha,line_style=line_style,line_width=line_width,line_color=line_color,line_fill=line_fill,color_res=color_res)
        ofile = os.path.join(dirsuffix,str(pro_name))
        fig.savefig(ofile,transparent=True,bbox_inches='tight',dpi=300)
        plt.close('all')

def plot_topfaces(pro_hots,dirsuffix='',column_num=4, row_num=4, font_size=6,title_posi=0.99,circle_alpha=0.6,line_alpha=0.6,line_style='dotted',line_width=0.6,line_color='none',line_fill=False,color_res=1):
    #pro_hots format: {pro_name:['RRR','KKK','YYYY',...],...}
    # plot several(default 16) topfaces on one figure
    pro_hots = pro_hots[::-1]
    fig_num = len(pro_hots)
    c_num = column_num
    r_num = row_num
    if fig_num%(c_num*r_num) == 0:
        p_num = fig_num//(c_num*r_num)
    else:
        p_num = fig_num//(c_num*r_num) + 1
    for p in range(p_num):
        fig = plt.figure()
        for i in range(c_num*r_num):
            try:
                pro_name,pro_blade = pro_hots.pop()
                ax = fig.add_subplot(r_num,c_num,i+1,aspect='equal')
                blade_num = len(pro_blade)
                title = str(pro_name) + ' ' + 'bn:' + str(blade_num)
                title = str(pro_name)
                patch = []
                for i,vi in enumerate(pro_blade):
                    patch.append(((i,0),vi[0]))
                    patch.append(((i,1),vi[1]))
                    patch.append(((i,2),vi[2]))
                plot_hotspot(ax,blade_num,patch,title,font_size=font_size,title_posi=title_posi,circle_alpha=circle_alpha,line_alpha=line_alpha,line_style=line_style,line_width=line_width,line_color=line_color,line_fill=line_fill,color_res=color_res)
            except:
                ofile_name = str(p+1)
                ofile = os.path.join(dirsuffix,ofile_name)
                fig.savefig(ofile,transparent=True,bbox_inches='tight',dpi=300)
                plt.close('all')
                return
        ofile_name = str(p+1)
        ofile = os.path.join(dirsuffix,ofile_name)
        fig.savefig(ofile,transparent=True,bbox_inches='tight',dpi=300)
        plt.close('all')

# pro_hots = [('empty',['   ','   ','   ','   ','   ','   ','   '])]
# pro_hots = [('patch',['  K','RRY','   ','   ','   ','   ','   '])]
# pro_hots = [('FBXW7_HUMAN',['   ','   ','R R','R Y','   ','   ','   ','   '])]
# pro_hots = [('FBW1A_HUMAN',['Y R',' SS','   ','   ','   ','   ','   '])]
# plot_topface(pro_hots,color_res=0)






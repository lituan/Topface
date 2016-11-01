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
def plot_topface(blade_num,title):
#patch format (((0,0),'R'),...)
    fig = plt.figure()
    ax = fig.add_subplot(111,aspect='equal')
    ax.axis('off')
    fontdict = {'fontsize':4}

    colors = {0:'blue',1:'red',2:'green',3:'white',4:'purple',5:'brown',6:'yellow',7:'cyan',8:'none'}

    color_in = {}
    color_out = {}
    for i in range(blade_num):
        color_in[i] = 'red'
    for i in range(blade_num*2):
        color_out[i*2] = 'blue'
        color_out[i*2+1] = 'purple'

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
    
    for i in range(num_in):
        a = center_in[i]
        b = center_out[i*2]
        c = center_out[i*2+1]
        vx=[(a[0],a[1]),(b[0],b[1]),(c[0],c[1])]
        trip=patches.Polygon(vx,alpha=0.9,ls='dotted',lw=1.2,fill=False,facecolor='none',transform=ax.transAxes)
        ax.add_patch(trip)
        #ax.triplot([a[0],b[0],c[0]],[a[1],b[1],c[1]],transform=ax.transAxes)

    fig.savefig('topface',transparent=True)
plot_topface(7,'test')

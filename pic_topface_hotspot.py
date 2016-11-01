import sys
import os
import itertools
import operator
import lt
from collections import defaultdict
from collections import OrderedDict

from pic_topface_base import plot_topface,plot_topfaces

def get_hotspot(hot_f):
    lines = hot_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line for line in lines if len(line) > 0]
    lines = lines[1:]

    hotspot_d = OrderedDict()
    for line in lines:
        ls = line.split()
        pro_name,similarity,hotspot = ls[0],ls[1],ls[5:]
        similarity = str(int(float(similarity)*1))
        hotspot_d[similarity+'_'+pro_name] = hotspot

    return hotspot_d

@lt.run_time
def main():
    wdsp_f = open(sys.argv[-1])
    hotspot_d = get_hotspot(wdsp_f)
    hotspot_d = [(k,v) for k,v in hotspot_d.items()]
    file_path,file_name = os.path.split(sys.argv[-1])
    script_short_name, script_extension = os.path.splitext(sys.argv[0])
    file_short_name, file_extension = os.path.splitext(file_name)
    result_path = os.path.join(file_path,file_short_name + '_' +script_short_name+'_'+'result')
    if not os.path.exists(result_path):
        os.makedirs(result_path)
    # hotspot_d = {'FBXW7_HUMAN':['   ','   ','R R','R Y','   ','   ','   ','   ']}
    # hotspot_d = {'FBW1A_HUMAN':['Y R',' SS','   ','   ','   ','   ','   ']}
    # plot each topface on one figure
    plot_topface(hotspot_d,dirsuffix=result_path,font_size=26,title_posi=0.99,circle_alpha=0.6,line_alpha=0.6,line_style='dotted',line_width=0.4,line_color='none',line_fill=False,color_res=1)
    # plot several topfaces on one figure
    plot_topfaces(hotspot_d,dirsuffix=result_path,column_num=4, row_num=4, font_size=6,title_posi=0.99,circle_alpha=0.6,line_alpha=0.6,line_style='dotted',line_width=0.6,line_color='none',line_fill=False,color_res=1)

main()





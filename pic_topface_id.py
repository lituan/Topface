import sys
import os
import itertools
import operator
import lt
from collections import defaultdict
from collections import OrderedDict

from pic_topface_base import plot_topface,plot_topfaces

def get_id(id_f):
    lines = id_f.readlines()
    id_l = []
    for line in lines:
        word = line.split()
        id_l.append(word[0])
    return id_l

def get_hotspot(wdsp_f):
    wdsp_lines = wdsp_f.readlines()
    wdsp_hotspot = OrderedDict()
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
    id_f = open(sys.argv[-2])
    wdsp_f = open(sys.argv[-1])
    id_l = get_id(id_f)
    hotspot_d = get_hotspot(wdsp_f)
    hotspot_d = [(k,v) for k,v in hotspot_d]
    found = {}
    not_found = []
    for ik in id_l:
        if not ik in hotspot_d.keys():
            not_found.append(ik)
        else:
            found[ik] = hotspot_d[ik]
    file_path,file_name = os.path.split(sys.argv[-1])
    script_short_name, script_extension = os.path.splitext(sys.argv[0])
    file_short_name, file_extension = os.path.splitext(file_name)
    result_path = os.path.join(file_path,file_short_name + '_' +script_short_name+'_'+'result')
    if not os.path.exists(result_path):
        os.makedirs(result_path)
    # plot each topface on one figure
    plot_topface(found,dirsuffix=result_path,font_size=26,title_posi=0.99,circle_alpha=0.6,line_alpha=0.6,line_style='dotted',line_width=0.4,line_color='none',line_fill=False,color_res=1)
    # plot several topfaces on one figure
    plot_topfaces(found,dirsuffix=result_path,column_num=4, row_num=4, font_size=6,title_posi=0.99,circle_alpha=0.6,line_alpha=0.6,line_style='dotted',line_width=0.6,line_color='none',line_fill=False,color_res=1)
    lt.print_list(not_found)

if __name__ == "__main__":
    main()





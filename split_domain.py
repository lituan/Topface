##############################################
## here give the WDSP outfile name
ifile = file('wd648_chordata.wdsp', 'r') 
##############################################

import os
WD40 = {}
Names = {}

#ofile = file('WDSP_domains_for_WDSP3D.out', 'w')
ofile = file('wd648_chordata_split.wdsp', 'w')

for line in ifile:
    words = line.split()
    if line[0] == '>' :
        ID = words[1]
        WD40[ID] = []
        Names[ID] = words[2:]
    elif len(words) > 0:
        WD40[ID].append(words)
print len(WD40)

def Output_blade(B):
    corr = float(B[-1])
    hotspot = []
    tetrad = []
    start = [int(B[1])]
    if corr == 44:
        mark = 'pentad'
    if corr == 28:
        mark = 'tetrad'
    if corr == 18:
        mark = 'triad'
    if corr == 19:
        mark = 'triad'
    if corr == 0:
        mark = 'NA'
        
    for i in B[:3]:
        print>>ofile, i,

    
    blade_seq = B[3:-1]
    for i in blade_seq:
        start.append(start[-1]+len(i))
    #print>>ofile, start

    
    #pentad
    if 5 <= len(blade_seq[1]) <= 8 and ('H' in blade_seq[1][-5:-3] or (len(blade_seq[1]) > 5 and 'H' in blade_seq[1][1:3])) and\
        blade_seq[4][3] in 'TS' and \
        blade_seq[5][1] == 'D' and len(blade_seq[5][1])<= 6 and\
        blade_seq[6][4] in 'W' and \
        blade_seq[6][0] in 'ST':
        if blade_seq[1][-4] == 'H':
            tetrad.append('H'+str(start[1]+len(blade_seq[1])-4))
        elif blade_seq[1][-5] == 'H':
            tetrad.append('H'+str(start[1]+len(blade_seq[1])-5))
        elif blade_seq[1][1] == 'H':
            tetrad.append('H'+str(start[1]+1))
        elif blade_seq[1][2] == 'H':
            tetrad.append('H'+str(start[1]+2))
        tetrad.append(blade_seq[4][3]+str(start[4]+3))
        tetrad.append('D'+str(start[5]+1))
        tetrad.append(blade_seq[6][0]+str(start[6]))
        tetrad.append('W'+str(start[6]+4))             
    #tetrad
    elif 5 <= len(blade_seq[1]) <= 8 and ('H' in blade_seq[1][-5:-3] or (len(blade_seq[1]) > 5 and 'H' in blade_seq[1][1:3])) and\
        blade_seq[4][3] in 'TS' and\
        (blade_seq[5][1] == 'D' and len(blade_seq[5][1])<= 6) and\
        blade_seq[6][4] in 'W':
        if blade_seq[1][-4] == 'H':
            tetrad.append('H'+str(start[1]+len(blade_seq[1])-4))
        elif blade_seq[1][-5] == 'H':
            tetrad.append('H'+str(start[1]+len(blade_seq[1])-5))
        elif blade_seq[1][1] == 'H':
            tetrad.append('H'+str(start[1]+1))
        elif blade_seq[1][2] == 'H':
            tetrad.append('H'+str(start[1]+2))
        tetrad.append(blade_seq[4][3]+str(start[4]+3))
        tetrad.append('D'+str(start[5]+1))
        tetrad.append('W'+str(start[6]+4))     

    #triad
    elif 5 <= len(blade_seq[1]) <= 8 and ('H' in blade_seq[1][-5:-3] or (len(blade_seq[1]) > 5 and 'H' in blade_seq[1][1:3])) and\
        blade_seq[4][3] in 'TS' and blade_seq[6][4] in 'W' :
        if blade_seq[1][-4] == 'H':
            tetrad.append('H'+str(start[1]+len(blade_seq[1])-4))
        elif blade_seq[1][-5] == 'H':
            tetrad.append('H'+str(start[1]+len(blade_seq[1])-5))
        elif blade_seq[1][1] == 'H':
            tetrad.append('H'+str(start[1]+1))
        elif blade_seq[1][2] == 'H':
            tetrad.append('H'+str(start[1]+2))
        tetrad.append(blade_seq[4][3]+str(start[4]+3))
        tetrad.append('W'+str(start[6]+4))     
    elif 5 <= len(blade_seq[1]) <= 8 and (blade_seq[5][1] == 'D' and len(blade_seq[5])<= 6) and blade_seq[4][3] in 'TS' and\
           ('H' in blade_seq[1][-5:-3] or (len(blade_seq[1]) > 5 and 'H' in blade_seq[1][1:3]) ) :
        if blade_seq[1][-4] == 'H':
            tetrad.append('H'+str(start[1]+len(blade_seq[1])-4))
        elif blade_seq[1][-5] == 'H':
            tetrad.append('H'+str(start[1]+len(blade_seq[1])-5))
        elif blade_seq[1][1] == 'H':
            tetrad.append('H'+str(start[1]+1))
        elif blade_seq[1][2] == 'H':
            tetrad.append('H'+str(start[1]+2))
        tetrad.append(blade_seq[4][3]+str(start[4]+3))
        tetrad.append('D'+str(start[5]+1))
 
    #hotspot
    if blade_seq[1][-1] in 'WRYFKILDEQNHM':
        hotspot.append(blade_seq[1][-1] + str(start[2]-1))
    if blade_seq[2][1] in 'WRYFKILDEQNHM':
        hotspot.append(blade_seq[2][1] + str(start[2]+1))
    if blade_seq[5][0] in 'WRYFKILDEQNHM' and blade_seq[5][1] in 'DN'  and len(blade_seq[5])<=6 :
        hotspot.append(blade_seq[5][0] + str(start[5]))
    
        
    for i in blade_seq:
        print>>ofile, i,
    out_1 = '['
    for i in tetrad[:-1]:
        out_1 += i+','
    if len(tetrad)>0:
        out_1 += tetrad[-1]+']'
    else:
        out_1 += ']'
    print>>ofile, out_1,
    
    out_2 = '['
    for i in hotspot[:-1]:
        out_2 += i+','
    if len(hotspot)>0:
        out_2 += hotspot[-1]+']'
    else:
        out_2 += ']'
    print>>ofile, out_2
    


for p in WD40:
    L = len(WD40[p])
    #print>>ofile, L
    if 4 <= L <= 8:
        print>>ofile, '>', p, Names[p][0], L
        #print>>ofile, 1
        for i in WD40[p]:
            Output_blade(i)

    elif L % 7 == 0:
        for i in range(L):
            if i % 7 == 0:
                print>>ofile, '>', p+'_'+str(i//7), Names[p][0], 7
                #print>>ofile, str(i//7+1)
            Output_blade(WD40[p][i])
    elif L % 6 == 0:
        for i in range(L):
            if i % 6 == 0:
                print>>ofile, '>', p+'_'+str(i//6), Names[p][0], 6
                #print>>ofile, str(i//6+1)
            Output_blade(WD40[p][i])
    elif L % 8 == 0:
        for i in range(L):
            if i % 8 == 0:
                print>>ofile, '>', p+'_'+str(i//8), Names[p][0], 8
                #print>>ofile, str(i//8+1)
            Output_blade(WD40[p][i])
            
    elif L == 13:
        for i in range(L):
            if i % 7 == 0:
                print>>ofile, '>', p+'_'+str(i//7), Names[p][0],
                if i == 7:
                    print>>ofile, 6
                else:
                    print>>ofile, 7
                #print>>ofile, str(i//7+1)
            Output_blade(WD40[p][i])
            
    elif L == 15:
        for i in range(L):
            if i % 8 == 0:
                print>>ofile, '>', p+'_'+str(i//8), Names[p][0]
                #print>>ofile, str(i//8+1)
            Output_blade(WD40[p][i])
    elif L == 9:
        if int(WD40[p][0][0]) > int(WD40[p][-1][0]) and int(WD40[p][-1][0]) < 48:
            del(WD40[p][-1])
        elif int(WD40[p][0][0]) < int(WD40[p][-1][0]) and int(WD40[p][0][0]) < 48:
            del(WD40[p][0])
        if len(WD40[p]) == 8:
            print>>ofile, '>', p, Names[p][0], 8
            #print>>ofile, 1
            for i in range(8):
                Output_blade(WD40[p][i])
            L = len(WD40[p])
            
    elif 9 <= L <= 11:
        print>>ofile, '>', p+'_'+str(0), Names[p][0], 7, 'overlapped'
        #print>>ofile, 1
        for i in range(7):
            Output_blade(WD40[p][i])
        print>>ofile, '>', p+'_'+str(1), Names[p][0], 7, 'overlapped'
        #print>>ofile, 2
        for i in range(-7,0):
            Output_blade(WD40[p][i])
    elif L >= 17:
        lens = L//7
        for i in range(lens*7):
            if i % 7 == 0:
                print>>ofile, '>', p+'_'+str(i//7), Names[p][0], 7
                #print>>ofile, str(i//7+1)
            Output_blade(WD40[p][i])
        print>>ofile, '>', p, Names[p][0],7, 'overlapped'
        #print>>ofile, str(lens+1)
        for i in range(-7,0):
            Output_blade(WD40[p][i])   
    else:
        raise      
ofile.close()


print 'end'


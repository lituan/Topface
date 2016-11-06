#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
creat logo for aligned hotspots
input formats be as follows:

input file has following formats:
A4SBD7              EVC KER KDY WTG LAG NDD

you can change read_align or read_hot for your formats
"""
import sys
import os
import operator
import Tkinter
from Tkinter import Tk, Canvas, Frame, BOTH, NW, SW
from PIL import Image, ImageTk
from math import log


class HOTLOGO:

    def __init__(self, align_f):
        self.seq = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',\
                    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        self.align_f = align_f
        self.title = os.path.splitext(align_f.name)[0] + ' Sequence Logo'
        self.position = {}  # format {0:['ARNAA...'],...}
        self.matrix = {}  # format {0:{'A':2,...},...}
        self.posibits = {}  # foramt {0:1.4,...}
        self.posiinfo = {}
        self.posiheight = {}
        self.sort_posiheight = {}
        self.logo_len = 0

    def set_size(self, x=0, y=0):
        self.logo_w = x
        self.logo_h = y

    def create_logo(self):
        self.read_hot(self.align_f)
        self.get_matrix()
        self.get_bits()
        self.get_info()
        self.get_posiheight()
        self.get_sort_posiheight()

        # self.print_matrix()

        self.set_parameters()
        self.makelogo()

    def read_align(self, seq_f):
        lines = seq_f.readlines()
        lines = filter(lambda x: len(x) > 0, lines)
        lines = [line.rstrip('\n\r') for line in lines]
        seq_lines = {}
        for line in lines:
            if '>' in line:
                name = line
                seq_lines[name] = []
            else:
                seq_lines[name].append(line)
        seq_line = {}
        for k, w in seq_lines.iteritems():
            seq_line[k] = ''.join(w)
        seqs = seq_line.values()

        self.logo_len = len(seqs[0])
        for i in range(self.logo_len):
            self.position[i] = ''.join(s[i] for s in seqs)

    def read_hot(self, hot_f):
        lines = hot_f.readlines()
        lines = filter(lambda x: len(x) > 0, lines)
        self.hot = [''.join(h) for h in [l.split()[1:] for l in lines]]
        self.logo_len = len(self.hot[0])
        self.hot_num = len(self.hot)
        for i in range(self.logo_len):
            self.position[i] = ''.join([h[i] for h in self.hot])

    def get_matrix(self):
        for i in range(self.logo_len):
            self.matrix[i] = {}
            for aa in self.seq:
                self.matrix[i][aa] = 0
        for i in range(self.logo_len):
            self.matrix[i] = {}
            for aa in self.seq:
                self.matrix[i][aa] = self.position[i].count(aa)/(len(self.position[i])*1.0)

    def get_bits(self):
        for posi,aa_freq in self.matrix.iteritems():
            freq = filter(lambda x: x > 0,aa_freq.values())
            freq_h = [x*log(x,2) for x in freq]
            hbits = 0 - sum(freq_h)
            self.posibits[posi] = hbits

    def get_info(self):
        for posi,bits in self.posibits.iteritems():
            self.posiinfo[posi] = log(20,2) - bits

    def get_posiheight(self):
        for posi,aa_freq in self.matrix.iteritems():
            self.posiheight[posi] = {}
            for aa,freq in aa_freq.iteritems():
                self.posiheight[posi][aa] = self.posiinfo[posi] * freq

    def get_sort_posiheight(self):
        for posi,aa_height in self.posiheight.iteritems():
            aa_h = []
            for aa,height in aa_height.items():
                if height > 0.001:
                    aa_h.append((aa,height))
            aa_h = sorted(aa_h,key=operator.itemgetter(1))
            self.sort_posiheight[posi]=aa_h

    def print_matrix(self):
        print '    '.join(self.seq)
        for posi,aa_freq in self.matrix.iteritems():
            freq = [aa_freq[aa] for aa in self.seq]
            freq = [str(i)[:4] for i in freq]
            freq = [i+' '*(5-len(i)) for i in freq]
            print ''.join(freq)

    def set_parameters(self):
        logo_w = self.logo_len*30
        logo_h = 200
        if self.logo_w == 0:
            self.logo_w = logo_w
        if self.logo_h == 0:
            self.logo_h = logo_h

        self.background = 'white'
        self.title_font = ('arial',24)
        self.label_font = ('arial',14)

        self.top = 20 + 50 # top pad, more pad for title, at least
        self.left = 20 + 50 # more pad in left for left label, at least
        self.bottom = 20 + 50 # more pad in bottom for bottom label, at least
        self.right = 20 # right pad will be greater than this value in the end

        self.column_pad = 6 # space between columns,should be even

        self.line_width = 2
        self.left_tick_length = 5
        self.left_label_offset = 2
        self.bottom_label_offset = 2

        self.left_axis_offset = 5
        self.bottom_axis_offset = 2

        self.win_w = self.logo_w + self.left + self.right
        self.win_h = self.logo_h + self.top + self.bottom
        self.aa_w = self.logo_w/self.logo_len - self.column_pad
        if self.aa_w % 2 == 0: # even aa_w is better
            pass
        else:
            self.aa_w -= 1

        self.bottom_origin = self.top + self.logo_h
        self.left_origin = self.left

        self.bottom_axis_origin = self.bottom_origin + self.bottom_axis_offset
        self.left_axis_origin = self.left_origin - self.left_axis_offset

    def makelogo(self):
        root = Tkinter.Tk()
        canvas = Tkinter.Canvas(root,width=self.win_w,height=self.win_h,bg=self.background)

        AA_list = ['G','A','P','V','I','L','F','C','M','W','Y','Q','T','S','N','H','D','E','K','R']
        AA_figures = {}
        for aa in AA_list:
            AA_figures[aa] = Image.open('./color_aa/'+aa+'.png')

        self.posi_figures = {}
        for posi in range(self.logo_len):
            x = self.left_origin + self.aa_w*posi + self.column_pad*posi
            self.posi_figures[posi] = []
            aa_height = self.sort_posiheight[posi]
            y = self.bottom_origin
            for i,(aa,height) in enumerate(aa_height):
                aa_h = int(round(self.logo_h*(height/5.0)))
                if aa_h > 0:
                    aa_fig = AA_figures[aa].resize((self.aa_w,aa_h),Image.ANTIALIAS)
                    aa_fig = ImageTk.PhotoImage(aa_fig)
                    self.posi_figures[posi].append(aa_fig) #this is necessary or the pcitrue only show part
                    #canvas.create_image(x,y,anchor=SW,image=self.posi_figures[posi][-1])
                    canvas.create_image(x,y,anchor='sw',image=aa_fig)
                    y -= aa_h

        canvas.create_line((self.left_axis_origin, self.bottom_axis_origin, self.left_axis_origin, self.bottom_origin-self.logo_h), width=self.line_width) #create left axis
        #canvas.create_line((self.left_axis_origin, self.bottom_axis_origin, self.left_origin+self.logo_w, self.bottom_axis_origin), width=self.line_width) #create bottom axis

        canvas.create_text((self.left_axis_origin-15,self.bottom_origin-self.logo_h),text='bit',anchor='ne',font=self.label_font) #create left legend
        for i in range(4):
            i += 1
            h = self.bottom_origin - int((self.logo_h/5.0)*i)
            canvas.create_line((self.left_axis_origin, h, self.left_axis_origin-self.left_tick_length, h),width=self.line_width) #create left tickle
            canvas.create_text((self.left_axis_origin-self.left_tick_length-self.left_label_offset,h),text=str(i),anchor='e',font=self.label_font) #create left label

        for i in range(self.logo_len): # create bottom label
            x = self.left_origin + (self.aa_w + self.column_pad)*i + self.aa_w/2
            y = self.bottom_axis_origin + self.bottom_label_offset
            canvas.create_text((x,y),text=str(i+1),anchor='n',font=self.label_font) #create left legend

        canvas.create_text((self.left+self.logo_w/2,self.top),text=self.title,anchor='center',font=self.title_font)

        canvas.pack()
        canvas.update()
        canvas.postscript(file='test.ps',colormode='color')
        root.mainloop()

a = HOTLOGO(open(sys.argv[-1],'r'))
a.set_size()
a.create_logo()

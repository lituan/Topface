"""
calculate amino acid bit for each position for a group of hotspots
not finished yet

input file has following formats:

A4SBD7              EVC KER KDY WTG LAG NDD
A7AM88              WRS RTF KER KDY WTA QRT

"""
import sys
import os
from math import log

class HOTLOGO:

    def __init__(self,hot_f):
        self.seq = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        self.matrix = {}
        self.position = {}
        self.posibits = {}
        self.posiinfo = {}
        self.posiheight = {}

        self.read_hot(hot_f)
        self.get_matrix()
        self.print_matrix()
        self.get_bits()
        self.get_info()
        self.get_posiheight()

    def read_hot(self,hot_f):
        lines = hot_f.readlines()
        lines = filter(lambda x: len(x) > 0,lines)
        self.hot = [''.join(h) for h in [l.split()[2:] for l in lines]]
        self.logo_len = len(self.hot[0])
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
                self.matrix[i][aa] = self.position[i].count(aa)/(len(self.hot)*1.0)

    def get_posiheight(self):
        # calculate bits per position
        for posi,aa_freq in self.matrix.iteritems():
            freq = filter(lambda x: x > 0,aa_freq.values())
            freq_h = [x*log(x,2) for x in freq]
            hbits = 0 - sum(freq_h)
            self.posibits[posi] = hbits

        # calculate info per position
        for posi,bits in self.posibits.iteritems():
            self.posiinfo[posi] = log(20,2) - bits

        # calculate height per residue per position
        for posi,aa_freq in self.matrix.iteritems():
            self.posiheight[posi] = {}
            for aa,freq in aa_freq.iteritems():
                self.posiheight[posi][aa] = self.posiinfo[posi] * freq


    def print_matrix(self):
        print '    '.join(self.seq)
        for posi,aa_freq in self.matrix.iteritems():
            freq = [aa_freq[aa] for aa in self.seq]
            freq = [str(i)[:4] for i in freq]
            freq = [i+' '*(5-len(i)) for i in freq]
            print ''.join(freq)


a = HOTLOGO(open(sys.argv[- 1]))































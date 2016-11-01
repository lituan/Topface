import sys
import os
import numpy as np

from Bio.SubsMat import MatrixInfo

"""
usage: python hotspot_compare.py hotspot1 hotspot2

compare two hotspot files
hotspot file format

2ABA_HUMAN          SIK NKN NHD TVS SSY FDY LKT
A8JEA1_CHLRE        APG KVE IAR GPM QAS TPS KYN
"""

SEQ_MATRIX = MatrixInfo.blosum62
UP = 1
LEFT = 2
DIAG = 0
INDEL = -7


def gap_penalty(i, j):
    return -3 * abs(i - j)

def score_hot(h1, h2):
    score = []
    for h11,h22 in zip(h1,h2):
        if (h11,h22) in SEQ_MATRIX.keys():
            score.append(SEQ_MATRIX[(h11,h22)])
        elif (h22,h11) in SEQ_MATRIX.keys():
            score.append(SEQ_MATRIX[(h22,h11)])
        else:
            score.append(-5)

    return sum(score)

def needleman_wunsch_matrix(hot1, hot2, score_fun):
    """
    calculate the matrix
    """

    m, n = len(hot1), len(hot2)
    score_matrix = np.zeros((m + 1, n + 1))
    pointer_matrix = np.zeros((m + 1, n + 1), dtype=int)

    for i in range(1, m + 1):
        score_matrix[i, 0] = INDEL * i + (-3) * i * (i + 1) * 0.5
    for j in range(1, n + 1):
        score_matrix[0, j] = INDEL * j + (-3) * j * (j + 1) * 0.5

    pointer_matrix[0, 1:] = LEFT
    pointer_matrix[1:, 0] = UP

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            up = score_matrix[i - 1, j] + INDEL + gap_penalty(i, j)
            left = score_matrix[i, j - 1] + INDEL + gap_penalty(i, j)
            diag = score_matrix[i - 1, j - 1] + \
                score_fun(hot1[i - 1], hot2[j - 1])

            max_score = max(up, left, diag)
            if max_score == diag:
                pointer_matrix[i, j] = DIAG
            elif max_score == up:
                pointer_matrix[i, j] = UP
            elif max_score == left:
                pointer_matrix[i, j] = LEFT

            score_matrix[i, j] = max_score
    return score_matrix, pointer_matrix

def needleman_wunsch_trace(hot1, hot2, score_matrix, pointer_matrix):
    align1 = []
    align2 = []
    m, n = len(hot1), len(hot2)

    curr = pointer_matrix[m, n]
    i, j = m, n
    while (i > 0 or j > 0):
        if curr == DIAG:
            align1.append(hot1[i - 1])
            align2.append(hot2[j - 1])
            i -= 1
            j -= 1
        elif curr == LEFT:
            align1.append('   ')
            align2.append(hot2[j - 1])
            j -= 1
        elif curr == UP:
            align1.append(hot1[i - 1])
            align2.append('   ')
            i -= 1
        curr = pointer_matrix[i, j]

    return align1[::-1], align2[::-1]


def needleman_wunsch(hot1, hot2):
    score_matrix, pointer_matrix = needleman_wunsch_matrix(
        hot1, hot2, score_hot)
    hot1, hot2 = needleman_wunsch_trace(
        hot1, hot2, score_matrix, pointer_matrix)

    hot11 = ''.join(hot1)
    hot22 = ''.join(hot2)
    identity = sum([1 for h1, h2 in zip(hot11, hot22) if h1 == h2])

    length = min(len(hot1), len(hot2)) * 3

    return hot1, hot2, identity, length


def read_hotspot(o_f):
    lines = o_f.readlines()
    lines = [line.strip('\n\r') for line in lines]
    lines = [line for line in lines if len(line.split()) > 0]
    pro_hotspot = [(line.split()[0], line.split()[1:]) for line in lines]
    return dict(pro_hotspot)


def main():
    with open(sys.argv[1]) as o_f:
        pro_hotspot1 = read_hotspot(o_f)
    with open(sys.argv[2]) as o_f:
        pro_hotspot2 = read_hotspot(o_f)

    result = []
    for pro1, hotspot1 in pro_hotspot1.iteritems():
        if pro1 in pro_hotspot2.keys():
            hotspot2 = pro_hotspot2[pro1]
            result.append(needleman_wunsch(hotspot1, hotspot2) + (pro1,))

    total_identity = sum(r[2] for r in result)
    total_length = sum(r[3] for r in result)
    accuracy = round(total_identity * 1.0*100 / total_length)

    with open('hotspot_compare_result.txt','w') as w_f:
        print >>w_f,'accuracy', '\t', '{0:<}'.format(int(accuracy))
        for r in result:
            print >> w_f, '{0:<20}{1:<10}{2:<5}{3:<5}{4}'.format(
                r[4], int(round(r[2]*1.0*100/r[3])),r[2], r[3], ' '.join(r[0]))
            print >> w_f, '{0:40}{1}'.format('', ' '.join(r[1]))

if __name__ == "__main__":
    main()

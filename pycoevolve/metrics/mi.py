#/usr/bin/env python
from itertools import combinations_with_replacement as car
from collections import Counter
import numpy as np

aminobet = "ACDEFGHIKLMNPQRSTVWY"

def _transpose_alignment(alignment):
    """
    Returns list of columns from alignment,
    which is a dict of sequences.
    """
    just_seqs = alignment.values()
    transpose_align = []
    for i in xrange(len(just_seqs[0])):
        transpose_align.append("".join([seq[i]
            for seq in just_seqs]))
    return transpose_align

def nmi(alignment,alphabet=aminobet,epsilon=1e-8):
    alpha_size = len(alphabet)

    trans_align = _transpose_alignment(alignment)
    position_size = len(trans_align)

    combo_alpha = {}
    for ind,i in enumerate(alphabet):
        for jnd,j in enumerate(alphabet):
            combo_alpha[i+j] = (ind,jnd)

    mi_matrix = np.zeros((position_size,position_size),dtype=float)

    for two_pos in car(range(position_size),2):
        pos1 = two_pos[0]
        pos2 = two_pos[1]
        pmf = np.zeros((alpha_size,alpha_size),dtype=int)
        two_letter_counts = Counter([i[0] + i[1] for i in
                zip(trans_align[pos1],trans_align[pos2])])
        for two_letter in combo_alpha:
            i,j = combo_alpha[two_letter]
            pmf[i,j] = two_letter_counts[two_letter]
        npmf = np.array(pmf,dtype=float) + epsilon
        npmf /= npmf.sum() 
        npmf += epsilon
        row_sum = npmf.sum(axis=1)
        col_sum = npmf.sum(axis=0)
        sum_ent = -(row_sum*np.log(row_sum)).sum() - (col_sum*np.log(col_sum)).sum()
        joint_ent = -(npmf*np.log(npmf)).sum()
        mi_matrix[pos1,pos2] = (sum_ent - joint_ent)/sum_ent
    return mi_matrix

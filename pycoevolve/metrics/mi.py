#/usr/bin/env python
from itertools import combinations as comb
from collections import Counter
import numpy as np

aminobet = "ACDEFGHIKLMNPQRSTVWY"

def transpose_alignment(alignment):
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

def nmi(alignment,alphabet=aminobet,epsilon=1e-9,gap_cut=.1,cnsrv_cut=.95):
    """
    Returns the joint entropy normalized mutual information matrix.
    Alphabet refers to the aminoacids to consider.
    Epsilon is there so the 0*log(0) is approximately 0.
    gap_cut refers to the frequency of gaps tolerated at scorable positions.
        gap_cut = 0 means that only positions with no gaps will be scored.
    consrv_cut refers the maximum frequency of a single aminoacid at scorable
        positions. cnsrv_cut = .95 means that if a column has any single aminoacid
        more than 95% of the time, the position will be ignored
        This is because calculation of NMI becomes unreliable at low entropy.
    """
    alpha_size = len(alphabet)

    trans_align = transpose_alignment(alignment)
    position_size = len(trans_align)
    seq_size = float(len(trans_align[0]))

    combo_alpha = {}
    for ind,i in enumerate(alphabet):
        for jnd,j in enumerate(alphabet):
            combo_alpha[i+j] = (ind,jnd)

    mi_matrix = np.zeros((position_size,position_size),dtype=float)

    for two_pos in comb(range(position_size),2):
        pos1 = two_pos[0]
        pos2 = two_pos[1]
        col1 = trans_align[pos1]
        col2 = trans_align[pos2]
        if ((col1.count("-")/seq_size) > gap_cut) \
                or ((col2.count("-")/seq_size) > gap_cut):
            pass
        else:
            pmf = np.zeros((alpha_size,alpha_size),dtype=int)
            two_letter_counts = Counter([i[0] + i[1] for i \
                in zip(col1,col2)])
            for two_letter in combo_alpha:
                i,j = combo_alpha[two_letter]
                pmf[i,j] = two_letter_counts[two_letter]
            npmf = np.array(pmf,dtype=float) + epsilon
            npmf /= npmf.sum() 
            npmf += epsilon
            row_sum = npmf.sum(axis=1)
            col_sum = npmf.sum(axis=0)
            if np.any(row_sum > cnsrv_cut) \
                    or np.any(col_sum > cnsrv_cut):
                pass
            else:
                sum_ent = -(row_sum*np.log(row_sum)).sum() \
                        - (col_sum*np.log(col_sum)).sum()
                joint_ent = -(npmf*np.log(npmf)).sum()
                mi_matrix[pos1,pos2] = 2*(sum_ent \
                        - joint_ent)/joint_ent
    return mi_matrix + mi_matrix.T

def mi(alignment,alphabet=aminobet,epsilon=1e-9,gap_cut=.1,cnsrv_cut=.95):
    """
    Returns the mutual information matrix.
    Alphabet refers to the aminoacids to consider.
    Epsilon is there so the 0*log(0) is approximately 0.
    gap_cut refers to the frequency of gaps tolerated at scorable positions.
        gap_cut = 0 means that only positions with no gaps will be scored.
    consrv_cut refers the maximum frequency of a single aminoacid at scorable
        positions. cnsrv_cut = .95 means that if a column has any single aminoacid
        more than 95% of the time, the position will be ignored
        This is because calculation of NMI becomes unreliable at low entropy.
    """
    alpha_size = len(alphabet)

    trans_align = transpose_alignment(alignment)
    position_size = len(trans_align)
    seq_size = float(len(trans_align[0]))

    combo_alpha = {}
    for ind,i in enumerate(alphabet):
        for jnd,j in enumerate(alphabet):
            combo_alpha[i+j] = (ind,jnd)

    mi_matrix = np.zeros((position_size,position_size),dtype=float)

    for two_pos in comb(range(position_size),2):
        pos1 = two_pos[0]
        pos2 = two_pos[1]
        col1 = trans_align[pos1]
        col2 = trans_align[pos2]
        if ((col1.count("-")/seq_size) > gap_cut) \
                or ((col2.count("-")/seq_size) > gap_cut):
            pass
        else:
            pmf = np.zeros((alpha_size,alpha_size),dtype=int)
            two_letter_counts = Counter([i[0] + i[1] for i \
                in zip(col1,col2)])
            for two_letter in combo_alpha:
                i,j = combo_alpha[two_letter]
                pmf[i,j] = two_letter_counts[two_letter]
            npmf = np.array(pmf,dtype=float) + epsilon
            npmf /= npmf.sum() 
            npmf += epsilon
            row_sum = npmf.sum(axis=1)
            col_sum = npmf.sum(axis=0)
            if np.any(row_sum > cnsrv_cut) \
                    or np.any(col_sum > cnsrv_cut):
                pass
            else:
                sum_ent = -(row_sum*np.log(row_sum)).sum() \
                        - (col_sum*np.log(col_sum)).sum()
                joint_ent = -(npmf*np.log(npmf)).sum()
                mi_matrix[pos1,pos2] = (sum_ent - joint_ent)
    return mi_matrix + mi_matrix.T

def entropy(alignment,alphabet=aminobet,epsilon=1e-9,gap_cut=.1,cnsrv_cut=.95):
    """
    Returns the sum of entropies matrix.
    Alphabet refers to the aminoacids to consider.
    Epsilon is there so the 0*log(0) is approximately 0.
    gap_cut refers to the frequency of gaps tolerated at scorable positions.
        gap_cut = 0 means that only positions with no gaps will be scored.
    consrv_cut refers the maximum frequency of a single aminoacid at scorable
        positions. cnsrv_cut = .95 means that if a column has any single aminoacid
        more than 95% of the time, the position will be ignored
        This is because calculation of NMI becomes unreliable at low entropy.
    """
    alpha_size = len(alphabet)

    trans_align = transpose_alignment(alignment)
    position_size = len(trans_align)

    combo_alpha = {}
    for ind,i in enumerate(alphabet):
        for jnd,j in enumerate(alphabet):
            combo_alpha[i+j] = (ind,jnd)

    mi_matrix = np.zeros((position_size,position_size),dtype=float)

    for two_pos in comb(range(position_size),2):
        pos1 = two_pos[0]
        pos2 = two_pos[1]
        col1 = trans_align[pos1]
        col2 = trans_align[pos2]
        pmf = np.zeros((alpha_size,alpha_size),dtype=int)
        two_letter_counts = Counter([i[0] + i[1] for i \
            in zip(col1,col2)])
        for two_letter in combo_alpha:
            i,j = combo_alpha[two_letter]
            pmf[i,j] = two_letter_counts[two_letter]
        npmf = np.array(pmf,dtype=float) + epsilon
        npmf /= npmf.sum() 
        npmf += epsilon
        row_sum = npmf.sum(axis=1)
        col_sum = npmf.sum(axis=0)
        sum_ent = -(row_sum*np.log(row_sum)).sum() \
                - (col_sum*np.log(col_sum)).sum()
        mi_matrix[pos1,pos2] = (sum_ent)
    return mi_matrix + mi_matrix.T


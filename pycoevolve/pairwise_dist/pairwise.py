#/usr/bin/env python
from __future__ import division
from sets import Set
import numpy as np

def calcSimilarity(aln,threshold):
    """Create a ndarray representing the adjacency
    matrix, for a given threshold."""
    assert threshold < 1
    assert 0 <= threshold
    similarity = np.zeros([len(aln),len(aln)],dtype=int)
    perm_order = aln.keys()
    for ind,i in enumerate(perm_order):
        for jnd,j in enumerate(perm_order[ind:]):
            iseq = aln[i]
            jseq = aln[j]
            sim = 0
            for a,b, in zip(iseq,jseq):
                if a==b and a!='-':
                    sim += 1
            l = len(iseq)
            assert l == len(jseq)
            if sim/len(iseq) >= threshold:
                similarity[jnd,ind] = 1
                similarity[ind,jnd] = 1
    return similarity,perm_order

def calcContinuousSimilarity(aln):
    """
    Create a ndarray representing the
    Levehdistance between sequences.
    Returns the ndarray, the order used.
    """
    similarity = np.zeros([len(aln),len(aln)],dtype=float)
    perm_order = aln.keys()
    for ind,i in enumerate(perm_order):
        for jnd,j in enumerate(perm_order):
            iseq = aln[i]
            jseq = aln[j]
            sim = 0
            sublen = 0
            for a,b, in zip(iseq,jseq):
                if a==b:
                    if a!='-':
                        sim += 1
                    else:
                        sublen += 1
            l = float(len(iseq))
            usedLen = float(l - sublen)
            similarity[ind,jnd] = sim/usedLen
    return similarity,perm_order


def removeHigh(aln,threshold):
    """Removes sequences in the alignment greater than
    the threshold of percent sequence similarity for
    ungapped positions."""
    similarity,perm_order = calcSimilarity(aln,threshold)
    keep = Set(range(len(perm_order)))
    while np.any(similarity == 1):
        remove = np.where(similarity.sum(axis=1))[0][0]
        similarity[remove,:] = 0
        similarity[:,remove] = 0
        keep.discard(remove)
    keep_seqs = [perm_order[k] for k in keep]
    n_aln = {}
    for i in keep_seqs:
        n_aln[i] = aln[i]
    return n_aln





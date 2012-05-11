#/usr/bin/env python

def nongapIndexer(str):
    pos = []
    for ind,c in enumerate(str):
        if c != '-':
            pos.append(ind)
    return pos

def seqTruncer(str,pos):
    truncChars = []
    for ind in pos:
        truncChars.append(str[ind])
    return ''.join([t for t in truncChars])


def truncBySeqId(aln,seqId):
    """ Truncates a loaded alignment by the string for a
    regular expression. Fails for multiple matches."""
    truncTarget = aln[seqId]
    positions = nongapIndexer(truncTarget)
    truncAlign = {}
    for seq in aln:
        truncAlign[seq] = seqTruncer(aln[seq],positions)
    return truncAlign

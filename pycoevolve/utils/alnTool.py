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


def truncBySeqId(aln):
    """ Truncates a loaded alignment by the string for a
    regular expression. Fails for multiple matches."""
    matches = [k for k in aln.keys() if k.find("-") != -1]
    if len(matches) == 0:
        raise ValueError('No matches to regStr')
    if len(matches) > 1:
        n_mat = []
        for mat in matches:
            n_mat.append(''.join([i for i in mat if i != '-']))
        match = [n_mat.index(max(n_mat))]
    match = matches[0]
    truncTarget = aln[match]
    positions = nongapIndexer(truncTarget)
    truncAlign = {}
    for seq in aln:
        truncAlign[seq] = seqTruncer(aln[seq],positions)
    return truncAlign

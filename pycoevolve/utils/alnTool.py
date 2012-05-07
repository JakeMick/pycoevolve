#/usr/bin/env python
import re

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


def truncBySeqId(aln, regStr):
    """ Truncates a loaded alignment by the string for a
    regular expression. Fails for multiple matches."""
    reg = re.compile(regStr)
    matches = [k for k in aln.keys() if reg.match(k)]
    if len(matches) == 0:
        raise ValueError('No matches to regStr')
    if len(matches) > 1:
        raise ValueError('More than one match to regStr')
    match = matches[0]
    truncTarget = aln[match]
    positions = nongapIndexer(truncTarget)
    truncAlign = {}
    for seq in aln:
        truncAlign[seq] = seqTruncer(aln[seq],positions)
    return truncAlign

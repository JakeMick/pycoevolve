#/usr/bin/env python

from pycoevolve.toy_data import load_toy

def test_two_seqs():
    align = load_toy(subset="two_seqs")
    return align

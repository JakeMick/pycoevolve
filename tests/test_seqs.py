#/usr/bin/env python
from pycoevolve import seqs

def test_fasta():
    toy_cycled = "toy_data/cycled.fa"

    try:
        open(toy_cycled)
    except IOError:
        print "Need to run gen_toy_data first."

    align = seqs.fasta.load(toy_cycled,"r")
    


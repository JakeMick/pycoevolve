#/usr/bin/env python
import os
from pycoevolve.seqs import load_fasta_align


#The following three letter abbreviations
#           are prepended to toy data IDs:
###########################################################
#   -'filetype = 'abbreviations'
###########################################################
#   - fasta = seq
fast_join = lambda fasta : os.path.join("toy_fastas",fasta)
data_names = {
            "cyc_seq" : fast_join("cycled.fa"),
            "pmm_seq" : fast_join("pmm_out.fa"),
            "two_seq" : fast_join("two_seqs.fa")
            }

def load_toy_data(subset="two_seq"):
    file_handle = data_names[subset]

    if subset[-3:] == 'seq':
        return load_fasta_align(file_handle)



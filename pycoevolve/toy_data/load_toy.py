#/usr/bin/env python
from os.path import dirname, join

from pycoevolve.seqs.fasta import load_fasta_align

def _seq_join(fasta):
    toy_path = join(dirname(__file__),"toy_fastas")
    return join(toy_path,fasta)

#The following three letter abbreviations
#           are prepended to toy data IDs:
###########################################################
#   -'filetype = 'abbreviations'
###########################################################
#   - fasta = seq
#   - pdb   = pdb

data_names = {
            "cyc_seq" : _seq_join("cycled.fa"),
            "pmm_seq" : _seq_join("pmm_out.fa"),
            "two_seq" : _seq_join("two_seqs.fa")
            }

def load_toy_data(subset="two_seq"):
    """ Valid subsets are:
        cyc_seq : A 21x210 fasta of all of the standard
                    aminoalphabet and a gap.

        two_seq : The first two sequences from cyc_seq.

        pmm_seq : A MUSCLE alignment of phosphomannomutase
                    retrieved from the PIR on Apr 3, 2012.
    """

    file_handle = data_names[subset]

    if subset[-3:] == 'seq':
        return load_fasta_align(file_handle)



#/usr/bin/env python
from os.path import dirname, join, pardir
from pycoevolve.seqs.fasta import load_fasta_align

toy_fastas = {"cyc" : "cycled.fa",
            "pmm" : "pmm_out.fa",
            "two" : "two_seqs.fa"}

def test_load_fasta_align():
    """
    This tests load_fasta_align using toy_data's fastas.
    """
    aminobet = "-ACDEFGHIKLMNPQRSTVWY"
    def ring_pop(str):
        out = str[-1]
        str = str[:-1]
        return out + str


    seq_loc = join(dirname(__file__),pardir)
    mod_loc = join(seq_loc,pardir)
    toy_loc = join(mod_loc,"toy_data")
    fas_loc = join(toy_loc,"toy_fastas")

    #test cycled
    cyc_loc = join(fas_loc,toy_fastas["cyc"])
    cyc_aln = load_fasta_align(cyc_loc)
    assert len(cyc_aln) == 210
    assert sum([int(seq_id) for seq_id in cyc_aln]) == \
            sum([i for i in xrange(210)])
    for i in xrange(210):
        assert cyc_aln[str(i)] == aminobet
        aminobet = ring_pop(aminobet)

if __name__ == "__main__":
    test_load_fasta_align()

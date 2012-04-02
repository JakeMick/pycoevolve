#/usr/bin/env python

def make_cycled_align(n_times=10,filename="toy_data/cycled.fa"):
    """This function will create a fake fasta file that
    cycles through the aminoacid alphabet, where the
    sequence labels just numeric indices associated with
    the loop.

    n_times defines the time number of times to cycle through
    the aminoacid alphabets. n_times=10 results in 210 sequences

    filename is the default location of the file
    """

    def cycle_str(s):
        out = s[-1]
        return out + s[:-1]

    file_handle = open(filename,"w")
    aminobet = "-ACDEFGHIKLMNPQRSTVWY"

    for i in xrange(21*n_times):
        file_handle.write(">%i\n" % i)
        file_handle.write("%s\n" % aminobet)
        aminobet = cycle_str(aminobet)

    file_handle.close()

if __name__ == "__main__":
    make_cycled_align()



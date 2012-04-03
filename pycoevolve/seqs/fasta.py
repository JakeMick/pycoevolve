#/usr/bin/env python

def load_fasta_align(file_location,trunc_prot_name=True,
                                            trunc_len=6):
    """
    Creates a dict of the alignment where the keys are
    sequence name.

    trunc_prot_name = True
        Will truncate the sequence names.

    trunc_len = 6
        Determines the length of truncation.
    """
    file_handle = open(file_location,"r")

    align_list = [row.rstrip("\n") for row in file_handle]

    align = {}
    for row in align_list:
        if row[0] == ">":
            if trunc_prot_name:
                prot_name = row[1:trunc_len+1]
            else:
                prot_name = row[1:]
            align[prot_name] = []
        else:
            align[prot_name].append(row)
    for prot in align:
        align[prot] = "".join([row for row in align[prot]])

    return align


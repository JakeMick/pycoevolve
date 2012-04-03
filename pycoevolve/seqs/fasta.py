#/usr/bin/env python

def load(file_location,trunc_prot_name=True,trunc_len=6):
    """
    >>>from pycoevolve import toy_data
    >>>load(two_seqs.fa)
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
            align[prot_name] = ""
        else:
            align[prot_name] += row

    return align


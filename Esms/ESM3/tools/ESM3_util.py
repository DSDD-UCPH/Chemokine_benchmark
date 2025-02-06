import os
from Bio import SeqIO



def fasta_to_dict(path_fasta, dict_ID):
    info_dict = {}
    fasta_file_list = os.listdir(path_fasta)
    for ffile in fasta_file_list:
        if ffile.endswith(".DS_Store"):
            pass
        else:
            ffile_path = path_fasta + ffile
            for record in SeqIO.parse(ffile_path, "fasta"):
                uniprot = ffile.replace(".fasta", "")
                Gene_ID = dict_ID[uniprot]
                info_dict[Gene_ID] = str(record.seq)
    return info_dict
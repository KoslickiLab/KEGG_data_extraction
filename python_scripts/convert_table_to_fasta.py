import os
import sys
from os import listdir
from os.path import isfile, join
import pandas as pd

# This is the base directory on the GPU server
base_directory = "/data/shared_data/KEGG_data/"
#base_directory = "/home/dkoslicki/Documents/KEGG_data_extraction/data"
KEGG_prot_directory = os.path.join(base_directory, "organisms/kegg_gene_info")
#KEGG_prot_directory = os.path.join(base_directory, "input")
file_names = [os.path.join(KEGG_prot_directory, f) for f in listdir(KEGG_prot_directory) if isfile(join(KEGG_prot_directory, f))]

# Some of the sequences do not have KO's associated with them, so let's separate those out
aa_KO_out_file = os.path.join(base_directory, "output", "kegg_genes_KO.faa")
aa_NoKO_out_file = os.path.join(base_directory, "output", "kegg_genes_No_KO.faa")
nt_KO_out_file = os.path.join(base_directory, "output", "kegg_genes_KO.fna")
nt_NoKO_out_file = os.path.join(base_directory, "output", "kegg_genes_No_KO.fna")

with open(aa_KO_out_file, 'w') as aa_KO_fid:
    with open(aa_NoKO_out_file, 'w') as aa_NoKO_fid:
        with open(nt_KO_out_file, 'w') as nt_KO_fid:
            with open(nt_NoKO_out_file, 'w') as nt_NoKO_fid:
                for filename in file_names:
                    # read in one of the gene tables
                    df = pd.read_csv(filename, sep='\t', lineterminator='\n', header=0, keep_default_na=False)
                    df = df.reset_index()
                    for _, row in df.iterrows():
                        index, kegg_gene_id, desc, koid, aaseq, ntseq = row
                        # header will be the concatenated (with deliminter "|") kegg gene id, description, and koid
                        header = "|".join([kegg_gene_id, desc, koid])
                        if koid:
                            if aaseq:
                                aa_KO_fid.write(f">{header}\n{aaseq}\n")
                            if ntseq:
                                nt_KO_fid.write(f">{header}\n{ntseq}\n")
                        else:
                            if aaseq:
                                aa_NoKO_fid.write(f">{header}\n{aaseq}\n")
                            if ntseq:
                                nt_NoKO_fid.write(f">{header}\n{aaseq}\n")



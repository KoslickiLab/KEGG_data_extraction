import os
import sys
from os import listdir
from os.path import isfile, join
import pandas as pd

# This is the base directory on the GPU server
base_directory = "/data/shared_data/KEGG_data/"
KEGG_prot_directory = os.path.join(base_directory, "organisms/kegg_gene_info")
file_names = [os.path.join(KEGG_prot_directory, f) for f in listdir(KEGG_prot_directory) if isfile(join(KEGG_prot_directory, f))]
aa_out_file = os.path.join(base_directory, "kegg_genes.faa")
nt_out_file = os.path.join(base_directory, "kegg_genes.fna")
with open(aa_out_file, 'w') as aa_fid:
    with open(nt_out_file, 'w') as nt_fid:
        for filename in file_names:
            # read in one of the gene tables
            df = pd.read_csv(filename, sep='\t', lineterminator='\n', header=0, keep_default_na=False)
            df = df.reset_index()
            for _, row in df.iterrows():
                index, kegg_gene_id, desc, koid, aaseq, ntseq = row
                # header will be the concatenated (with deliminter "|") kegg gene id, description, and koid
                header = "|".join([kegg_gene_id, desc, koid])
                aa_fid.write(f">{header}\n{aaseq}\n")
                nt_fid.write(f">{header}\n{ntseq}\n")



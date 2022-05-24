#!/usr/bin/env python
import os
import sys
from os import listdir
from os.path import isfile, join
import pandas as pd
import argparse

def convert_table_to_FASTA(file_names, aa_KO_out_file, aa_NoKO_out_file, nt_KO_out_file, nt_NoKO_out_file):
    """
    Convert the tables of genes into FASTA sequences
    :param file_names: a list of all the KEGG gene tables to convert
    :param aa_KO_out_file: FASTA amino acid sequences that have associated KO IDs
    :param aa_NoKO_out_file: FASTA amino acid sequences that do NOY have associated KO IDs
    :param nt_KO_out_file: FASTA nucleotide sequences that have associated KO IDs
    :param nt_NoKO_out_file: FASTA nucleotide sequences that do NOY have associated KO IDs
    :return: None
    """
    with open(aa_KO_out_file, 'w') as aa_KO_fid:
        with open(aa_NoKO_out_file, 'w') as aa_NoKO_fid:
            with open(nt_KO_out_file, 'w') as nt_KO_fid:
                with open(nt_NoKO_out_file, 'w') as nt_NoKO_fid:
                    for filename in file_names:
                        print(f"converting file: {filename}")
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
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene_dir", type=str,
                        help="The full path of the directory that contains all of the .txt gene tables (and nothing else)",
                        default="/data/shared_data/KEGG_data/organisms/kegg_gene_info")
    parser.add_argument("--out_dir", type=str, help="The full path to the directory that the files will be written",
                        default="/data/shared_data/KEGG_data/")
    args = parser.parse_args()
    # parse the args
    KEGG_prot_directory = args.gene_dir
    if not os.path.exists(KEGG_prot_directory):
        raise Exception(f"Folder {KEGG_prot_directory} does not exist")
    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        print(f"Output folder {out_dir} does not exist, making one now.")
        os.makedirs(out_dir)
    # get the file names to convert
    file_names = [os.path.join(KEGG_prot_directory, f) for f in listdir(KEGG_prot_directory) if
                  isfile(join(KEGG_prot_directory, f))]
    # name the output files
    aa_KO_out_file = os.path.join(out_dir, "kegg_genes_KO.faa")
    aa_NoKO_out_file = os.path.join(out_dir, "kegg_genes_No_KO.faa")
    nt_KO_out_file = os.path.join(out_dir, "kegg_genes_KO.fna")
    nt_NoKO_out_file = os.path.join(out_dir, "kegg_genes_No_KO.fna")
    # then do the conversion
    convert_table_to_FASTA(file_names, aa_KO_out_file, aa_NoKO_out_file, nt_KO_out_file, nt_NoKO_out_file)
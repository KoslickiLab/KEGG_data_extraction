#!/bin/bash

## set up current path
here=$(pwd)

## extract organism data from kegg: organism table, organism associated genes
python ${here}/python_scripts/extract_kegg_organism_data.py --organisms 'Archaea' 'Bacteria' 'Fungi' --outdir ${here}/out_results/kegg_organisms

## extract virus data from kegg: virus table, virus associated genes
python ${here}/python_scripts/extract_kegg_virus_data.py --outdir ${here}/out_results/kegg_viruses

## download sequences
python ${here}/python_scripts/download_seq_fasta.py --table ${here}/out_results/kegg_viruses/virus_table.txt --col 'rs_ncbi_seq_ids' --outfile ${here}/out_results/kegg_viruses/rs_ncbi_virus.fasta
python ${here}/python_scripts/download_seq_fasta.py --table ${here}/out_results/kegg_viruses/virus_table.txt --col 'gb_ncbi_seq_id' --outfile ${here}/out_results/kegg_viruses/gb_ncbi_virus.fasta
python ${here}/python_scripts/download_seq_fasta.py --table ${here}/out_results/kegg_organisms/organism_table.txt --col 'rs_ncbi_seq_ids' --organisms 'Archaea' 'Bacteria' 'Fungi' --outfile ${here}/out_results/kegg_organisms/rs_ncbi_organism.fasta
python ${here}/python_scripts/download_seq_fasta.py --table ${here}/out_results/kegg_organisms/organism_table.txt --col 'gb_ncbi_seq_id' --organisms 'Archaea' 'Bacteria' 'Fungi' --outfile ${here}/out_results/kegg_organisms/gb_ncbi_organism.fasta

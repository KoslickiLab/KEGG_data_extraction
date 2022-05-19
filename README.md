# KEGG Data Extraction Scripts

## Description
This repository contains the python scripts used to download sequence data of organisms (e.g., Eukaryotes, Prokaryotes) or viruses from KEGG via its [APIs](https://www.kegg.jp/kegg/rest/keggapi.html). There is a shell script `main.sh` to summarize the steps to downalod KEGG sequence data.

## Implementation
You can simply `git clone` this repo to your local computer, and then run:
```shell
bash main.sh
```
Or

You can run the specific python scripts based on your need. There are three python scripts under `./python_scripts` folder:
##### extract_kegg_organism_data.py
This script is used to download the organism table and the associated origanisms' RefSeq and GeneBank genomes based on KEGG information. It has the following two parameters:

- \--organisms: Specify specific [organisms](http://rest.kegg.jp/list/organism) (e.g., 'Archaea' 'Bacteria' 'Fungi') for which you want to extract sequence information.
- \--outdir: Specify your output folder

Example: `python ${your_current_path}/python_scripts/extract_kegg_organism_data.py --organisms 'Archaea' 'Bacteria' 'Fungi' --outdir ${your_current_path}/out_results/kegg_organisms`

##### extract_kegg_virus_data.py
This script is used to download the viruses table and their associated RefSeq and GeneBank genomes based on KEGG information. It has only one parameter:

- \--outdir: Specify your output folder

Example: `python ${your_current_path}/python_scripts/extract_kegg_virus_data.py --outdir ${your_current_path}/out_results/kegg_viruses`

##### download_seq_fasta.py
This script needs to run after either/both of the above two scripts have been implemented. Iis used to download the gene sequences into a fasta-format file. It has only one parameter:

- \--table: Specify the full path of organisms/viruses table that is generated from the above two scripts.
- \--organisms: Specify specific [organisms](http://rest.kegg.jp/list/organism) (e.g., 'Archaea' 'Bacteria' 'Fungi') for which you want to extract sequence information.
- \--col: Specify which gene database `RefSeq` (e.g., "rs_ncbi_seq_ids") or `GeneBank` (e.g., "gb_ncbi_seq_id") you want to use
- \--outfile: Specify the full path of output file.

Examples: 
1. `python ${your_current_path}/python_scripts/download_seq_fasta.py --table ${here}/out_results/kegg_organisms/organism_table.txt --col 'rs_ncbi_seq_ids' --organisms 'Archaea' 'Bacteria' 'Fungi' --outfile ${your_current_path}/out_results/kegg_organisms/rs_ncbi_organism.fasta`

2. `python ${your_current_path}/python_scripts/download_seq_fasta.py --table ${here}/out_results/kegg_organisms/organism_table.txt --col 'gb_ncbi_seq_id' --organisms 'Archaea' 'Bacteria' 'Fungi' --outfile ${your_current_path}/out_results/kegg_organisms/gb_ncbi_organism.fasta`

## Data
You can find the data (only for Archaea' 'Bacteria' 'Fungi' and 'Viruses') that I have already downloaded previously from our GPU server. The data locates `/data/shared_data/KEGG_data`.
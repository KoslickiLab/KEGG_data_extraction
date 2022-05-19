
import os
import sys
import argparse
import logging
from tqdm import tqdm, trange
import pickle
from itertools import chain
import pandas as pd
import pytaxonkit
import requests
from bs4 import BeautifulSoup
from multiprocessing import Pool, cpu_count
import time
import re

def get_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s  [%(levelname)s]  %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

def extract_taxaid_seq(inlist):
    temp = '|'.join(inlist)
    if 'TAXONOMY' in temp:
        taxaid = [re.sub('\s.*','',re.sub('TAXONOMY\s*TAX:','',line)) for line in inlist if 'TAXONOMY' in line][0]
    else:
        taxaid = None
    if 'ORTHOLOGY' in temp:
        koid = 'ko:'+[re.sub('\s.*','',re.sub('ORTHOLOGY\s*','',line)) for line in inlist if 'ORTHOLOGY' in line][0]
    else:
        koid = None
    if 'AASEQ' in temp:
        aaseq = re.sub('\d*','','|'.join(inlist).split('AASEQ       ')[1].split('|COMMENT     ')[0].split('|NTSEQ     ')[0]).replace('|            ','').replace('|','')
    else:
        aaseq = None
    if 'NTSEQ' in temp:
        ntseq = re.sub('\d*','','|'.join(inlist).split('NTSEQ       ')[1].split('|COMMENT     ')[0]).replace('|            ','').replace('|','')
    else:
        ntseq = None
    return taxaid, koid, aaseq, ntseq

def process_query(instr):
    link = KEGG_api_link + f"/get/{instr}"
    res = requests.get(link)
    if res.status_code == 200:
        print(f"Successuflly extract info from {link}")
        res = [tuple([a])+b for a, b in zip(instr.split('+'),list(map(extract_taxaid_seq,[x.split('\n') for x in res.text.split('///')])))]
        return res
    else:
        print(f"Error: Fail to extract info from {link}")
        return []

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, help="The output dir")
    args = parser.parse_args()

    KEGG_api_link = 'http://rest.kegg.jp'

    logger = get_logger()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    ## download KEGG organism table
    if not os.path.exists(os.path.join(args.outdir,'virus_table.txt')):
        link = KEGG_api_link + '/get/br:br08620'
        res = requests.get(link)
        if res.status_code == 200:
            logger.info('Virus table is sucessfully downloaded!')
            temp_list = []
            for line in res.text.split('\n'):
                if 'TAX' in line:
                    line = re.sub('^\w\s*','',line).split('  ')[-1]
                    item1 = line.split(' [')[0] if len(line.split(' [')) != 0 else None
                    item2 = re.findall('TAX:\d+',line)[0] if len(re.findall('TAX:\d+',line)) != 0 else None
                    item3 = re.findall('RS:[a-zA-Z0-9_ ]*',line)[0] if len(re.findall('RS:[a-zA-Z0-9_ ]*',line)) != 0 else None
                    item4 = re.findall('GN:[a-zA-Z0-9_ ]*',line)[0] if len(re.findall('GN:[a-zA-Z0-9_ ]*',line)) != 0 else None
                    temp_list += [(item1,item2,item3,item4)]
            virus_table = pd.DataFrame(temp_list)
            virus_table.columns = ['name','taxaid','rs_ncbi_seq_ids','gb_ncbi_seq_id']
            virus_table['rs_ncbi_seq_ids'] = virus_table['rs_ncbi_seq_ids'].str.replace('RS:','').str.split(' ')
            virus_table['gb_ncbi_seq_id'] = virus_table['gb_ncbi_seq_id'].str.replace('GN:','').str.split(' ')
            virus_table['taxaid'] = virus_table.taxaid.str.replace('TAX:','')
            virus_table['taxaid'] = virus_table['taxaid'].astype(int)
            virus_table.to_csv(os.path.join(args.outdir,'virus_table.txt'),sep='\t',index=None)
    else:
        virus_table = pd.read_csv(os.path.join(args.outdir,'virus_table.txt'), sep='\t', header=0)
        virus_table.loc[virus_table['rs_ncbi_seq_ids'].isnull(),'rs_ncbi_seq_ids'] = 'None'
        virus_table['rs_ncbi_seq_ids'] = virus_table['rs_ncbi_seq_ids'].apply(lambda row: eval(row) if type(row) is str else row)
        virus_table.loc[virus_table['gb_ncbi_seq_id'].isnull(),'gb_ncbi_seq_id'] = 'None'
        virus_table['gb_ncbi_seq_id'] = virus_table['gb_ncbi_seq_id'].apply(lambda row: eval(row) if type(row) is str else row)
        virus_table['taxaid'] = virus_table['taxaid'].astype(int)
    
    ## extract gene/protein information from KEGG
    if not os.path.exists(os.path.join(args.outdir,'kegg_gene_info')):
        os.makedirs(os.path.join(args.outdir,'kegg_gene_info'))

    gene_table = pd.DataFrame(columns=['kegg_gene_id','desc'])
    for i in ['vg','vp']:
        link = KEGG_api_link + f'/list/{i}'
        res = requests.get(link)
        if res.status_code == 200:
            res.text.split('\n')
            temp = pd.DataFrame([line.split('\t') for line in res.text.split('\n') if line.split('\t')[0]])
            temp.columns = ['kegg_gene_id','desc']
            gene_table = pd.concat([gene_table,temp]).reset_index(drop=True)
        else:
            logger.error(f"Fail to access data from {link}")

    kegg_gene_id_list = list(gene_table['kegg_gene_id'])
    batch =list(range(0,len(kegg_gene_id_list),10))
    batch.append(len(kegg_gene_id_list))
    kegg_gene_id_list = ['+'.join(kegg_gene_id_list[batch[i-1]:batch[i]]) for i in range(1, len(batch))]
    batch =list(range(0,len(kegg_gene_id_list),20))
    batch.append(len(kegg_gene_id_list))

    final_res = []
    for i in trange(1, len(batch)):
        start = batch[i-1]
        end = batch[i]
        with Pool(processes=20) as excutator:
            res = excutator.map(process_query, kegg_gene_id_list[start:end])
        final_res += [y for x in res for y in x]
        time.sleep(5)
    final_res = pd.DataFrame(final_res)
    final_res.columns = ['kegg_gene_id','taxaid','koid','aaseq','ntseq']

    final_gene_table = final_res.merge(gene_table, on='kegg_gene_id').reset_index(drop=True)
    final_gene_table.to_csv(os.path.join(args.outdir,'kegg_gene_info','gene_table.txt'), sep='\t', index=None)

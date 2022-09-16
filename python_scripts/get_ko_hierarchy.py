import os
import sys
import argparse
import logging
from tqdm import tqdm, trange
import pickle
from itertools import chain
import pandas as pd
import requests
from bs4 import BeautifulSoup
import time
import re
from glob import glob
import json

def get_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s  [%(levelname)s]  %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

def organize_hierarchy(hiearchy_json, regex='^\w?\w?\d{5} ', prefix='', stop_level=None):

    def _iterate_multidimensional(res, hiearchy_json, res_list):
        if isinstance(hiearchy_json,dict):
            for k,v in hiearchy_json.items():
                if k == 'name':
                    temp_name = re.sub(' \[.*\]','',hiearchy_json['name'])
                    res += f"|{temp_name}"
                    if re.search(regex, hiearchy_json['name']) is not None:
                        res_list += [res]
                elif k == 'children':
                    for elem in hiearchy_json['children']:
                        _iterate_multidimensional(res, elem, res_list)
        else:
            self.logger.error(f"{hiearchy_json} is not dictionary")
            raise

    res_list = []
    _iterate_multidimensional('', hiearchy_json, res_list)
    if stop_level is None:
        res_dict = {prefix+string.split('|')[-1].split(' ')[0]:'|'.join(string.split('|')[1:]) for string in res_list}
        return res_dict
    else:
        res_dict = {prefix+string.split('|')[-1].split(' ')[0].split('\t')[0]:'|'.join(string.split('|')[1:(stop_level+1)]) for string in res_list}
        return res_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, help="The output dir")
    args = parser.parse_args()

    KEGG_api_link = 'http://rest.kegg.jp'

    logger = get_logger()

    ## get brite table
    link = f"{KEGG_api_link}/list/brite"
    res = requests.get(link)
    if res.status_code == 200:
        brite_table = pd.DataFrame([x.split('\t') for x in res.text.split('\n') if x.split('\t')[0]])
        brite_table.columns = ['kegg_brite_id','desc']
    else:
        logger.error(f"Fail to download KEGG brite information from {link}")
        exit()

    ## download KEGG KO associated hierarchy and process hierarchy
    logger.info(f"Download KEGG KO associated hierarchy")
    ko_hierarchy_dict = dict()
    for brite_id, desc in brite_table.to_numpy():
        logger.info(f"Processing brite id {brite_id}")
        check_link = f"{KEGG_api_link}/get/{brite_id}"
        res = requests.get(check_link)
        if res.status_code == 200:
            m = re.search('K\d{5}', res.text)
        else:
            logger.error(f"Fail to download KEGG brite information from {check_link}")
            continue
        if m is None:
            logger.warning(f"Brite ID {brite_id} doesn't contain KO ids and thus skip it.")
            continue
        link = f"{KEGG_api_link}/get/{brite_id}/json"
        json_res = requests.get(link)
        if json_res.status_code == 200:
            temp_dict = organize_hierarchy(json_res.json(), regex='^K\d{5} ')
            for key in temp_dict:
                if key in ko_hierarchy_dict:
                    ko_hierarchy_dict[key] += [temp_dict[key]]
                else:
                    ko_hierarchy_dict[key] = [temp_dict[key]]
        else:
            logger.error(f"Fail to download KEGG brite information from {link}")

    ## set up identifier mapping
    id_mapping = {f"{re.sub('^[a-z]*:[a-z]*','',x[0])} {x[1]}":x[0].split(':')[1] for x in brite_table.to_numpy()}

    ## convert hierarchy to edge list
    ko_edge_list = []
    for ko_id, hierarchy_list in ko_hierarchy_dict.items():
        for item in hierarchy_list:
            temp_list = [id_mapping.get(x,x) for x in item.split('|')[:-1]]
            temp_list += [ko_id]
            ko_edge_list += [(temp_list[index-1], temp_list[index]) for index in range(1, len(temp_list))]
                
    kegg_ko_edge_df = pd.DataFrame(ko_edge_list)
    kegg_ko_edge_df.columns = ['parent','child']
    kegg_ko_edge_df = kegg_ko_edge_df.drop_duplicates().reset_index(drop=True)
    kegg_ko_edge_df.to_csv(os.path.join(args.outdir, 'kegg_ko_edge_df.txt'), sep='\t', index=None)




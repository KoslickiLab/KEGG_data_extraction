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

    ## download KEGG KO brite hierarchy and process hierarchy
    logger.info(f"Download KEGG KO brite hierarchy")
    brite_id = 'br:ko00001'
    ko_hierarchy_dict = dict()
    link = f"{KEGG_api_link}/get/{brite_id}/json"
    res = requests.get(link)
    if res.status_code == 200:
        temp_dict = organize_hierarchy(res.json(), regex='^K\d{5} ')
        for key in temp_dict:
            if key in ko_hierarchy_dict:
                ko_hierarchy_dict[key] += [temp_dict[key]]
            else:
                ko_hierarchy_dict[key] = [temp_dict[key]]
    else:
        logger.error(f"Fail to download KEGG brite information from {link}")

    ## convert hierarchy to edge list
    ko_edge_list = []
    for ko_id, hierarchy_list in ko_hierarchy_dict.items():
        for item in hierarchy_list:
            temp_list = item.split('|')[:-1]
            temp_list += [ko_id]
            ko_edge_list += [(temp_list[index-1], temp_list[index]) for index in range(1, len(temp_list))]
                
    kegg_ko_edge_df = pd.DataFrame(ko_edge_list)
    kegg_ko_edge_df.columns = ['parent','child']
    kegg_ko_edge_df = kegg_ko_edge_df.drop_duplicates().reset_index(drop=True)
    kegg_ko_edge_df.to_csv(os.path.join(args.outdir, 'kegg_ko_edge_df.txt'), sep='\t', index=None)

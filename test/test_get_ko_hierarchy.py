import subprocess
import networkx as nx


def test_brite_subtree():
    cmd = "python ../python_scripts/get_ko_hierarchy.py --brite ko00001 --outdir test_data"
    res = subprocess.run(cmd, shell=True, check=True)
    assert res.returncode == 0
    edge_list_file = "test_data/kegg_ko_edge_df_br:ko00001.txt"
    # import as a networkx directed graph and check if it's a tree
    G = nx.read_edgelist(edge_list_file, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
    assert nx.is_tree(G)

import subprocess
import networkx as nx


def test_brite_subtree():
    cmd = "python ../python_scripts/get_ko_hierarchy.py --brite ko00001 --outdir test_data"
    res = subprocess.run(cmd, shell=True, check=True)
    assert res.returncode == 0
    edge_list_file = "test_data/kegg_ko_edge_df_br:ko00001.txt"
    # import as a networkx directed graph and check if it's a tree
    with open(edge_list_file, 'r') as fid:
        next(fid, '')
        G = nx.read_edgelist(fid, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
    if not nx.is_tree(G):
        print("Not a tree. Examples include:")
        itr = 0
        for node in G.nodes():
            ancestors = list(G.predecessors(node))
            if len(ancestors) > 1:
                print(f"{node}: {ancestors}")
                itr += 1
            if itr > 3:
                break
        print("Orphaned nodes are:")
        for node in G.nodes():
            if len(list(G.neighbors(node))) == 0:
                print(node)

    assert nx.is_tree(G)

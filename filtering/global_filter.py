import networkx as nx

def filter_consensus_matrix(edge_list_path, threshold):
    fh=open(edge_list_path, 'r')
    G=nx.read_weighted_edgelist(fh, nodetype=str)
    fh.close()
    
    # if weight < threshold, remove edge
    for node in G:
        for neighbor in G[node]:
            if G[node][neighbor]['weight'] < threshold:
                G.remove_edge(node, neighbor)
    
    nx.write_edgelist(G, f"globally_filtered_{edge_list_path}")
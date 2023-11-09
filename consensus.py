import networkx as nx
import community as community_louvain
import numpy as np
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
from sklearn.metrics.cluster import normalized_mutual_info_score


def apply_clustering_algorithm(G, clustering_algo, nP):
    partitions = []
    for _ in range(nP):
        partition = clustering_algo(G)
        # partition = [list(community) for community in partition]
        partitions.append(partition)
    return partitions

def compute_nmi(p1,p2):
    labels1 = list(p1.values())
    labels2 = list(p2.values())
    nmi = normalized_mutual_info_score(labels1, labels2)
    return nmi

def compute_consensus_matrix(partitions, n):
    consensus_matrix = np.zeros((n, n))
    nodes=[p for p in partitions[0].keys()]
    nodes.sort()

    # for partition in partitions:
        # for cluster in partition:
            # for i in range(n):
                # for j in range(n):
                    # if i in cluster and j in cluster:
                        # consensus_matrix[i][j] += 1

    node_to_index = {node: i for i, node in enumerate(nodes)}
    for partition in partitions:
        for n1, n2 in itertools.combinations(nodes, 2):
            if partition[n1] == partition[n2]:
                i, j = node_to_index[n1], node_to_index[n2]
                consensus_matrix[i][j] += 1
                consensus_matrix[j][i] += 1
                
    consensus_matrix /= len(partitions)
    return consensus_matrix, nodes 

def filter_consensus_matrix(consensus_matrix, threshold):
    consensus_matrix[consensus_matrix < threshold] = 0
    return consensus_matrix

def find_consensus_partition(G, clustering_algo, nP, threshold,name):
    
    n = len(G.nodes)
    partitions = apply_clustering_algorithm(G, clustering_algo, nP)
    
    with open(f'initial_{name}_cluster.txt', 'w') as f:
        for community in set(partitions[0].values()):
            f.write("%s\n" % ' '.join([node for node, c in partitions[0].items() if c == community]))
    
    ### 1. draw initial partitioned graph
         
    # parts=[[k for k, v in partitions[0].items() if v == value] for value in set(partitions[0].values())]
    # parts.sort(key=len, reverse=True)
    
    ## only draw top 10 largest clusters
    
    # global nodes_to_draw
    # nodes_to_draw = [node for cluster in parts[:10] for node in cluster]
    # G_sub = G.subgraph(nodes_to_draw)
    # G_sub = nx.Graph(G_sub)
    # G_sub.remove_edges_from(nx.selfloop_edges(G_sub))
    # color_map = [partitions[0][node] for node in G_sub.nodes()]
    
    ## draw all clusters
    
    # color_map = []
    # for node in G:
        # color_map.append(partitions[0][node])
    # most_common = [item[0] for item in Counter(color_map).most_common(3)]
    # color_map = [x if x in most_common else max(color_map)+1 for x in color_map]
    # G.remove_edges_from(nx.selfloop_edges(G))
    # pos = nx.spring_layout(G)
    # pos=nx.kamada_kawai_layout(G)
    # nx.draw(G,pos, node_color=color_map, with_labels=False, node_size=5, width=0.1, cmap=plt.cm.jet)
    # plt.savefig("initial_communities.png", dpi=1000)
    # plt.show()
    # plt.close()
    
    print("Initial partitions completed")
    
    global nodes
    consensus_matrix, nodes = compute_consensus_matrix(partitions, n)
    consensus_matrix = filter_consensus_matrix(consensus_matrix, threshold)
    
    i=0
    while True:
        
        print(f"iteration {i}")
        new_partitions = apply_clustering_algorithm(nx.Graph(consensus_matrix), clustering_algo, nP)
        
        p1 = set(tuple(k for k, v in new_partitions[-1].items() if v == value) for value in set(new_partitions[-1].values()))
        p2 = set(tuple(k for k, v in new_partitions[-2].items() if v == value) for value in set(new_partitions[-2].values()))
        print(f"Length of last 2 partitions: {len(p1)}, {len(p2)}")
        print(f"length of intersection:{len(p1 & p2)}")
        print(f"NMI = {compute_nmi(new_partitions[-1],new_partitions[-2])}")
        if p1 == p2:
            break
        else:
            consensus_matrix,_ = compute_consensus_matrix(new_partitions, n)
            
            ### 2. Draw heatmap for consensus_cluster to check if it is a block diagonal matrix
            
            # new_order = []
            # for cluster in p2:
            #     new_order.extend(cluster)
            # consensus_matrix_ordered = consensus_matrix[new_order, :]
            # heatmap = sns.heatmap(consensus_matrix_ordered, cmap="YlGnBu")
            # heatmap.set(xlabel='Node', ylabel='Node')
            # plt.savefig(f"consensus_matrix{i}.png", dpi=1000)
            # plt.clf()
            # plt.close()
            
        i+=1
    
    ### 3. Draw heatmap for final consensus matrix
    
    # new_order = []
    # for cluster in p2:
    #     new_order.extend(cluster)
    # consensus_matrix = consensus_matrix[new_order, :]
    # heatmap = sns.heatmap(consensus_matrix, cmap="YlGnBu")
    # heatmap.set(xlabel='Node', ylabel='Node')
    # plt.savefig("consensus_matrix_final.png", dpi=1000)
    # plt.close()
    
    ### 4. Draw NMI (final consensus compared each of the intial partitions) value in term of partition number
    
    # NMIs=[]
    # for partition in partitions:
    #     NMIs.append(compute_nmi(partition,new_partitions[-1]))
    # plot NMI value in term of partition number
    # plt.plot(range(1,11),NMIs)
    # plt.xlabel('Partition number')
    # plt.ylabel('NMI value')
    # plt.savefig("NMI_value.png", dpi=1000)
    # plt.close()
    print("Consensus partitions completed\n")
    return new_partitions[-1]


def run_consensus_alg(edge_list_path, nP, threshold, output_name):
    G = nx.read_weighted_edgelist(edge_list_path, nodetype=str)
    n = len(G.nodes)   
    consensus_partition = find_consensus_partition(G, community_louvain.best_partition, nP, threshold,output_name)
    with open(f'consensus_{output_name}_cluster.txt', 'w') as f:
        for community in set(consensus_partition.values()):
            f.write("%s\n" % ' '.join([nodes[node] for node, c in consensus_partition.items() if c == community]))
  
run_consensus_alg('../final_layers/reduced_mirna_backbone.txt',10,0.1,'mirna')          
# run_consensus_alg('../final_layers/ppi_final_net_selfless.txt',10,0.1,'ppi')
# run_consensus_alg('../final_layers/reduced_averaged_tumor_exp_backbone_0.02.txt',10,0.1,'exp')
# run_consensus_alg('../final_layers/reduced_tf_backbone.txt',10,0.1,'tf')
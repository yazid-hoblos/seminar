import networkx as nx
import community as community_louvain
import numpy as np
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
from sklearn.metrics.cluster import normalized_mutual_info_score
import os

def apply_clustering_algorithm(G, clustering_algo, nP):
    partitions = []
    for _ in range(nP):
        partition = clustering_algo(G)
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

    node_to_index = {node: i for i, node in enumerate(nodes)}
    for partition in partitions:
        for n1, n2 in itertools.combinations(nodes, 2):
            if n1 in partition and n2 in partition and partition[n1] == partition[n2]:
                i, j = node_to_index[n1], node_to_index[n2]
                consensus_matrix[i][j] += 1
                consensus_matrix[j][i] += 1
                
    consensus_matrix /= len(partitions)
    return consensus_matrix, nodes 

def filter_consensus_matrix(consensus_matrix, threshold):
    consensus_matrix[consensus_matrix < threshold] = 0
    return consensus_matrix

def read_partition(partition_dir):
    partitions = {}
    for file in os.listdir(partition_dir):
        partition={}
        i=0
        with open(os.path.join(partition_dir, file), 'r') as f:
            for line in f:
                for node in line.strip().split():
                    partition[node] = i
                i+=1
            if 'exp' in file:
                partitions['exp'] = partition
            elif 'ppi' in file:
                partitions['ppi'] = partition
            elif 'mirna' in file:
                partitions['mirna'] = partition
            else:
                partitions['tf'] = partition
    return partitions 

def find_consensus_partition(extracted_partitions,clustering_algo, n, nP, threshold,name):
    
    labels = ['exp', 'ppi', 'mirna', 'tf']
    # labels = ['exp', 'ppi', 'tf']
    partitions = [extracted_partitions[p] for p in labels]
    
    print("Layer partitions Opened")
    
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
        i+=1
    
    return new_partitions[-1]


def run_consensus_alg(layer_partitions_dir, n, nP, threshold, output_name=None):
    layer_partitions = read_partition(layer_partitions_dir)
    consensus_partition = find_consensus_partition(layer_partitions,community_louvain.best_partition, n, nP, threshold,output_name)
    with open(f'consensus_cluster_across_layers.txt', 'w') as f:
        for community in set(consensus_partition.values()):
            f.write("%s\n" % ' '.join([nodes[node] for node, c in consensus_partition.items() if c == community]))
            
run_consensus_alg("../final_partitions", 3924, 10, 0.1)
import networkx as nx
import matplotlib.pyplot as plt
import os

def plot_dist(edge_list_path):
    G=nx.read_weighted_edgelist(edge_list_path, nodetype=str)
    fig, ax = plt.subplots()
    ax.hist([d for n, d in G.degree()], bins=20, color='skyblue', edgecolor='black')
    ax.set_xlabel('Degree', fontsize=15)
    ax.set_ylabel('Frequency', fontsize=15)
    ax.set_title('Degree Distribution', fontsize=15)
    ax.grid(True)
    
    base_name = os.path.basename(edge_list_path)
    file_name_without_ext = os.path.splitext(base_name)[0]
    plt.savefig(f"{file_name_without_ext}_dist.png", dpi=1000)
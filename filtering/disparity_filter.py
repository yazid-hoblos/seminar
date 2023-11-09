
import networkx as nx
from scipy import integrate

def extract_backbone(g, alpha):
  backbone_graph = nx.Graph()
  print("Extracting backbone...")
  i=0
  for node in g:
    if i%1000==0:
        print(f"i {nodes} completed")
    k_n = len(g[node])
    if k_n > 1:
        sum_w = sum( g[node][neighbor]['weight'] for neighbor in g[node] )
        for neighbor in g[node]:
            edgeWeight = g[node][neighbor]['weight']
            pij = float(edgeWeight)/sum_w
            f = lambda x: (1-x)**(k_n-2) 
            alpha_ij =  1 - (k_n-1)*integrate.quad(f, 0, pij)[0] 
            if alpha_ij < alpha: 
                backbone_graph.add_edge(node,neighbor, weight = edgeWeight)
    i+=1
  return backbone_graph


def disparity_filter(edge_list_path,alpha):
  fh=open(edge_list_path, 'r')
  G=nx.read_weighted_edgelist(fh, nodetype=str)
  fh.close()
  nx.write_edgelist(extract_backbone(G, alpha), f"backbone_{alpha}_{edge_list_path}")

  # for i,a in enumerate(alpha):
  #   nx.write_edgelist(extract_backbone(G, a), "backbone" + str(i) + ".txt")
    #print(extract_backbone(G, a))

# vet=[0.005,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.4,0.5]
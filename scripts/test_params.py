from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np
from src.features.structure_parser import ProteinStructureParser
from src.features.distance_matrix import compute_distance_matrix
from src.features.sequential_graph_builder import build_sequential_knn_graph, build_contact_map_graph
from src.models.louvain_clustering import louvain_clustering, get_cluster_sizes

# Test on one protein
parser = ProteinStructureParser()
coords, _ = parser.parse_structure('6p5a', 'B')
D = compute_distance_matrix(coords)

print("Testing Sequential KNN Graph")
print("True domains: 1\n")

for k in [3, 5, 7, 10]:
    for seq_w in [2.0, 3.0, 5.0]:
        for res in [0.3, 0.5]:
            G = build_sequential_knn_graph(D, k=k, seq_weight=seq_w)
            partition = louvain_clustering(G, resolution=res)
            n_pred = len(get_cluster_sizes(partition))
            
            error = abs(n_pred - 1)
            print(f"k={k:2d}, seq_w={seq_w:.1f}, res={res:.1f} -> {n_pred:2d} domains (error={error})")

print("\n" + "="*60)
print("Testing Contact Map Graph")
print("="*60 + "\n")

for threshold in [6.0, 8.0, 10.0]:
    for res in [0.3, 0.5, 0.7]:
        G = build_contact_map_graph(D, threshold=threshold)
        partition = louvain_clustering(G, resolution=res)
        n_pred = len(get_cluster_sizes(partition))
        
        error = abs(n_pred - 1)
        print(f"threshold={threshold:4.1f}, res={res:.1f} -> {n_pred:2d} domains (error={error})")
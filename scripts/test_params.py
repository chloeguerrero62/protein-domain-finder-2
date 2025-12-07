from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np
from src.features.structure_parser import ProteinStructureParser
from src.features.distance_matrix import compute_distance_matrix
from src.features.graph_builder import build_knn_graph
from src.models.louvain_clustering import louvain_clustering, get_cluster_sizes

# Test on one protein
parser = ProteinStructureParser()
coords, _ = parser.parse_structure('6p5a', 'B')  # 1 domain protein
D = compute_distance_matrix(coords)

print("True domains: 1\n")

for k in [3, 5, 7, 10]:
    for res in [0.3, 0.5, 0.7, 1.0]:
        G = build_knn_graph(D, k=k)
        partition = louvain_clustering(G, resolution=res)
        n_pred = len(get_cluster_sizes(partition))
        
        error = abs(n_pred - 1)
        print(f"k={k:2d}, res={res:.1f} → {n_pred:2d} domains (error={error})")
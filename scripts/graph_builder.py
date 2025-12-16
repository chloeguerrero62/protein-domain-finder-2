# scripts/graph_builder.py (UPDATED)

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import networkx as nx
import pickle  # Add this

from src.features.graph_builder import build_knn_graph

# Setup
output_dir = Path('data/graphs')
output_dir.mkdir(parents=True, exist_ok=True)

# Load distance matrices
dm_dir = Path('data/distance_matrices')

for dm_file in dm_dir.glob('*_distmat.npy'):
    print(f"Processing {dm_file.stem}")
    
    # Load distance matrix
    D = np.load(dm_file)
    
    # Build graph
    G = build_knn_graph(D, k=10)
    
    # Save graph using pickle (modern way)
    output_file = output_dir / f"{dm_file.stem}_graph.pkl"
    with open(output_file, 'wb') as f:
        pickle.dump(G, f)
    
    print(f"Saved graph with {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

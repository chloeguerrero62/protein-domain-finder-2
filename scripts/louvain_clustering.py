"""
Apply Louvain clustering to all protein graphs
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import pickle

from src.models.louvain_clustering import (
    louvain_clustering,
    partition_to_labels,
    get_cluster_sizes,
    get_modularity
)


def main():
    # Setup
    graph_dir = Path('data/graphs')
    output_dir = Path('data/clusters')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    graph_files = list(graph_dir.glob('*_graph.pkl'))
    
    print(f"Clustering {len(graph_files)} proteins...\n")
    
    results = []
    
    for i, graph_file in enumerate(graph_files, 1):
        pdb_chain = graph_file.stem.replace('_graph', '')
        print(f"[{i}/{len(graph_files)}] {pdb_chain}...", end=' ')
        
        try:
            # Load graph
            with open(graph_file, 'rb') as f:
                G = pickle.load(f)
            
            # Run Louvain clustering
            partition = louvain_clustering(G, resolution=1.0, random_state=42)
            
            # Convert to labels array
            labels = partition_to_labels(partition)
            
            # Get statistics
            cluster_sizes = get_cluster_sizes(partition)
            n_clusters = len(cluster_sizes)
            modularity = get_modularity(G, partition)
            
            # Save partition
            output_file = output_dir / f"{pdb_chain}_partition.npy"
            np.save(output_file, labels)
            
            # Save partition dict as well (for debugging)
            partition_file = output_dir / f"{pdb_chain}_partition.pkl"
            with open(partition_file, 'wb') as f:
                pickle.dump(partition, f)
            
            results.append({
                'pdb_chain': pdb_chain,
                'n_residues': G.number_of_nodes(),
                'n_domains': n_clusters,
                'modularity': modularity,
                'cluster_sizes': str(cluster_sizes),
                'success': True
            })
            
            print(f"✓ {n_clusters} domains (Q={modularity:.3f})")
            
        except Exception as e:
            results.append({
                'pdb_chain': pdb_chain,
                'n_residues': 0,
                'n_domains': 0,
                'modularity': 0.0,
                'cluster_sizes': '',
                'success': False,
                'error': str(e)
            })
            print(f"✗ {e}")
    
    # Summary
    print(f"\n{'='*60}")
    successful = sum(1 for r in results if r['success'])
    print(f"Success: {successful}/{len(graph_files)}")
    
    if successful > 0:
        avg_domains = np.mean([r['n_domains'] for r in results if r['success']])
        avg_modularity = np.mean([r['modularity'] for r in results if r['success']])
        print(f"Average domains per protein: {avg_domains:.1f}")
        print(f"Average modularity: {avg_modularity:.3f}")
    
    print(f"{'='*60}")
    
    # Save results
    results_df = pd.DataFrame(results)
    results_file = Path('data/results/clustering_results.csv')
    results_file.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(results_file, index=False)
    
    print(f"\nResults saved to: {results_file}")
    print(f"Partitions saved to: {output_dir}/")


if __name__ == "__main__":
    main()
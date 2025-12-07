"""
Run complete domain detection pipeline on all proteins
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np
from tqdm import tqdm

from src.features.structure_parser import ProteinStructureParser
from src.features.distance_matrix import compute_distance_matrix
from src.features.sequential_graph_builder import build_sequential_knn_graph
from src.models.louvain_clustering import (
    louvain_clustering,
    partition_to_labels,
    get_cluster_sizes,
    get_modularity
)


def process_single_protein(row, parser, k=5, seq_weight=3.0, resolution=0.5):
    """Process one protein through entire pipeline"""
    
    pdb_id = row['pdb_id']
    chain_id = row['chain']
    
    try:
        coords, res_ids = parser.parse_structure(pdb_id, chain_id)
        D = compute_distance_matrix(coords)
        G = build_sequential_knn_graph(D, k=k, seq_weight=seq_weight, weighting='inverse')
        partition = louvain_clustering(G, resolution=resolution, random_state=42)
        labels = partition_to_labels(partition)
        cluster_sizes = get_cluster_sizes(partition)
        n_predicted = len(cluster_sizes)
        modularity = get_modularity(G, partition)
        
        return {
            'pdb_chain': f"{pdb_id}_{chain_id}",
            'pdb_id': pdb_id,
            'chain': chain_id,
            'n_residues': len(coords),
            'n_predicted_domains': n_predicted,
            'n_true_domains': row['n_domains'],
            'modularity': modularity,
            'cluster_sizes': str(cluster_sizes),
            'success': True,
            'error': None
        }
        
    except Exception as e:
        return {
            'pdb_chain': f"{pdb_id}_{chain_id}",
            'pdb_id': pdb_id,
            'chain': chain_id,
            'n_residues': 0,
            'n_predicted_domains': 0,
            'n_true_domains': row.get('n_domains', 0),
            'modularity': 0.0,
            'cluster_sizes': '',
            'success': False,
            'error': str(e)
        }


def main():
    df = pd.read_csv('data/pdb_clustering/selected_proteins_mmseqs2.csv')
    print(f"Processing {len(df)} proteins with sequential KNN graph\n")
    
    parser = ProteinStructureParser(pdb_dir='data/selected_structures')
    
    results = []
    
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing proteins"):
        result = process_single_protein(row, parser, k=5, seq_weight=3.0, resolution=0.5)
        results.append(result)
    
    results_df = pd.DataFrame(results)
    
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    
    successful = results_df['success'].sum()
    print(f"Successfully processed: {successful}/{len(df)}")
    
    if successful > 0:
        correct = results_df[results_df['success']]
        
        exact_match = (correct['n_predicted_domains'] == correct['n_true_domains']).sum()
        mae = (correct['n_predicted_domains'] - correct['n_true_domains']).abs().mean()
        avg_modularity = correct['modularity'].mean()
        
        print(f"\nAccuracy:")
        print(f"  Exact matches: {exact_match}/{successful} ({100*exact_match/successful:.1f}%)")
        print(f"  Mean Absolute Error: {mae:.2f} domains")
        print(f"  Average modularity: {avg_modularity:.3f}")
        
        print(f"\nPredicted domain distribution:")
        print(correct['n_predicted_domains'].value_counts().sort_index())
        
        print(f"\nTrue domain distribution:")
        print(correct['n_true_domains'].value_counts().sort_index())
    
    if len(results_df) - successful > 0:
        print(f"\nFailed proteins:")
        failed = results_df[~results_df['success']]
        for _, row in failed.iterrows():
            print(f"  {row['pdb_chain']}: {row['error']}")
    
    Path('data/results').mkdir(parents=True, exist_ok=True)
    results_df.to_csv('data/results/pipeline_results_sequential.csv', index=False)
    print(f"\nResults saved to: data/results/pipeline_results_sequential.csv")


if __name__ == "__main__":
    main()
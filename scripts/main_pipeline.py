"""
Run complete domain detection pipeline on all proteins
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np
from tqdm import tqdm  # pip install tqdm if needed

# ✅ CORRECT IMPORTS (functions, not classes)
from src.features.structure_parser import ProteinStructureParser  # ← This IS a class
from src.features.distance_matrix import compute_distance_matrix  # ← Function
from src.features.graph_builder import build_knn_graph  # ← Function
from src.models.louvain_clustering import (
    louvain_clustering,
    partition_to_labels,
    get_cluster_sizes,
    get_modularity
)


def process_single_protein(row, parser):
    """Process one protein through entire pipeline"""
    
    pdb_id = row['pdb_id']
    chain_id = row['chain']
    
    try:
        # Step 1: Parse structure
        coords, res_ids = parser.parse_structure(pdb_id, chain_id)
        
        # Step 2: Compute distance matrix
        D = compute_distance_matrix(coords)
        
        # Step 3: Build graph
        G = build_knn_graph(D, k=5, weighting='inverse')
        
        # Step 4: Detect domains using Louvain
        partition = louvain_clustering(G, resolution=0.7, random_state=42)
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
    # Load dataset
    df = pd.read_csv('data/pdb_clustering/selected_proteins_mmseqs2.csv')
    print(f"Processing {len(df)} proteins...\n")
    
    # Initialize parser (only one that's a class)
    parser = ProteinStructureParser(pdb_dir='data/selected_structures')
    
    # Process all proteins
    results = []
    
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing proteins"):
        result = process_single_protein(row, parser)
        results.append(result)
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    
    successful = results_df['success'].sum()
    print(f"Successfully processed: {successful}/{len(df)}")
    
    if successful > 0:
        correct = results_df[results_df['success']]
        
        # Exact match accuracy
        exact_match = (correct['n_predicted_domains'] == correct['n_true_domains']).sum()
        
        # Mean Absolute Error
        mae = (correct['n_predicted_domains'] - correct['n_true_domains']).abs().mean()
        
        # Average modularity
        avg_modularity = correct['modularity'].mean()
        
        print(f"\nAccuracy:")
        print(f"  Exact matches: {exact_match}/{successful} ({100*exact_match/successful:.1f}%)")
        print(f"  Mean Absolute Error: {mae:.2f} domains")
        print(f"  Average modularity: {avg_modularity:.3f}")
        
        print(f"\nPredicted domain distribution:")
        print(correct['n_predicted_domains'].value_counts().sort_index())
        
        print(f"\nTrue domain distribution:")
        print(correct['n_true_domains'].value_counts().sort_index())
    
    # Failed proteins
    if len(results_df) - successful > 0:
        print(f"\nFailed proteins:")
        failed = results_df[~results_df['success']]
        for _, row in failed.iterrows():
            print(f"  {row['pdb_chain']}: {row['error']}")
    
    # Save results
    Path('data/results').mkdir(parents=True, exist_ok=True)
    results_df.to_csv('data/results/pipeline_results.csv', index=False)
    print(f"\nResults saved to: data/results/pipeline_results.csv")


if __name__ == "__main__":
    main()
"""
Comprehensive comparison of all clustering methods
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np
from tqdm import tqdm

from src.features.structure_parser import ProteinStructureParser
from src.features.distance_matrix import compute_distance_matrix
from src.features.graph_builder import build_knn_graph
from src.evaluation.clustering_comparison import (
    apply_louvain,
    apply_spectral_distance,
    apply_spectral_graph,
    apply_hierarchical,
    apply_two_stage_spectral
)
from src.evaluation.metrics import compute_all_metrics


def process_protein_all_methods(row, parser):
    """Apply all clustering methods to one protein"""
    
    pdb_id = row['pdb_id']
    chain_id = row['chain']
    n_true = row['n_domains']
    
    try:
        coords, _ = parser.parse_structure(pdb_id, chain_id)
        D = compute_distance_matrix(coords)
        G = build_knn_graph(D, k=10)
        
        methods_results = []
        
        # Method 1: Louvain (unsupervised)
        labels = apply_louvain(D, G)
        metrics = compute_all_metrics(labels, None, n_true)
        metrics['method'] = 'Louvain'
        metrics['supervised'] = False
        methods_results.append(metrics)
        
        # Method 2: Spectral on distance matrix (supervised)
        labels = apply_spectral_distance(D, n_true)
        metrics = compute_all_metrics(labels, None, n_true)
        metrics['method'] = 'Spectral-Distance'
        metrics['supervised'] = True
        methods_results.append(metrics)
        
        # Method 3: Spectral on graph (supervised)
        labels = apply_spectral_graph(G, n_true)
        metrics = compute_all_metrics(labels, None, n_true)
        metrics['method'] = 'Spectral-Graph'
        metrics['supervised'] = True
        methods_results.append(metrics)
        
        # Method 4: Hierarchical (supervised)
        labels = apply_hierarchical(D, n_true)
        metrics = compute_all_metrics(labels, None, n_true)
        metrics['method'] = 'Hierarchical'
        metrics['supervised'] = True
        methods_results.append(metrics)
        
        # Method 5: Two-stage spectral (unsupervised)
        labels, n_estimated = apply_two_stage_spectral(D, G)
        metrics = compute_all_metrics(labels, None, n_true)
        metrics['method'] = 'Two-Stage-Spectral'
        metrics['supervised'] = False
        metrics['n_estimated'] = n_estimated
        methods_results.append(metrics)
        
        # Add common info
        for m in methods_results:
            m['pdb_chain'] = f"{pdb_id}_{chain_id}"
            m['n_residues'] = len(coords)
            m['success'] = True
        
        return methods_results
        
    except Exception as e:
        return [{
            'pdb_chain': f"{pdb_id}_{chain_id}",
            'method': 'ALL',
            'success': False,
            'error': str(e)
        }]


def main():
    df = pd.read_csv('data/pdb_clustering/selected_proteins_mmseqs2.csv')
    parser = ProteinStructureParser()
    
    print("Comparing all clustering methods on all proteins\n")
    
    all_results = []
    
    for _, row in tqdm(df.iterrows(), total=len(df)):
        results = process_protein_all_methods(row, parser)
        all_results.extend(results)
    
    results_df = pd.DataFrame(all_results)
    
    # Summary by method
    print("\n" + "="*70)
    print("SUMMARY BY METHOD")
    print("="*70)
    
    successful = results_df[results_df['success']].copy()
    
    for method in successful['method'].unique():
        method_df = successful[successful['method'] == method]
        exact = (method_df['exact_match'] == 1).sum()
        mae = method_df['absolute_error'].mean()
        
        print(f"\n{method}:")
        print(f"  Exact matches: {exact}/{len(method_df)} ({100*exact/len(method_df):.1f}%)")
        print(f"  Mean Absolute Error: {mae:.2f}")
        
        if 'ari' in method_df.columns and method_df['ari'].notna().any():
            print(f"  Mean ARI: {method_df['ari'].mean():.3f}")
            print(f"  Mean NMI: {method_df['nmi'].mean():.3f}")
    
    # Save results
    Path('data/results').mkdir(parents=True, exist_ok=True)
    results_df.to_csv('data/results/method_comparison.csv', index=False)
    
    # Create summary table
    summary = successful.groupby('method').agg({
        'exact_match': 'sum',
        'absolute_error': 'mean',
        'n_predicted': 'mean'
    }).round(2)
    summary.to_csv('data/results/method_summary.csv')
    
    print(f"\nDetailed results: data/results/method_comparison.csv")
    print(f"Summary table: data/results/method_summary.csv")


if __name__ == "__main__":
    main()

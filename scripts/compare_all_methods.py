"""
Comprehensive comparison of clustering methods for protein domain detection

Methods compared:
1. Louvain (unsupervised) - discovers domain count automatically
2. Spectral-Graph (supervised) - requires oracle knowledge of true domain count
3. Two-Stage-Spectral (unsupervised) - estimates domain count via silhouette score
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
    apply_spectral_graph,
    apply_two_stage_spectral
)
from src.evaluation.metrics import compute_all_metrics


def process_protein_all_methods(row, parser):
    """Apply all three clustering methods to one protein"""
    
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
        
        # Method 2: Spectral on graph (supervised - uses oracle n_domains)
        labels = apply_spectral_graph(G, n_true)
        metrics = compute_all_metrics(labels, None, n_true)
        metrics['method'] = 'Spectral-Graph'
        metrics['supervised'] = True
        methods_results.append(metrics)
        
        # Method 3: Two-stage spectral (unsupervised)
        labels, n_estimated = apply_two_stage_spectral(D, G)
        metrics = compute_all_metrics(labels, None, n_true)
        metrics['method'] = 'Two-Stage-Spectral'
        metrics['supervised'] = False
        metrics['n_estimated'] = n_estimated
        methods_results.append(metrics)
        
        # Add common info to all results
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
    
    print("="*70)
    print("COMPARING CLUSTERING METHODS")
    print("="*70)
    print("\nMethods:")
    print("  1. Louvain (unsupervised)")
    print("  2. Spectral-Graph (supervised - oracle n_domains)")
    print("  3. Two-Stage-Spectral (unsupervised)")
    print(f"\nDataset: {len(df)} proteins\n")
    
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
    
    for method in ['Louvain', 'Spectral-Graph', 'Two-Stage-Spectral']:
        method_df = successful[successful['method'] == method]
        if len(method_df) == 0:
            continue
            
        exact = (method_df['exact_match'] == 1).sum()
        mae = method_df['absolute_error'].mean()
        supervised = method_df['supervised'].iloc[0]
        
        print(f"\n{method} ({'supervised' if supervised else 'unsupervised'}):")
        print(f"  Exact matches: {exact}/{len(method_df)} ({100*exact/len(method_df):.1f}%)")
        print(f"  Mean Absolute Error: {mae:.2f}")
        print(f"  Mean predicted domains: {method_df['n_predicted'].mean():.2f}")
    
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
    
    print(f"\n" + "="*70)
    print("FILES SAVED")
    print("="*70)
    print(f"Detailed results: data/results/method_comparison.csv")
    print(f"Summary table: data/results/method_summary.csv")


if __name__ == "__main__":
    main()
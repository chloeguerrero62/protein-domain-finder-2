"""
Apply Louvain clustering to all protein graphs

This script generates predictions WITHOUT knowing the true number of domains.
Louvain is expected to OVER-CLUSTER severely (predicting ~10 domains when true is ~2.5)
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm

from src.models.louvain_clustering import (
    louvain_clustering,
    partition_to_labels,
    get_cluster_sizes,
    get_modularity
)


def main():
    # Load dataset to get ground truth
    dataset_file = Path('data/pdb_clustering/selected_proteins_mmseqs2.csv')
    if not dataset_file.exists():
        print(f"ERROR: Dataset file not found: {dataset_file}")
        print("Cannot compute accuracy without ground truth.")
        print("Running without evaluation...")
        df = None
    else:
        df = pd.read_csv(dataset_file)
        # Create pdb_chain column if not exists
        if 'pdb_chain' not in df.columns:
            df['pdb_chain'] = df['pdb_id'] + '_' + df['chain']
        df = df.set_index('pdb_chain')
    
    # Setup
    graph_dir = Path('data/graphs')
    output_dir = Path('data/clusters')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    graph_files = list(graph_dir.glob('*_graph.pkl'))
    
    print("="*70)
    print("LOUVAIN CLUSTERING (UNSUPERVISED)")
    print("="*70)
    print(f"Processing {len(graph_files)} proteins")
    print(f"Resolution: 1.0")
    print(f"Expected: SEVERE over-clustering (~10 domains vs true ~2.5)")
    print()
    
    results = []
    
    for graph_file in tqdm(graph_files, desc="Clustering"):
        # Extract pdb_chain from filename
        # Format: pdb_id_chain_distmat_graph.pkl
        filename = graph_file.stem  # e.g., "6p5a_B_distmat_graph"
        # Remove _distmat_graph suffix
        pdb_chain = filename.replace('_distmat_graph', '')
        
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
            n_predicted = len(cluster_sizes)
            modularity = get_modularity(G, partition)
            
            # Save partition
            output_file = output_dir / f"{pdb_chain}_partition.npy"
            np.save(output_file, labels)
            
            # Save partition dict as well (for debugging)
            partition_file = output_dir / f"{pdb_chain}_partition.pkl"
            with open(partition_file, 'wb') as f:
                pickle.dump(partition, f)
            
            # Get ground truth if available
            n_true = None
            exact_match = None
            absolute_error = None
            relative_error = None
            
            if df is not None and pdb_chain in df.index:
                n_true = int(df.loc[pdb_chain, 'n_domains'])
                exact_match = int(n_predicted == n_true)
                absolute_error = abs(n_predicted - n_true)
                relative_error = absolute_error / n_true if n_true > 0 else 0
            
            results.append({
                'pdb_chain': pdb_chain,
                'n_residues': G.number_of_nodes(),
                'n_predicted': n_predicted,
                'n_true': n_true,
                'exact_match': exact_match,
                'absolute_error': absolute_error,
                'relative_error': relative_error,
                'modularity': modularity,
                'cluster_sizes': str(cluster_sizes),
                'success': True
            })
            
        except Exception as e:
            results.append({
                'pdb_chain': pdb_chain,
                'n_residues': 0,
                'n_predicted': 0,
                'n_true': None,
                'exact_match': None,
                'absolute_error': None,
                'relative_error': None,
                'modularity': 0.0,
                'cluster_sizes': '',
                'success': False,
                'error': str(e)
            })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Summary statistics
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    successful = results_df[results_df['success'] == True]
    failed = results_df[results_df['success'] == False]
    
    print(f"Successful: {len(successful)}/{len(results_df)}")
    print(f"Failed: {len(failed)}/{len(results_df)}")
    
    if len(successful) > 0:
        print(f"\nPredicted domain counts:")
        print(f"  Mean: {successful['n_predicted'].mean():.2f}")
        print(f"  Median: {successful['n_predicted'].median():.0f}")
        print(f"  Range: {successful['n_predicted'].min():.0f} - {successful['n_predicted'].max():.0f}")
        print(f"  Std Dev: {successful['n_predicted'].std():.2f}")
        
        print(f"\nModularity:")
        print(f"  Mean: {successful['modularity'].mean():.3f}")
        print(f"  Median: {successful['modularity'].median():.3f}")
        
        # If ground truth available
        if successful['n_true'].notna().any():
            with_truth = successful[successful['n_true'].notna()].copy()
            
            print(f"\n" + "="*70)
            print("ACCURACY (vs Pfam ground truth)")
            print("="*70)
            
            exact = with_truth['exact_match'].sum()
            print(f"Exact matches: {exact}/{len(with_truth)} ({100*exact/len(with_truth):.1f}%)")
            print(f"Mean Absolute Error: {with_truth['absolute_error'].mean():.2f}")
            print(f"Median Absolute Error: {with_truth['absolute_error'].median():.2f}")
            print(f"Std Deviation: {with_truth['absolute_error'].std():.2f}")
            print(f"Mean Relative Error: {with_truth['relative_error'].mean():.2%}")
            
            print(f"\nComparison:")
            print(f"  True domain count (mean): {with_truth['n_true'].mean():.2f}")
            print(f"  Predicted domain count (mean): {with_truth['n_predicted'].mean():.2f}")
            print(f"  Over-prediction factor: {with_truth['n_predicted'].mean() / with_truth['n_true'].mean():.2f}x")
            
            # Distribution comparison
            print(f"\nError distribution:")
            error_bins = pd.cut(with_truth['absolute_error'], 
                              bins=[-0.5, 0.5, 2.5, 5.5, 10.5, 100],
                              labels=['0', '1-2', '3-5', '6-10', '>10'])
            print(error_bins.value_counts().sort_index())
            
            # Show worst predictions
            print(f"\nWorst 5 predictions (highest error):")
            worst = with_truth.nlargest(5, 'absolute_error')[
                ['pdb_chain', 'n_true', 'n_predicted', 'absolute_error', 'n_residues']
            ]
            print(worst.to_string(index=False))
            
    
    if len(failed) > 0:
        print(f"\nFailed proteins:")
        for _, row in failed.iterrows():
            print(f"  {row['pdb_chain']}: {row.get('error', 'Unknown error')}")
    
    # Save results
    results_file = Path('data/results/louvain_results.csv')
    results_file.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(results_file, index=False)
    
    print(f"\n" + "="*70)
    print("FILES SAVED")
    print("="*70)
    print(f"Results: {results_file}")
    print(f"Partitions: {output_dir}/*_partition.npy")
    print(f"Partition dicts: {output_dir}/*_partition.pkl")


if __name__ == "__main__":
    main()
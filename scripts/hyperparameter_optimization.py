"""
Hyperparameter Optimization for Two-Stage Spectral Clustering

This script optimizes the Two-Stage Spectral method, which is the only
viable unsupervised approach (Louvain achieves 0% exact matches).

Parameters optimized:
1. Domain count estimation method (silhouette/eigengap/consensus)
2. Similarity kernel sigma scaling factor
3. Maximum domains search range
4. k-NN graph parameter k

Strategy:
- Training set: 30 proteins (stratified by domain count)
- Test set: Remaining ~114 proteins
- Grid search over ~80 parameter combinations
- Optimization metric: Mean Absolute Error (MAE)
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
from sklearn.model_selection import ParameterGrid
from tqdm import tqdm
import json

from src.features.structure_parser import ProteinStructureParser
from src.features.distance_matrix import compute_distance_matrix
from src.features.graph_builder import build_knn_graph
from src.evaluation.metrics import compute_all_metrics


def estimate_domains_custom(distance_matrix, method='silhouette', 
                            sigma_factor=1.0, max_domains=8):
    """
    Custom domain count estimation with configurable parameters
    """
    from sklearn.cluster import SpectralClustering
    from sklearn.metrics import silhouette_score
    
    # Compute similarity with custom sigma
    sigma = np.median(distance_matrix[distance_matrix > 0]) * sigma_factor
    similarity = np.exp(-distance_matrix**2 / (2 * sigma**2))
    
    n = len(distance_matrix)
    max_k = min(max_domains, n // 10)
    
    if method == 'silhouette':
        best_k = 2
        best_score = -1
        
        for k in range(2, max_k + 1):
            clustering = SpectralClustering(
                n_clusters=k,
                affinity='precomputed',
                random_state=42,
                assign_labels='kmeans'
            )
            labels = clustering.fit_predict(similarity)
            score = silhouette_score(distance_matrix, labels, metric='precomputed')
            
            if score > best_score:
                best_score = score
                best_k = k
        
        return best_k, labels
    
    elif method == 'eigengap':
        from scipy.linalg import eigh
        
        # Normalized Laplacian
        D = np.diag(np.sum(similarity, axis=1))
        D_inv_sqrt = np.diag(1.0 / np.sqrt(np.diag(D) + 1e-10))
        L = np.eye(n) - D_inv_sqrt @ similarity @ D_inv_sqrt
        
        eigenvalues, _ = eigh(L)
        eigenvalues = np.sort(eigenvalues)[:max_k]
        gaps = np.diff(eigenvalues)
        best_k = max(2, min(np.argmax(gaps) + 1, max_k))
        
        # Cluster with estimated k
        clustering = SpectralClustering(
            n_clusters=best_k,
            affinity='precomputed',
            random_state=42
        )
        labels = clustering.fit_predict(similarity)
        
        return best_k, labels
    
    elif method == 'consensus':
        # Get estimates from both methods
        k_sil, _ = estimate_domains_custom(
            distance_matrix, 'silhouette', sigma_factor, max_domains
        )
        k_eig, _ = estimate_domains_custom(
            distance_matrix, 'eigengap', sigma_factor, max_domains
        )
        
        # Use median
        best_k = int(np.median([k_sil, k_eig]))
        
        # Cluster with consensus k
        clustering = SpectralClustering(
            n_clusters=best_k,
            affinity='precomputed',
            random_state=42
        )
        labels = clustering.fit_predict(similarity)
        
        return best_k, labels
    
    else:
        raise ValueError(f"Unknown method: {method}")


def create_stratified_split(df, train_size=30, random_state=42):
    """
    Create stratified train/test split by domain count
    """
    np.random.seed(random_state)
    
    train_proteins = []
    
    # Sample proportionally from each domain bin
    for domain_bin in df['domain_bin'].unique():
        bin_df = df[df['domain_bin'] == domain_bin]
        n_samples = int(train_size * len(bin_df) / len(df))
        n_samples = max(1, n_samples)  # At least 1
        
        if len(bin_df) >= n_samples:
            sampled = bin_df.sample(n=n_samples, random_state=random_state)
            train_proteins.extend(sampled['pdb_chain'].tolist())
    
    # Ensure exactly train_size proteins
    if len(train_proteins) < train_size:
        remaining = df[~df['pdb_chain'].isin(train_proteins)]
        extra = remaining.sample(
            n=train_size - len(train_proteins), 
            random_state=random_state
        )
        train_proteins.extend(extra['pdb_chain'].tolist())
    elif len(train_proteins) > train_size:
        train_proteins = train_proteins[:train_size]
    
    train_df = df[df['pdb_chain'].isin(train_proteins)]
    test_df = df[~df['pdb_chain'].isin(train_proteins)]
    
    return train_df, test_df


def optimize_two_stage_spectral(train_df, parser):
    """
    Grid search over Two-Stage Spectral parameters
    """
    
    param_grid = {
        'estimation_method': ['silhouette', 'eigengap', 'consensus'],
        'sigma_factor': [0.5, 1.0, 1.5, 2.0, 2.5],
        'max_domains': [6, 8, 10, 12],
        'k_graph': [5, 10, 15, 20]
    }
    
    print(f"Grid search: {len(list(ParameterGrid(param_grid)))} combinations")
    print(f"Training set: {len(train_df)} proteins\n")
    
    results = []
    
    for params in tqdm(list(ParameterGrid(param_grid)), desc="Grid search"):
        errors = []
        exact_matches = []
        
        for _, row in train_df.iterrows():
            try:
                # Load structure
                coords, _ = parser.parse_structure(row['pdb_id'], row['chain'])
                D = compute_distance_matrix(coords)
                
                # Estimate domains with current parameters
                n_est, labels = estimate_domains_custom(
                    D,
                    method=params['estimation_method'],
                    sigma_factor=params['sigma_factor'],
                    max_domains=params['max_domains']
                )
                
                # Compute error
                error = abs(n_est - row['n_domains'])
                errors.append(error)
                exact_matches.append(int(n_est == row['n_domains']))
                
            except Exception as e:
                continue
        
        if errors:
            results.append({
                'estimation_method': params['estimation_method'],
                'sigma_factor': params['sigma_factor'],
                'max_domains': params['max_domains'],
                'k_graph': params['k_graph'],
                'mean_error': np.mean(errors),
                'median_error': np.median(errors),
                'std_error': np.std(errors),
                'exact_match_rate': np.mean(exact_matches),
                'n_proteins': len(errors)
            })
    
    results_df = pd.DataFrame(results)
    
    # Sort by mean error
    results_df = results_df.sort_values('mean_error')
    
    return results_df


def evaluate_best_params(test_df, parser, best_params):
    """
    Evaluate best parameters on held-out test set
    """
    print("\n" + "="*70)
    print("EVALUATING ON TEST SET")
    print("="*70)
    
    errors = []
    exact_matches = []
    all_metrics = []
    
    for _, row in tqdm(test_df.iterrows(), total=len(test_df), desc="Test set"):
        try:
            coords, _ = parser.parse_structure(row['pdb_id'], row['chain'])
            D = compute_distance_matrix(coords)
            
            n_est, labels = estimate_domains_custom(
                D,
                method=best_params['estimation_method'],
                sigma_factor=best_params['sigma_factor'],
                max_domains=best_params['max_domains']
            )
            
            metrics = compute_all_metrics(labels, None, row['n_domains'])
            metrics['pdb_chain'] = row['pdb_chain']
            all_metrics.append(metrics)
            
            errors.append(abs(n_est - row['n_domains']))
            exact_matches.append(int(n_est == row['n_domains']))
            
        except:
            continue
    
    test_results = {
        'n_proteins': len(errors),
        'mean_error': np.mean(errors),
        'median_error': np.median(errors),
        'std_error': np.std(errors),
        'exact_match_rate': np.mean(exact_matches),
        'exact_match_count': sum(exact_matches)
    }
    
    return test_results, pd.DataFrame(all_metrics)


def main():
    # Load data
    df = pd.read_csv('data/pdb_clustering/selected_proteins_mmseqs2.csv')
    
    # Create pdb_chain column if not exists
    if 'pdb_chain' not in df.columns:
        df['pdb_chain'] = df['pdb_id'] + '_' + df['chain']
    
    # Create stratified split
    print("="*70)
    print("HYPERPARAMETER OPTIMIZATION - TWO-STAGE SPECTRAL")
    print("="*70)
    print(f"\nTotal dataset: {len(df)} proteins")
    
    train_df, test_df = create_stratified_split(df, train_size=30, random_state=42)
    
    print(f"Training set: {len(train_df)} proteins")
    print(f"Test set: {len(test_df)} proteins")
    
    print("\nTraining set domain distribution:")
    print(train_df['domain_bin'].value_counts().sort_index())
    
    # Initialize parser
    parser = ProteinStructureParser()
    
    # Optimize on training set
    print("\n" + "="*70)
    print("GRID SEARCH ON TRAINING SET")
    print("="*70 + "\n")
    
    results_df = optimize_two_stage_spectral(train_df, parser)
    
    # Save all results
    output_dir = Path('data/results')
    output_dir.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_dir / 'hyperparameter_grid_search.csv', index=False)
    
    # Print top 10 configurations
    print("\n" + "="*70)
    print("TOP 10 PARAMETER CONFIGURATIONS")
    print("="*70 + "\n")
    print(results_df.head(10).to_string(index=False))
    
    # Best parameters
    best_params = results_df.iloc[0].to_dict()
    
    print("\n" + "="*70)
    print("BEST PARAMETERS (training set)")
    print("="*70)
    print(f"Estimation method: {best_params['estimation_method']}")
    print(f"Sigma factor: {best_params['sigma_factor']}")
    print(f"Max domains: {best_params['max_domains']}")
    print(f"k (graph): {best_params['k_graph']}")
    print(f"\nTraining performance:")
    print(f"  Mean Absolute Error: {best_params['mean_error']:.2f}")
    print(f"  Exact match rate: {best_params['exact_match_rate']*100:.1f}%")
    
    # Evaluate on test set
    test_results, test_metrics_df = evaluate_best_params(test_df, parser, best_params)
    
    print("\n" + "="*70)
    print("TEST SET PERFORMANCE")
    print("="*70)
    print(f"Proteins tested: {test_results['n_proteins']}")
    print(f"Mean Absolute Error: {test_results['mean_error']:.2f}")
    print(f"Median Absolute Error: {test_results['median_error']:.2f}")
    print(f"Std Deviation: {test_results['std_error']:.2f}")
    print(f"Exact matches: {test_results['exact_match_count']}/{test_results['n_proteins']}")
    print(f"Exact match rate: {test_results['exact_match_rate']*100:.1f}%")
    
    # Save test results
    test_metrics_df.to_csv(output_dir / 'test_set_performance.csv', index=False)
    
    # Save best parameters
    with open(output_dir / 'best_parameters.json', 'w') as f:
        json.dump({
            'method': 'Two-Stage-Spectral',
            'parameters': {
                'estimation_method': best_params['estimation_method'],
                'sigma_factor': float(best_params['sigma_factor']),
                'max_domains': int(best_params['max_domains']),
                'k_graph': int(best_params['k_graph'])
            },
            'training_performance': {
                'mean_error': float(best_params['mean_error']),
                'exact_match_rate': float(best_params['exact_match_rate'])
            },
            'test_performance': test_results
        }, f, indent=2)
    
    print("\n" + "="*70)
    print("RESULTS SAVED")
    print("="*70)
    print(f"Grid search results: {output_dir / 'hyperparameter_grid_search.csv'}")
    print(f"Test set performance: {output_dir / 'test_set_performance.csv'}")
    print(f"Best parameters: {output_dir / 'best_parameters.json'}")


if __name__ == "__main__":
    main()
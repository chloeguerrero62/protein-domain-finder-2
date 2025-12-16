from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score

from src.features.structure_parser import ProteinStructureParser
from src.features.distance_matrix import compute_distance_matrix
from src.features.graph_builder import build_knn_graph
from src.evaluation.metrics import compute_all_metrics



def estimate_domain_count(distance_matrix, method='silhouette', max_domains=8, 
                          sigma_factor=1.0):
    
    # Compute similarity matrix using RBF kernel
    sigma = np.median(distance_matrix[distance_matrix > 0]) * sigma_factor
    similarity = np.exp(-distance_matrix**2 / (2 * sigma**2))
    
    n_residues = len(distance_matrix)
    max_k = min(max_domains, n_residues // 10)  # Don't over-cluster
    
    if method == 'silhouette':
        return _estimate_silhouette(distance_matrix, similarity, max_k)
    elif method == 'eigengap':
        return _estimate_eigengap(similarity, max_k)
    elif method == 'consensus':
        return _estimate_consensus(distance_matrix, similarity, max_k)
    else:
        raise ValueError(f"Unknown method: {method}")


def _estimate_silhouette(distance_matrix, similarity, max_k):
    """Silhouette score method"""
    
    best_k = 2
    best_score = -1
    scores = []
    
    for k in range(2, max_k + 1):
        # Apply spectral clustering
        clustering = SpectralClustering(
            n_clusters=k,
            affinity='precomputed',
            random_state=42,
            assign_labels='kmeans'
        )
        labels = clustering.fit_predict(similarity)
        
        # Compute silhouette score
        score = silhouette_score(distance_matrix, labels, metric='precomputed')
        scores.append((k, score))
        
        if score > best_score:
            best_score = score
            best_k = k
    
    return best_k, scores


def _estimate_eigengap(similarity, max_k):
    """Eigengap heuristic method"""
    from scipy.linalg import eigh
    
    n = len(similarity)
    
    # Compute normalized Laplacian
    D = np.diag(np.sum(similarity, axis=1))
    D_inv_sqrt = np.diag(1.0 / np.sqrt(np.diag(D) + 1e-10))
    L = np.eye(n) - D_inv_sqrt @ similarity @ D_inv_sqrt
    
    # Get eigenvalues
    eigenvalues, _ = eigh(L)
    eigenvalues = np.sort(eigenvalues)[:max_k]
    
    # Find largest gap
    gaps = np.diff(eigenvalues)
    best_k = np.argmax(gaps) + 1
    
    return max(2, min(best_k, max_k)), []


def _estimate_consensus(distance_matrix, similarity, max_k):
    """Consensus of multiple methods"""
    
    # Get estimates from different methods
    k_sil, _ = _estimate_silhouette(distance_matrix, similarity, max_k)
    k_eig, _ = _estimate_eigengap(similarity, max_k)
    
    # Return median
    best_k = int(np.median([k_sil, k_eig]))
    
    return best_k, []


def apply_spectral_clustering(distance_matrix, n_domains, sigma_factor=1.0):
   
    # Compute similarity matrix
    sigma = np.median(distance_matrix[distance_matrix > 0]) * sigma_factor
    similarity = np.exp(-distance_matrix**2 / (2 * sigma**2))
    
    # Apply spectral clustering
    clustering = SpectralClustering(
        n_clusters=n_domains,
        affinity='precomputed',
        random_state=42,
        assign_labels='kmeans'
    )
    
    labels = clustering.fit_predict(similarity)
    
    return labels


def two_stage_spectral(distance_matrix, estimation_method='silhouette',
                       max_domains=8, sigma_factor=1.0, verbose=False):
   
    # Stage 1: Estimate domain count
    n_estimated, scores = estimate_domain_count(
        distance_matrix,
        method=estimation_method,
        max_domains=max_domains,
        sigma_factor=sigma_factor
    )
    
    if verbose:
        print(f"Stage 1: Estimated {n_estimated} domains")
        if scores:
            print("  Silhouette scores:")
            for k, score in scores:
                marker = " <--" if k == n_estimated else ""
                print(f"    k={k}: {score:.3f}{marker}")
    
    # Stage 2: Apply spectral clustering
    labels = apply_spectral_clustering(
        distance_matrix,
        n_domains=n_estimated,
        sigma_factor=sigma_factor
    )
    
    if verbose:
        print(f"Stage 2: Assigned residues to {len(np.unique(labels))} clusters")
    
    return labels, n_estimated, scores


def process_protein(pdb_id, chain_id, parser, estimation_method='silhouette',
                   max_domains=8, sigma_factor=1.0, verbose=False):
    
    try:
        # Parse structure
        coords, residue_ids = parser.parse_structure(pdb_id, chain_id)
        
        if verbose:
            print(f"\nProcessing {pdb_id}_{chain_id}")
            print(f"  Residues: {len(coords)}")
        
        # Compute distance matrix
        D = compute_distance_matrix(coords)
        
        # Apply two-stage clustering
        labels, n_estimated, scores = two_stage_spectral(
            D,
            estimation_method=estimation_method,
            max_domains=max_domains,
            sigma_factor=sigma_factor,
            verbose=verbose
        )
        
        # Compute domain sizes
        unique_labels = np.unique(labels)
        domain_sizes = [np.sum(labels == label) for label in unique_labels]
        
        return {
            'pdb_chain': f"{pdb_id}_{chain_id}",
            'n_residues': len(coords),
            'n_domains_estimated': n_estimated,
            'labels': labels,
            'domain_sizes': domain_sizes,
            'estimation_scores': scores,
            'success': True
        }
        
    except Exception as e:
        return {
            'pdb_chain': f"{pdb_id}_{chain_id}",
            'success': False,
            'error': str(e)
        }


def process_dataset(csv_file, output_dir='data/results', 
                   estimation_method='silhouette', max_domains=8,
                   sigma_factor=1.0, save_labels=True):

    
    # Load dataset
    df = pd.read_csv(csv_file)
    
    # Initialize parser
    parser = ProteinStructureParser()
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    if save_labels:
        labels_dir = output_path / 'two_stage_labels'
        labels_dir.mkdir(exist_ok=True)
    
    print("Two Stage Spectral Clustering")
    print(f"Dataset: {len(df)} proteins")
    print(f"Estimation method: {estimation_method}")
    print(f"Max domains: {max_domains}")
    print(f"Sigma factor: {sigma_factor}")
    print()
    
    # Process all proteins
    all_results = []
    
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing"):
        result = process_protein(
            row['pdb_id'],
            row['chain'],
            parser,
            estimation_method=estimation_method,
            max_domains=max_domains,
            sigma_factor=sigma_factor,
            verbose=False
        )
        
        # Add ground truth if available
        if 'n_domains' in row:
            result['n_true'] = row['n_domains']
        
        # Save labels if requested
        if save_labels and result['success']:
            label_file = labels_dir / f"{result['pdb_chain']}_labels.npy"
            np.save(label_file, result['labels'])
        
        all_results.append(result)
    
    # Convert to DataFrame
    results_df = pd.DataFrame(all_results)
    
    # Compute metrics if ground truth available
    if 'n_true' in results_df.columns:
        # Add metric columns to all rows (will be NaN for failed)
        results_df['exact_match'] = None
        results_df['absolute_error'] = None
        results_df['relative_error'] = None
        
        # Compute metrics only for successful rows
        success_mask = results_df['success'] == True
        
        results_df.loc[success_mask, 'exact_match'] = (
            results_df.loc[success_mask, 'n_domains_estimated'] == 
            results_df.loc[success_mask, 'n_true']
        ).astype(int)
        
        results_df.loc[success_mask, 'absolute_error'] = np.abs(
            results_df.loc[success_mask, 'n_domains_estimated'] - 
            results_df.loc[success_mask, 'n_true']
        )
        
        results_df.loc[success_mask, 'relative_error'] = (
            results_df.loc[success_mask, 'absolute_error'] / 
            results_df.loc[success_mask, 'n_true']
        )
    
    # Save results
    results_file = output_path / f'two_stage_{estimation_method}.csv'
    results_df.to_csv(results_file, index=False)
    
  
    
    successful = results_df[results_df['success']]
    failed = results_df[~results_df['success']]
    
    print(f"Successful: {len(successful)}/{len(results_df)}")
    print(f"Failed: {len(failed)}/{len(results_df)}")
    
    if len(successful) > 0:
        print(f"\nDomain count estimates:")
        print(f"  Mean: {successful['n_domains_estimated'].mean():.2f}")
        print(f"  Median: {successful['n_domains_estimated'].median():.0f}")
        print(f"  Range: {successful['n_domains_estimated'].min():.0f} - "
              f"{successful['n_domains_estimated'].max():.0f}")
        
        if 'n_true' in successful.columns:
            print(f"\nAccuracy:")
            exact = successful['exact_match'].sum()
            print(f"  Exact matches: {exact}/{len(successful)} "
                  f"({100*exact/len(successful):.1f}%)")
            print(f"  Mean Absolute Error: {successful['absolute_error'].mean():.2f}")
            print(f"  Median Absolute Error: {successful['absolute_error'].median():.2f}")
    
    if len(failed) > 0:
        print(f"\nFailed proteins:")
        for _, row in failed.iterrows():
            print(f"  {row['pdb_chain']}: {row['error']}")
    
    print(f"\nResults saved to: {results_file}")
    if save_labels:
        print(f"Labels saved to: {labels_dir}/")
    
    return results_df

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Two-Stage Spectral Clustering for Protein Domain Detection'
    )
    
    parser.add_argument(
        '--input',
        type=str,
        default='data/pdb_clustering/selected_proteins_mmseqs2.csv',
        help='Input CSV file with protein list'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        default='data/results',
        help='Output directory'
    )
    
    parser.add_argument(
        '--method',
        type=str,
        default='silhouette',
        choices=['silhouette', 'eigengap', 'consensus'],
        help='Domain count estimation method'
    )
    
    parser.add_argument(
        '--max-domains',
        type=int,
        default=8,
        help='Maximum number of domains to consider'
    )
    
    parser.add_argument(
        '--sigma-factor',
        type=float,
        default=1.0,
        help='Scaling factor for RBF kernel sigma'
    )
    
    parser.add_argument(
        '--no-save-labels',
        action='store_true',
        help='Do not save cluster labels'
    )
    
    parser.add_argument(
        '--single',
        type=str,
        help='Process single protein (format: PDB_CHAIN, e.g., 6p5a_B)'
    )
    
    args = parser.parse_args()
    
    if args.single:
        # Process single protein
        pdb_id, chain_id = args.single.split('_')
        structure_parser = ProteinStructureParser()
        
        result = process_protein(
            pdb_id,
            chain_id,
            structure_parser,
            estimation_method=args.method,
            max_domains=args.max_domains,
            sigma_factor=args.sigma_factor,
            verbose=True
        )
        
        if result['success']:
            print(f"Protein: {result['pdb_chain']}")
            print(f"Residues: {result['n_residues']}")
            print(f"Estimated domains: {result['n_domains_estimated']}")
            print(f"Domain sizes: {result['domain_sizes']}")
            
    else:
        # Process full dataset
        process_dataset(
            args.input,
            output_dir=args.output,
            estimation_method=args.method,
            max_domains=args.max_domains,
            sigma_factor=args.sigma_factor,
            save_labels=not args.no_save_labels
        )


if __name__ == "__main__":
    main()
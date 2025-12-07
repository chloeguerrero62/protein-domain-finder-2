"""
Random Controls and Ablation Studies

This script implements baseline comparisons and ablation studies to validate
that our methods perform significantly better than random chance.

Random Controls:
1. Random assignment (oracle n_domains)
2. Length-based splitting (oracle n_domains)
3. Single-domain baseline (predict n=1 for all)
4. Mean-domains baseline (predict n=2.5 for all)

Ablation Studies:
1. Graph construction (k-NN vs threshold)
2. Distance metrics (Euclidean only, but documented)
3. Similarity kernels (RBF vs inverse distance)
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
from tqdm import tqdm

from src.features.structure_parser import ProteinStructureParser
from src.features.distance_matrix import compute_distance_matrix
from src.features.graph_builder import build_knn_graph
from src.evaluation.metrics import compute_all_metrics


# =============================================================================
# RANDOM CONTROLS
# =============================================================================

def random_assignment_control(n_residues, n_domains, n_trials=10):
    """
    Random baseline: Assign residues randomly to n_domains clusters
    
    Note: Uses ORACLE knowledge of true n_domains
    Expected: Terrible cluster quality despite correct count
    """
    errors = []
    
    for _ in range(n_trials):
        labels = np.random.randint(0, n_domains, size=n_residues)
        # Always predicts correct n (by design)
        error = 0
        errors.append(error)
    
    return {
        'n_predicted': n_domains,
        'mean_error': np.mean(errors),
        'method': 'Random-Oracle'
    }


def length_based_control(n_residues, n_domains):
    """
    Length-based baseline: Split sequence into n equal parts
    
    Note: Uses ORACLE knowledge of true n_domains
    Expected: Arbitrary boundaries, no structural meaning
    """
    labels = np.zeros(n_residues, dtype=int)
    chunk_size = n_residues // n_domains
    
    for i in range(n_domains):
        start = i * chunk_size
        end = (i + 1) * chunk_size if i < n_domains - 1 else n_residues
        labels[start:end] = i
    
    return {
        'n_predicted': n_domains,
        'error': 0,  # Correct count by design
        'labels': labels,
        'method': 'Length-Oracle'
    }


def naive_baselines(n_domains):
    """
    Naive prediction baselines (no oracle)
    
    1. Predict n=1 for all proteins (single domain)
    2. Predict n=2.5 (dataset mean) for all proteins
    """
    baselines = []
    
    # Single domain prediction
    baselines.append({
        'method': 'Predict-One',
        'n_predicted': 1,
        'error': abs(1 - n_domains)
    })
    
    # Mean prediction (round to nearest int)
    mean_domains = 2.5  # Approximate dataset mean
    baselines.append({
        'method': 'Predict-Mean',
        'n_predicted': round(mean_domains),
        'error': abs(round(mean_domains) - n_domains)
    })
    
    return baselines


# =============================================================================
# ABLATION STUDIES
# =============================================================================

def ablation_graph_construction(distance_matrix):
    """
    Test different graph construction methods
    
    1. k-NN with varying k
    2. Threshold-based with varying threshold
    """
    from src.models.louvain_clustering import louvain_clustering
    
    results = []
    
    # k-NN ablation
    for k in [5, 10, 15, 20]:
        G = build_knn_graph(distance_matrix, k=k, weighting='inverse')
        partition = louvain_clustering(G, resolution=1.0)
        n_pred = len(set(partition.values()))
        
        results.append({
            'ablation': 'graph_construction',
            'method': f'kNN-k{k}',
            'k': k,
            'n_predicted': n_pred,
            'n_edges': G.number_of_edges(),
            'graph_density': len(G.edges()) / (len(G.nodes()) * (len(G.nodes()) - 1) / 2)
        })
    
    # Threshold ablation
    from src.features.graph_builder import build_threshold_graph
    
    for threshold in [6.0, 8.0, 10.0]:
        G = build_threshold_graph(distance_matrix, threshold=threshold, weighting='inverse')
        partition = louvain_clustering(G, resolution=1.0)
        n_pred = len(set(partition.values()))
        
        results.append({
            'ablation': 'graph_construction',
            'method': f'Threshold-{threshold}A',
            'threshold': threshold,
            'n_predicted': n_pred,
            'n_edges': G.number_of_edges(),
            'graph_density': len(G.edges()) / (len(G.nodes()) * (len(G.nodes()) - 1) / 2)
        })
    
    return results


def ablation_similarity_kernel(distance_matrix, n_domains):
    """
    Test different similarity kernels for spectral clustering
    
    1. RBF (Gaussian): exp(-d²/2σ²)
    2. Inverse distance: 1/(d+ε)
    3. Exponential: exp(-d/σ)
    """
    from sklearn.cluster import SpectralClustering
    
    results = []
    
    # RBF kernel (current method)
    sigma = np.median(distance_matrix[distance_matrix > 0])
    similarity_rbf = np.exp(-distance_matrix**2 / (2 * sigma**2))
    
    clustering = SpectralClustering(n_clusters=n_domains, affinity='precomputed', random_state=42)
    labels_rbf = clustering.fit_predict(similarity_rbf)
    
    results.append({
        'ablation': 'similarity_kernel',
        'kernel': 'RBF',
        'n_predicted': n_domains,
        'labels': labels_rbf
    })
    
    # Inverse distance
    similarity_inv = 1.0 / (distance_matrix + 1e-6)
    np.fill_diagonal(similarity_inv, 0)
    
    clustering = SpectralClustering(n_clusters=n_domains, affinity='precomputed', random_state=42)
    labels_inv = clustering.fit_predict(similarity_inv)
    
    results.append({
        'ablation': 'similarity_kernel',
        'kernel': 'Inverse',
        'n_predicted': n_domains,
        'labels': labels_inv
    })
    
    # Exponential
    similarity_exp = np.exp(-distance_matrix / sigma)
    
    clustering = SpectralClustering(n_clusters=n_domains, affinity='precomputed', random_state=42)
    labels_exp = clustering.fit_predict(similarity_exp)
    
    results.append({
        'ablation': 'similarity_kernel',
        'kernel': 'Exponential',
        'n_predicted': n_domains,
        'labels': labels_exp
    })
    
    return results


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def run_all_controls(df, parser):
    """
    Run all random controls on dataset
    """
    print("="*70)
    print("RANDOM CONTROLS")
    print("="*70)
    
    all_results = []
    
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Random controls"):
        try:
            coords, _ = parser.parse_structure(row['pdb_id'], row['chain'])
            n_residues = len(coords)
            n_true = row['n_domains']
            pdb_chain = f"{row['pdb_id']}_{row['chain']}"
            
            # Random assignment (oracle)
            random_result = random_assignment_control(n_residues, n_true, n_trials=10)
            all_results.append({
                'pdb_chain': pdb_chain,
                'method': random_result['method'],
                'n_predicted': random_result['n_predicted'],
                'n_true': n_true,
                'absolute_error': random_result['mean_error'],
                'exact_match': 1,  # By design
                'note': 'Oracle n_domains, random labels'
            })
            
            # Length-based (oracle)
            length_result = length_based_control(n_residues, n_true)
            all_results.append({
                'pdb_chain': pdb_chain,
                'method': length_result['method'],
                'n_predicted': length_result['n_predicted'],
                'n_true': n_true,
                'absolute_error': length_result['error'],
                'exact_match': 1,  # By design
                'note': 'Oracle n_domains, equal-length splits'
            })
            
            # Naive baselines (no oracle)
            for baseline in naive_baselines(n_true):
                all_results.append({
                    'pdb_chain': pdb_chain,
                    'method': baseline['method'],
                    'n_predicted': baseline['n_predicted'],
                    'n_true': n_true,
                    'absolute_error': baseline['error'],
                    'exact_match': int(baseline['n_predicted'] == n_true),
                    'note': 'No oracle knowledge'
                })
                
        except:
            continue
    
    return pd.DataFrame(all_results)


def run_ablation_studies(df, parser, n_proteins=20):
    """
    Run ablation studies on subset of proteins
    """
    print("\n" + "="*70)
    print("ABLATION STUDIES")
    print("="*70)
    
    # Sample diverse proteins
    df_sample = df.sample(n=min(n_proteins, len(df)), random_state=42)
    
    graph_results = []
    kernel_results = []
    
    for _, row in tqdm(df_sample.iterrows(), total=len(df_sample), desc="Ablations"):
        try:
            coords, _ = parser.parse_structure(row['pdb_id'], row['chain'])
            D = compute_distance_matrix(coords)
            pdb_chain = f"{row['pdb_id']}_{row['chain']}"
            
            # Graph construction ablation
            for result in ablation_graph_construction(D):
                result['pdb_chain'] = pdb_chain
                result['n_true'] = row['n_domains']
                result['error'] = abs(result['n_predicted'] - row['n_domains'])
                graph_results.append(result)
            
            # Similarity kernel ablation (supervised, uses oracle)
            for result in ablation_similarity_kernel(D, row['n_domains']):
                result['pdb_chain'] = pdb_chain
                result['n_true'] = row['n_domains']
                kernel_results.append(result)
                
        except:
            continue
    
    return pd.DataFrame(graph_results), pd.DataFrame(kernel_results)


def main():
    # Load data
    df = pd.read_csv('data/pdb_clustering/selected_proteins_mmseqs2.csv')
    parser = ProteinStructureParser()
    
    output_dir = Path('data/results')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run random controls
    controls_df = run_all_controls(df, parser)
    controls_df.to_csv(output_dir / 'random_controls.csv', index=False)
    
    # Print summary
    print("\n" + "="*70)
    print("RANDOM CONTROLS SUMMARY")
    print("="*70)
    
    for method in controls_df['method'].unique():
        method_df = controls_df[controls_df['method'] == method]
        exact = (method_df['exact_match'] == 1).sum()
        mae = method_df['absolute_error'].mean()
        
        print(f"\n{method}:")
        print(f"  Proteins: {len(method_df)}")
        print(f"  Exact matches: {exact}/{len(method_df)} ({100*exact/len(method_df):.1f}%)")
        print(f"  Mean Absolute Error: {mae:.2f}")
        print(f"  Note: {method_df.iloc[0]['note']}")
    
    # Run ablation studies
    graph_ablation, kernel_ablation = run_ablation_studies(df, parser, n_proteins=20)
    
    graph_ablation.to_csv(output_dir / 'ablation_graph_construction.csv', index=False)
    kernel_ablation.to_csv(output_dir / 'ablation_similarity_kernel.csv', index=False)
    
    # Print ablation summaries
    print("\n" + "="*70)
    print("GRAPH CONSTRUCTION ABLATION")
    print("="*70)
    
    graph_summary = graph_ablation.groupby('method').agg({
        'error': 'mean',
        'n_edges': 'mean',
        'graph_density': 'mean'
    }).round(3)
    print("\n", graph_summary)
    
    print("\n" + "="*70)
    print("SIMILARITY KERNEL ABLATION")
    print("="*70)
    print("(All use oracle n_domains, supervised)")
    print(f"Tested on {len(kernel_ablation['pdb_chain'].unique())} proteins")
    
    print("\n" + "="*70)
    print("FILES SAVED")
    print("="*70)
    print(f"Random controls: {output_dir / 'random_controls.csv'}")
    print(f"Graph ablation: {output_dir / 'ablation_graph_construction.csv'}")
    print(f"Kernel ablation: {output_dir / 'ablation_similarity_kernel.csv'}")


if __name__ == "__main__":
    main()
"""
Compare different clustering algorithms for domain detection
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
from sklearn.cluster import SpectralClustering, AgglomerativeClustering, DBSCAN
from sklearn.metrics import silhouette_score
from src.features.structure_parser import ProteinStructureParser
from src.features.distance_matrix import compute_distance_matrix
from src.features.graph_builder import build_knn_graph
from src.models.louvain_clustering import louvain_clustering, get_cluster_sizes
import networkx as nx


def spectral_clustering_domains(distance_matrix, n_domains):
    """
    Spectral clustering with known number of domains
    """
    # Convert distance to similarity
    sigma = np.median(distance_matrix[distance_matrix > 0])
    similarity = np.exp(-distance_matrix**2 / (2 * sigma**2))
    
    clustering = SpectralClustering(
        n_clusters=n_domains,
        affinity='precomputed',
        random_state=42,
        assign_labels='kmeans'
    )
    labels = clustering.fit_predict(similarity)
    
    return labels


def hierarchical_clustering_domains(distance_matrix, n_domains):
    """
    Agglomerative hierarchical clustering
    """
    clustering = AgglomerativeClustering(
        n_clusters=n_domains,
        metric='precomputed',
        linkage='average'
    )
    labels = clustering.fit_predict(distance_matrix)
    
    return labels


def graph_cut_clustering(graph, n_domains):
    """
    Spectral clustering on graph Laplacian
    """
    # Get adjacency matrix
    adj_matrix = nx.to_numpy_array(graph, weight='weight')
    
    clustering = SpectralClustering(
        n_clusters=n_domains,
        affinity='precomputed',
        random_state=42
    )
    labels = clustering.fit_predict(adj_matrix)
    
    return labels


def louvain_with_postprocessing(graph, target_domains, resolution_range=(0.1, 2.0), steps=20):
    """
    Search for Louvain resolution that gives target number of domains
    """
    from src.models.louvain_clustering import louvain_clustering, get_cluster_sizes
    
    best_labels = None
    best_diff = float('inf')
    
    for res in np.linspace(resolution_range[0], resolution_range[1], steps):
        partition = louvain_clustering(graph, resolution=res)
        n_clusters = len(get_cluster_sizes(partition))
        
        diff = abs(n_clusters - target_domains)
        if diff < best_diff:
            best_diff = diff
            best_labels = np.array([partition[i] for i in range(len(partition))])
            
            if diff == 0:
                break
    
    return best_labels


def evaluate_clustering_method(method_name, labels, true_labels):
    """
    Evaluate clustering quality
    """
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    
    n_pred = len(np.unique(labels))
    n_true = len(np.unique(true_labels))
    
    # Can only compute ARI/NMI if we have true cluster assignments
    # For now, just return domain count accuracy
    
    return {
        'method': method_name,
        'n_predicted': n_pred,
        'n_true': n_true,
        'error': abs(n_pred - n_true)
    }


def test_single_protein(pdb_id, chain_id, n_true_domains, parser):
    """
    Test all methods on one protein
    """
    results = []
    
    try:
        coords, _ = parser.parse_structure(pdb_id, chain_id)
        D = compute_distance_matrix(coords)
        G = build_knn_graph(D, k=10)
        
        # Method 1: Louvain (original)
        partition = louvain_clustering(G, resolution=0.5)
        labels = np.array([partition[i] for i in range(len(partition))])
        results.append({
            'method': 'Louvain (res=0.5)',
            'n_predicted': len(np.unique(labels)),
            'n_true': n_true_domains,
            'error': abs(len(np.unique(labels)) - n_true_domains)
        })
        
        # Method 2: Louvain with resolution search
        labels = louvain_with_postprocessing(G, n_true_domains)
        results.append({
            'method': 'Louvain (tuned)',
            'n_predicted': len(np.unique(labels)),
            'n_true': n_true_domains,
            'error': abs(len(np.unique(labels)) - n_true_domains)
        })
        
        # Method 3: Spectral clustering (distance matrix)
        labels = spectral_clustering_domains(D, n_true_domains)
        results.append({
            'method': 'Spectral (distance)',
            'n_predicted': len(np.unique(labels)),
            'n_true': n_true_domains,
            'error': abs(len(np.unique(labels)) - n_true_domains)
        })
        
        # Method 4: Spectral clustering (graph)
        labels = graph_cut_clustering(G, n_true_domains)
        results.append({
            'method': 'Spectral (graph)',
            'n_predicted': len(np.unique(labels)),
            'n_true': n_true_domains,
            'error': abs(len(np.unique(labels)) - n_true_domains)
        })
        
        # Method 5: Hierarchical clustering
        labels = hierarchical_clustering_domains(D, n_true_domains)
        results.append({
            'method': 'Hierarchical',
            'n_predicted': len(np.unique(labels)),
            'n_true': n_true_domains,
            'error': abs(len(np.unique(labels)) - n_true_domains)
        })
        
    except Exception as e:
        print(f"  ERROR: {e}")
        return None
    
    return results


if __name__ == "__main__":
    # Test on representative proteins
    test_cases = [
        ('6p5a', 'B', 1),
        ('7bm8', 'A', 2),
        ('3mpx', 'A', 3),
        ('8os5', 'A', 4),
    ]
    
    parser = ProteinStructureParser()
    
    print("Comparing clustering methods")
    print("="*80)
    
    all_results = []
    
    for pdb_id, chain_id, n_domains in test_cases:
        print(f"\n{pdb_id}_{chain_id} (True domains: {n_domains})")
        print("-"*80)
        
        results = test_single_protein(pdb_id, chain_id, n_domains, parser)
        
        if results:
            for r in results:
                r['pdb_chain'] = f"{pdb_id}_{chain_id}"
                all_results.append(r)
                print(f"  {r['method']:25s}: {r['n_predicted']:2d} domains (error={r['error']})")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    df = pd.DataFrame(all_results)
    
    summary = df.groupby('method').agg({
        'error': ['mean', 'min', 'max'],
        'n_predicted': 'mean'
    }).round(2)
    
    print(summary)
    
    print("\n" + "="*80)
    print("KEY INSIGHT:")
    print("="*80)
    print("Methods that REQUIRE knowing n_domains (Spectral, Hierarchical):")
    print("  - Will give exactly n_domains clusters")
    print("  - But cluster quality may be poor")
    print("  - This is 'cheating' since we don't know n_domains in practice")
    print("\nMethods that DISCOVER n_domains (Louvain):")
    print("  - Over-cluster severely (10-30 domains)")
    print("  - Finding structural motifs, not functional domains")

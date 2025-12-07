"""
Compare multiple clustering algorithms for protein domain detection
"""

import numpy as np
from sklearn.cluster import SpectralClustering, AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import networkx as nx
from src.models.louvain_clustering import louvain_clustering, get_cluster_sizes
from src.models.domain_count_estimator import estimate_domain_count_silhouette


def apply_louvain(distance_matrix, graph, n_domains=None):
    """Louvain clustering (unsupervised)"""
    partition = louvain_clustering(graph, resolution=0.5)
    labels = np.array([partition[i] for i in range(len(partition))])
    return labels


def apply_spectral_distance(distance_matrix, n_domains):
    """Spectral clustering on distance matrix (supervised)"""
    sigma = np.median(distance_matrix[distance_matrix > 0])
    similarity = np.exp(-distance_matrix**2 / (2 * sigma**2))
    
    clustering = SpectralClustering(
        n_clusters=n_domains,
        affinity='precomputed',
        random_state=42
    )
    return clustering.fit_predict(similarity)


def apply_spectral_graph(graph, n_domains):
    """Spectral clustering on graph (supervised)"""
    adj_matrix = nx.to_numpy_array(graph, weight='weight')
    
    clustering = SpectralClustering(
        n_clusters=n_domains,
        affinity='precomputed',
        random_state=42
    )
    return clustering.fit_predict(adj_matrix)


def apply_hierarchical(distance_matrix, n_domains):
    """Hierarchical clustering (supervised)"""
    clustering = AgglomerativeClustering(
        n_clusters=n_domains,
        metric='precomputed',
        linkage='average'
    )
    return clustering.fit_predict(distance_matrix)


def apply_two_stage_spectral(distance_matrix, graph):
    """Two-stage: estimate n_domains, then spectral clustering"""
    n_estimated, _ = estimate_domain_count_silhouette(distance_matrix, max_domains=8)
    return apply_spectral_distance(distance_matrix, n_estimated), n_estimated


def evaluate_clustering(labels_pred, labels_true):
    """Compute evaluation metrics"""
    if labels_true is None:
        return {
            'n_clusters': len(np.unique(labels_pred)),
            'ari': None,
            'nmi': None
        }
    
    return {
        'n_clusters': len(np.unique(labels_pred)),
        'ari': adjusted_rand_score(labels_true, labels_pred),
        'nmi': normalized_mutual_info_score(labels_true, labels_pred)
    }

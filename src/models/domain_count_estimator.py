"""
Estimate number of protein domains from distance matrix
"""

import numpy as np
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score


def estimate_domain_count_silhouette(distance_matrix, max_domains=10):
    """
    Estimate number of domains using silhouette score
    
    Try clustering with k=2 to max_domains, select k with best silhouette score
    """
    # Convert distance to similarity
    sigma = np.median(distance_matrix[distance_matrix > 0])
    similarity = np.exp(-distance_matrix**2 / (2 * sigma**2))
    
    best_k = 2
    best_score = -1
    
    scores = []
    
    for k in range(2, min(max_domains + 1, len(distance_matrix) // 10)):
        clustering = SpectralClustering(
            n_clusters=k,
            affinity='precomputed',
            random_state=42,
            assign_labels='kmeans'
        )
        labels = clustering.fit_predict(similarity)
        
        # Silhouette score (higher is better)
        score = silhouette_score(distance_matrix, labels, metric='precomputed')
        scores.append((k, score))
        
        if score > best_score:
            best_score = score
            best_k = k
    
    return best_k, scores


def estimate_domain_count_eigengap(distance_matrix, max_domains=10):
    """
    Estimate using eigengap heuristic
    
    Look for largest gap in eigenvalues of graph Laplacian
    """
    from scipy.linalg import eigh
    
    # Convert to similarity
    sigma = np.median(distance_matrix[distance_matrix > 0])
    similarity = np.exp(-distance_matrix**2 / (2 * sigma**2))
    
    # Compute normalized Laplacian
    D = np.diag(np.sum(similarity, axis=1))
    D_inv_sqrt = np.diag(1.0 / np.sqrt(np.diag(D) + 1e-10))
    L = np.eye(len(similarity)) - D_inv_sqrt @ similarity @ D_inv_sqrt
    
    # Get eigenvalues
    eigenvalues, _ = eigh(L)
    eigenvalues = np.sort(eigenvalues)[:max_domains]
    
    # Find largest gap
    gaps = np.diff(eigenvalues)
    best_k = np.argmax(gaps) + 1
    
    return max(2, min(best_k, max_domains))


def estimate_domain_count_elbow(distance_matrix, max_domains=10):
    """
    Estimate using elbow method on inertia/variance
    """
    n = len(distance_matrix)
    
    # Convert distance matrix to similarity
    sigma = np.median(distance_matrix[distance_matrix > 0])
    similarity = np.exp(-distance_matrix**2 / (2 * sigma**2))
    
    inertias = []
    
    for k in range(2, min(max_domains + 1, n // 10)):
        clustering = SpectralClustering(
            n_clusters=k,
            affinity='precomputed',
            random_state=42
        )
        labels = clustering.fit_predict(similarity)
        
        # Compute within-cluster variance
        variance = 0
        for cluster_id in range(k):
            mask = labels == cluster_id
            if mask.sum() > 1:
                cluster_dists = distance_matrix[mask][:, mask]
                variance += cluster_dists.sum()
        
        inertias.append((k, variance))
    
    # Find elbow (simplified: largest second derivative)
    if len(inertias) < 3:
        return 2
    
    variances = np.array([v for _, v in inertias])
    second_diff = np.diff(variances, n=2)
    elbow_idx = np.argmax(second_diff) + 2
    
    return inertias[elbow_idx][0]


def estimate_domain_count_consensus(distance_matrix, max_domains=10):
    """
    Consensus estimate using multiple methods
    """
    estimates = []
    
    # Method 1: Silhouette
    k_sil, _ = estimate_domain_count_silhouette(distance_matrix, max_domains)
    estimates.append(k_sil)
    
    # Method 2: Eigengap
    k_eig = estimate_domain_count_eigengap(distance_matrix, max_domains)
    estimates.append(k_eig)
    
    # Method 3: Elbow
    k_elbow = estimate_domain_count_elbow(distance_matrix, max_domains)
    estimates.append(k_elbow)
    
    # Return median
    return int(np.median(estimates))

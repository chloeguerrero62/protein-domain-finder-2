"""
Build graphs with sequential constraints for protein domain detection
"""

import numpy as np
import networkx as nx


def build_sequential_knn_graph(distance_matrix, k=10, seq_weight=2.0, weighting='inverse'):
    """
    Build k-nearest neighbors graph with sequential neighbor bias
    
    This addresses over-clustering by ensuring sequential residues
    are preferentially connected, preventing artificial fragmentation.
    
    Parameters:
    -----------
    distance_matrix : np.ndarray
        Pairwise distance matrix (n × n)
    k : int
        Number of nearest neighbors (default: 10)
    seq_weight : float
        Factor to reduce distance to sequential neighbors (default: 2.0)
        Higher values = stronger bias toward sequence neighbors
    weighting : str
        Edge weighting scheme:
        - 'inverse': weight = 1 / distance
        - 'exponential': weight = exp(-distance / 10)
        - 'binary': weight = 1.0
    
    Returns:
    --------
    G : networkx.Graph
        Undirected graph with weighted edges
    """
    n = len(distance_matrix)
    
    G = nx.Graph()
    G.add_nodes_from(range(n))
    
    for i in range(n):
        distances = distance_matrix[i].copy()
        distances[i] = np.inf
        
        # Reduce distance to sequential neighbors
        if i > 0:
            distances[i-1] /= seq_weight
        if i < n-1:
            distances[i+1] /= seq_weight
        
        # Find k nearest neighbors (with biased distances)
        nearest_indices = np.argsort(distances)[:k]
        
        for j in nearest_indices:
            # Use ORIGINAL distance for weight calculation
            dist = distance_matrix[i, j]
            
            if weighting == 'inverse':
                weight = 1.0 / (dist + 1e-6)
            elif weighting == 'exponential':
                weight = np.exp(-dist / 10.0)
            elif weighting == 'binary':
                weight = 1.0
            else:
                raise ValueError(f"Unknown weighting scheme: {weighting}")
            
            if not G.has_edge(i, j):
                G.add_edge(i, j, weight=weight, distance=dist)
    
    return G


def build_contact_map_graph(distance_matrix, threshold=8.0, weighting='inverse'):
    """
    Build graph using contact distance threshold with sequential backbone
    
    Parameters:
    -----------
    distance_matrix : np.ndarray
        Pairwise distance matrix
    threshold : float
        Contact distance threshold in Angstroms (default: 8.0)
    weighting : str
        Edge weighting scheme
    
    Returns:
    --------
    G : networkx.Graph
    """
    n = len(distance_matrix)
    
    G = nx.Graph()
    G.add_nodes_from(range(n))
    
    # Always connect sequential neighbors with high weight
    for i in range(n-1):
        G.add_edge(i, i+1, weight=10.0, distance=0.0)
    
    # Add spatial contacts
    for i in range(n):
        for j in range(i+2, n):
            dist = distance_matrix[i, j]
            
            if dist < threshold:
                if weighting == 'inverse':
                    weight = 1.0 / (dist + 1e-6)
                elif weighting == 'exponential':
                    weight = np.exp(-dist / 10.0)
                elif weighting == 'binary':
                    weight = 1.0
                else:
                    raise ValueError(f"Unknown weighting: {weighting}")
                
                G.add_edge(i, j, weight=weight, distance=dist)
    
    return G

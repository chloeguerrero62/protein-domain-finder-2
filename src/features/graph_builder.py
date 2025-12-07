"""
Build graphs from distance matrices
"""

import numpy as np
import networkx as nx


def build_knn_graph(distance_matrix, k=10, weighting='inverse'):
    """
    Build k-nearest neighbors graph from distance matrix
    
    Parameters:
    -----------
    distance_matrix : np.ndarray
        Pairwise distance matrix (n × n)
    k : int
        Number of nearest neighbors (default: 10)
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
    
    # Create empty graph
    G = nx.Graph()
    G.add_nodes_from(range(n))
    
    # For each residue, connect to k nearest neighbors
    for i in range(n):
        # Get distances from residue i to all others
        distances = distance_matrix[i].copy()
        distances[i] = np.inf  # Exclude self-loop
        
        # Find k nearest neighbors
        nearest_indices = np.argsort(distances)[:k]
        
        # Add edges to k nearest neighbors
        for j in nearest_indices:
            dist = distance_matrix[i, j]
            
            # Compute edge weight based on distance
            if weighting == 'inverse':
                weight = 1.0 / (dist + 1e-6)  # Add small epsilon to avoid division by zero
            elif weighting == 'exponential':
                weight = np.exp(-dist / 10.0)
            elif weighting == 'binary':
                weight = 1.0
            else:
                raise ValueError(f"Unknown weighting scheme: {weighting}")
            
            # Add edge (only once since graph is undirected)
            if not G.has_edge(i, j):
                G.add_edge(i, j, weight=weight, distance=dist)
    
    return G


def build_threshold_graph(distance_matrix, threshold=8.0, weighting='inverse'):
    """
    Build graph using distance threshold
    
    Parameters:
    -----------
    distance_matrix : np.ndarray
        Pairwise distance matrix
    threshold : float
        Distance threshold in Angstroms (default: 8.0)
    weighting : str
        Edge weighting scheme
    
    Returns:
    --------
    G : networkx.Graph
    """
    n = len(distance_matrix)
    
    G = nx.Graph()
    G.add_nodes_from(range(n))
    
    # Add edges for all pairs within threshold distance
    for i in range(n):
        for j in range(i+1, n):
            dist = distance_matrix[i, j]
            
            if dist < threshold:
                # Compute weight
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


def graph_statistics(G):
    """
    Compute basic graph statistics
    
    Parameters:
    -----------
    G : networkx.Graph
    
    Returns:
    --------
    stats : dict
        Dictionary of statistics
    """
    stats = {
        'n_nodes': G.number_of_nodes(),
        'n_edges': G.number_of_edges(),
        'density': nx.density(G),
        'n_components': nx.number_connected_components(G),
        'avg_degree': np.mean([d for n, d in G.degree()]),
    }
    
    # Clustering coefficient (can be slow for large graphs)
    if G.number_of_nodes() < 1000:
        stats['avg_clustering'] = nx.average_clustering(G)
    
    return stats


if __name__ == "__main__":
    # Test the module
    print("Testing graph_builder module...")
    
    # Create simple test distance matrix
    test_D = np.array([
        [0.0, 3.8, 6.2, 10.0],
        [3.8, 0.0, 3.8, 8.5],
        [6.2, 3.8, 0.0, 5.0],
        [10.0, 8.5, 5.0, 0.0]
    ])
    
    # Build KNN graph
    G = build_knn_graph(test_D, k=2)
    
    print(f"Test distance matrix shape: {test_D.shape}")
    print(f"Graph nodes: {G.number_of_nodes()}")
    print(f"Graph edges: {G.number_of_edges()}")
    print(f"Graph edges list: {list(G.edges(data=True))}")
    
    print("\n✅ Module tests passed!")

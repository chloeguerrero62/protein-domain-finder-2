"""
Louvain community detection for protein domain identification
"""

import numpy as np
import networkx as nx
from community import community_louvain


def louvain_clustering(graph, resolution=1.0, random_state=42):
    """
    Apply Louvain community detection to identify protein domains
    
    Parameters:
    -----------
    graph : networkx.Graph
        Protein residue graph with weighted edges
    resolution : float
        Resolution parameter (higher = more clusters, default: 1.0)
    random_state : int
        Random seed for reproducibility
    
    Returns:
    --------
    partition : dict
        Dictionary mapping node ID to cluster ID
        Example: {0: 0, 1: 0, 2: 1, 3: 1, ...}
    
    Notes:
    ------
    Uses the python-louvain library (community package)
    """
    
    # Run Louvain algorithm
    partition = community_louvain.best_partition(
        graph,
        resolution=resolution,
        random_state=random_state
    )
    
    return partition


def partition_to_labels(partition):
    """
    Convert partition dict to numpy array of labels
    
    Parameters:
    -----------
    partition : dict
        Node ID -> cluster ID mapping
    
    Returns:
    --------
    labels : np.ndarray
        Array of cluster labels in node order
    """
    n_nodes = len(partition)
    labels = np.zeros(n_nodes, dtype=int)
    
    for node_id, cluster_id in partition.items():
        labels[node_id] = cluster_id
    
    return labels


def get_cluster_sizes(partition):
    """
    Get size of each cluster
    
    Parameters:
    -----------
    partition : dict
        Node ID -> cluster ID mapping
    
    Returns:
    --------
    sizes : dict
        Cluster ID -> size mapping
    """
    from collections import Counter
    return dict(Counter(partition.values()))


def get_modularity(graph, partition):
    """
    Compute modularity of partition
    
    Parameters:
    -----------
    graph : networkx.Graph
    partition : dict
        Node ID -> cluster ID mapping
    
    Returns:
    --------
    modularity : float
        Modularity score (-1 to 1, higher is better)
    """
    return community_louvain.modularity(partition, graph)


def filter_small_clusters(partition, min_size=10):
    """
    Merge clusters smaller than min_size into nearest larger cluster
    
    Parameters:
    -----------
    partition : dict
        Node ID -> cluster ID mapping
    min_size : int
        Minimum cluster size (default: 10)
    
    Returns:
    --------
    filtered_partition : dict
        Updated partition with small clusters merged
    """
    from collections import Counter
    
    cluster_sizes = Counter(partition.values())
    
    # If all clusters are large enough, return as-is
    if all(size >= min_size for size in cluster_sizes.values()):
        return partition.copy()
    
    # Otherwise, merge small clusters
    # (Simple approach: assign to cluster 0)
    filtered = {}
    for node, cluster in partition.items():
        if cluster_sizes[cluster] >= min_size:
            filtered[node] = cluster
        else:
            # Assign to largest cluster (simple heuristic)
            largest_cluster = max(cluster_sizes, key=cluster_sizes.get)
            filtered[node] = largest_cluster
    
    return filtered


def partition_to_ranges(partition):
    """
    Convert partition to domain ranges
    
    Parameters:
    -----------
    partition : dict
        Node ID -> cluster ID mapping
    
    Returns:
    --------
    domains : list of tuples
        List of (start, end, cluster_id) for each contiguous domain
    """
    # Convert to sorted array
    labels = partition_to_labels(partition)
    
    domains = []
    current_cluster = labels[0]
    start = 0
    
    for i in range(1, len(labels)):
        if labels[i] != current_cluster:
            # Domain boundary
            domains.append((start, i-1, current_cluster))
            start = i
            current_cluster = labels[i]
    
    # Add final domain
    domains.append((start, len(labels)-1, current_cluster))
    
    return domains


if __name__ == "__main__":
    # Test the module
    print("Testing Louvain clustering module...")
    
    # Create a simple test graph with 2 communities
    G = nx.Graph()
    
    # Community 1: nodes 0-4
    for i in range(5):
        for j in range(i+1, 5):
            G.add_edge(i, j, weight=1.0)
    
    # Community 2: nodes 5-9
    for i in range(5, 10):
        for j in range(i+1, 10):
            G.add_edge(i, j, weight=1.0)
    
    # Weak connection between communities
    G.add_edge(4, 5, weight=0.1)
    
    # Run clustering
    partition = louvain_clustering(G, resolution=1.0)
    
    print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"Partition: {partition}")
    print(f"Cluster sizes: {get_cluster_sizes(partition)}")
    print(f"Modularity: {get_modularity(G, partition):.3f}")

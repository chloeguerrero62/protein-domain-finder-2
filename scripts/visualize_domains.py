"""
Visualize domain assignments for a single protein
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import matplotlib.pyplot as plt
import pickle


def visualize_partition(pdb_chain, output_dir='data/visualizations'):
    """Visualize domain assignments"""
    
    # Load partition
    labels = np.load(f'data/clusters/{pdb_chain}_partition.npy')
    
    # Load graph for context
    with open(f'data/graphs/{pdb_chain}_graph.pkl', 'rb') as f:
        G = pickle.load(f)
    
    # Create figure
    fig, axes = plt.subplots(2, 1, figsize=(12, 6))
    
    # Plot 1: Domain assignments along sequence
    ax1 = axes[0]
    colors = plt.cm.Set3(labels / max(labels))
    ax1.scatter(range(len(labels)), labels, c=colors, s=20, alpha=0.7)
    ax1.set_xlabel('Residue Position')
    ax1.set_ylabel('Domain ID')
    ax1.set_title(f'{pdb_chain} - Domain Assignments')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Domain sizes
    ax2 = axes[1]
    domain_sizes = np.bincount(labels)
    ax2.bar(range(len(domain_sizes)), domain_sizes, color='steelblue')
    ax2.set_xlabel('Domain ID')
    ax2.set_ylabel('Number of Residues')
    ax2.set_title('Domain Sizes')
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    # Save
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path / f'{pdb_chain}_domains.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved visualization: {output_path / f'{pdb_chain}_domains.png'}")


if __name__ == "__main__":
    # Visualize first few proteins
    cluster_files = sorted(Path('data/clusters').glob('*_partition.npy'))[:5]
    
    for cluster_file in cluster_files:
        pdb_chain = cluster_file.stem.replace('_partition', '')
        visualize_partition(pdb_chain)
    
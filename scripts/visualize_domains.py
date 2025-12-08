"""
Visualize domain assignments for multiple clustering methods

This script creates comparison visualizations showing:
1. Louvain clustering results
2. Two-Stage Spectral clustering results
3. Ground truth (Pfam domains)
4. Side-by-side comparisons
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import pickle


def load_ground_truth(pdb_chain):
    """Load ground truth domain count from CSV"""
    csv_file = Path('data/pdb_clustering/selected_proteins_mmseqs2.csv')
    if not csv_file.exists():
        return None
    
    df = pd.read_csv(csv_file)
    if 'pdb_chain' not in df.columns:
        df['pdb_chain'] = df['pdb_id'] + '_' + df['chain']
    
    row = df[df['pdb_chain'] == pdb_chain]
    if len(row) == 0:
        return None
    
    return int(row.iloc[0]['n_domains'])


def load_clustering_results(pdb_chain):
    """
    Load clustering results from all methods
    
    Returns dict with keys: 'louvain', 'two_stage', 'n_true'
    """
    results = {'n_true': load_ground_truth(pdb_chain)}
    
    # Load Louvain results
    louvain_file = Path(f'data/clusters/{pdb_chain}_partition.npy')
    if louvain_file.exists():
        results['louvain'] = np.load(louvain_file)
    
    # Load Two-Stage Spectral results
    two_stage_file = Path(f'data/results/two_stage_labels/{pdb_chain}_labels.npy')
    if two_stage_file.exists():
        results['two_stage'] = np.load(two_stage_file)
    
    return results


def visualize_single_method(ax, labels, method_name, n_true=None):
    """Visualize a single clustering method"""
    
    n_predicted = len(np.unique(labels))
    n_residues = len(labels)
    
    # Use consistent colormap
    colors = plt.cm.tab20(labels / max(labels))
    
    # Scatter plot of domain assignments
    ax.scatter(range(n_residues), labels, c=colors, s=30, alpha=0.8, edgecolors='black', linewidth=0.5)
    
    # Title with accuracy info
    title = f'{method_name}\n'
    if n_true is not None:
        error = abs(n_predicted - n_true)
        match = "✓" if n_predicted == n_true else "✗"
        title += f'Predicted: {n_predicted}, True: {n_true}, Error: {error} {match}'
    else:
        title += f'Predicted: {n_predicted} domains'
    
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.set_xlabel('Residue Position', fontsize=9)
    ax.set_ylabel('Domain ID', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-5, n_residues + 5)


def visualize_comparison(pdb_chain, output_dir='data/visualizations'):
    """
    Create comprehensive comparison visualization for a protein
    """
    
    # Load all results
    results = load_clustering_results(pdb_chain)
    n_true = results.get('n_true')
    
    # Check what data we have
    has_louvain = 'louvain' in results
    has_two_stage = 'two_stage' in results
    
    if not has_louvain and not has_two_stage:
        print(f"No clustering results found for {pdb_chain}")
        return
    
    # Determine number of subplots needed
    n_plots = sum([has_louvain, has_two_stage])
    
    # Create figure
    fig = plt.figure(figsize=(14, 4 * n_plots))
    gs = GridSpec(n_plots + 1, 2, figure=fig, height_ratios=[4]*n_plots + [1])
    
    plot_idx = 0
    
    # Plot Louvain if available
    if has_louvain:
        ax = fig.add_subplot(gs[plot_idx, :])
        visualize_single_method(ax, results['louvain'], 'Louvain Clustering (Unsupervised)', n_true)
        plot_idx += 1
    
    # Plot Two-Stage Spectral if available
    if has_two_stage:
        ax = fig.add_subplot(gs[plot_idx, :])
        visualize_single_method(ax, results['two_stage'], 'Two-Stage Spectral (Unsupervised)', n_true)
        plot_idx += 1
    
    # Add domain size comparison bar chart
    ax_left = fig.add_subplot(gs[plot_idx, 0])
    ax_right = fig.add_subplot(gs[plot_idx, 1])
    
    if has_louvain:
        domain_sizes_louvain = np.bincount(results['louvain'])
        ax_left.bar(range(len(domain_sizes_louvain)), domain_sizes_louvain, 
                    color='coral', alpha=0.7, edgecolor='black')
        ax_left.set_xlabel('Domain ID', fontsize=9)
        ax_left.set_ylabel('Number of Residues', fontsize=9)
        ax_left.set_title('Louvain: Domain Sizes', fontsize=10, fontweight='bold')
        ax_left.grid(True, alpha=0.3, axis='y')
    
    if has_two_stage:
        domain_sizes_two_stage = np.bincount(results['two_stage'])
        ax_right.bar(range(len(domain_sizes_two_stage)), domain_sizes_two_stage, 
                     color='steelblue', alpha=0.7, edgecolor='black')
        ax_right.set_xlabel('Domain ID', fontsize=9)
        ax_right.set_ylabel('Number of Residues', fontsize=9)
        ax_right.set_title('Two-Stage: Domain Sizes', fontsize=10, fontweight='bold')
        ax_right.grid(True, alpha=0.3, axis='y')
    
    # Main title
    fig.suptitle(f'Protein {pdb_chain} - Domain Prediction Comparison', 
                 fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    # Save
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    output_file = output_path / f'{pdb_chain}_comparison.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_file}")


def visualize_method_only(pdb_chain, method='louvain', output_dir='data/visualizations'):
    """
    Visualize results from a single method
    
    Parameters:
    -----------
    method : str
        'louvain' or 'two_stage'
    """
    
    n_true = load_ground_truth(pdb_chain)
    
    # Load labels
    if method == 'louvain':
        label_file = Path(f'data/clusters/{pdb_chain}_partition.npy')
        method_name = 'Louvain Clustering'
    elif method == 'two_stage':
        label_file = Path(f'data/results/two_stage_labels/{pdb_chain}_labels.npy')
        method_name = 'Two-Stage Spectral'
    else:
        raise ValueError(f"Unknown method: {method}")
    
    if not label_file.exists():
        print(f"Results not found: {label_file}")
        return
    
    labels = np.load(label_file)
    
    # Create figure
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    
    # Plot 1: Domain assignments
    visualize_single_method(axes[0], labels, method_name, n_true)
    
    # Plot 2: Domain sizes
    domain_sizes = np.bincount(labels)
    colors = 'coral' if method == 'louvain' else 'steelblue'
    axes[1].bar(range(len(domain_sizes)), domain_sizes, color=colors, alpha=0.7, edgecolor='black')
    axes[1].set_xlabel('Domain ID', fontsize=10)
    axes[1].set_ylabel('Number of Residues', fontsize=10)
    axes[1].set_title('Domain Sizes', fontsize=11, fontweight='bold')
    axes[1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    # Save
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    output_file = output_path / f'{pdb_chain}_{method}.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved: {output_file}")


def create_summary_figure(pdb_chains, output_dir='data/visualizations'):
    """
    Create a summary figure showing multiple proteins side-by-side
    
    Shows both Louvain and Two-Stage predictions for each protein
    """
    
    n_proteins = len(pdb_chains)
    fig, axes = plt.subplots(n_proteins, 2, figsize=(14, 3 * n_proteins))
    
    if n_proteins == 1:
        axes = axes.reshape(1, -1)
    
    for i, pdb_chain in enumerate(pdb_chains):
        results = load_clustering_results(pdb_chain)
        n_true = results.get('n_true')
        
        # Louvain
        if 'louvain' in results:
            labels = results['louvain']
            n_pred = len(np.unique(labels))
            colors = plt.cm.tab20(labels / max(labels))
            axes[i, 0].scatter(range(len(labels)), labels, c=colors, s=20, alpha=0.7)
            title = f'{pdb_chain}\nLouvain: {n_pred}'
            if n_true:
                title += f' (true: {n_true})'
            axes[i, 0].set_title(title, fontsize=9)
            axes[i, 0].set_ylabel('Domain ID', fontsize=8)
            axes[i, 0].grid(True, alpha=0.3)
        
        # Two-Stage
        if 'two_stage' in results:
            labels = results['two_stage']
            n_pred = len(np.unique(labels))
            colors = plt.cm.tab20(labels / max(labels))
            axes[i, 1].scatter(range(len(labels)), labels, c=colors, s=20, alpha=0.7)
            title = f'{pdb_chain}\nTwo-Stage: {n_pred}'
            if n_true:
                title += f' (true: {n_true})'
            axes[i, 1].set_title(title, fontsize=9)
            axes[i, 1].set_ylabel('Domain ID', fontsize=8)
            axes[i, 1].grid(True, alpha=0.3)
        
        # X-labels only on bottom row
        if i == n_proteins - 1:
            axes[i, 0].set_xlabel('Residue Position', fontsize=8)
            axes[i, 1].set_xlabel('Residue Position', fontsize=8)
    
    plt.tight_layout()
    
    # Save
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    output_file = output_path / 'summary_comparison.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved summary: {output_file}")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Visualize protein domain predictions'
    )
    
    parser.add_argument(
        '--protein',
        type=str,
        help='Single protein to visualize (e.g., 6p5a_B)'
    )
    
    parser.add_argument(
        '--method',
        type=str,
        choices=['louvain', 'two_stage', 'comparison', 'all'],
        default='comparison',
        help='Which method(s) to visualize'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        default='data/visualizations',
        help='Output directory'
    )
    
    parser.add_argument(
        '--batch',
        type=int,
        default=5,
        help='Number of proteins to visualize in batch mode'
    )
    
    parser.add_argument(
        '--summary',
        action='store_true',
        help='Create summary figure with multiple proteins'
    )
    
    args = parser.parse_args()
    
    if args.protein:
        # Visualize single protein
        if args.method == 'comparison' or args.method == 'all':
            visualize_comparison(args.protein, args.output)
        if args.method == 'louvain' or args.method == 'all':
            visualize_method_only(args.protein, 'louvain', args.output)
        if args.method == 'two_stage' or args.method == 'all':
            visualize_method_only(args.protein, 'two_stage', args.output)
    
    elif args.summary:
        # Create summary figure
        # Get first N proteins that have results
        results_df = pd.read_csv('data/results/two_stage_silhouette.csv')
        successful = results_df[results_df['success'] == True]
        pdb_chains = successful['pdb_chain'].head(args.batch).tolist()
        
        create_summary_figure(pdb_chains, args.output)
    
    else:
        # Batch mode - visualize first N proteins
        print(f"Batch mode: visualizing {args.batch} proteins")
        
        # Find proteins with clustering results
        louvain_files = sorted(Path('data/clusters').glob('*_partition.npy'))
        
        count = 0
        for cluster_file in louvain_files:
            if count >= args.batch:
                break
            
            pdb_chain = cluster_file.stem.replace('_partition', '')
            
            if args.method == 'comparison' or args.method == 'all':
                visualize_comparison(pdb_chain, args.output)
            if args.method == 'louvain' or args.method == 'all':
                visualize_method_only(pdb_chain, 'louvain', args.output)
            if args.method == 'two_stage' or args.method == 'all':
                visualize_method_only(pdb_chain, 'two_stage', args.output)
            
            count += 1
        
        print(f"\nVisualized {count} proteins")
        print(f"Saved to: {args.output}/")


if __name__ == "__main__":
    main()
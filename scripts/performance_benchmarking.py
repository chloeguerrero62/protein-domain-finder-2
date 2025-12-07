"""
Performance Benchmarking and Computational Requirements

This script measures runtime, memory usage, and scalability for all methods.

Metrics tracked:
1. Runtime per protein (by size)
2. Peak memory usage
3. Scalability with protein size
4. Disk space requirements
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import time
import tracemalloc
from tqdm import tqdm

from src.features.structure_parser import ProteinStructureParser
from src.features.distance_matrix import compute_distance_matrix
from src.features.graph_builder import build_knn_graph
from src.models.louvain_clustering import louvain_clustering
from src.evaluation.clustering_comparison import (
    apply_spectral_distance,
    apply_hierarchical,
    apply_two_stage_spectral
)


def measure_runtime_memory(func, *args, **kwargs):
    """
    Measure runtime and peak memory usage of a function
    """
    # Start memory tracking
    tracemalloc.start()
    
    # Measure runtime
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    
    # Get peak memory
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    return {
        'runtime_seconds': end_time - start_time,
        'peak_memory_mb': peak / 1024 / 1024,
        'result': result
    }


def benchmark_single_protein(pdb_id, chain_id, parser):
    """
    Benchmark all processing steps for one protein
    """
    results = {}
    
    # 1. Structure parsing
    perf = measure_runtime_memory(parser.parse_structure, pdb_id, chain_id)
    coords = perf['result'][0]
    results['parsing'] = {
        'runtime': perf['runtime_seconds'],
        'memory': perf['peak_memory_mb']
    }
    
    # 2. Distance matrix computation
    perf = measure_runtime_memory(compute_distance_matrix, coords)
    D = perf['result']
    results['distance_matrix'] = {
        'runtime': perf['runtime_seconds'],
        'memory': perf['peak_memory_mb'],
        'matrix_size_mb': D.nbytes / 1024 / 1024
    }
    
    # 3. Graph construction
    perf = measure_runtime_memory(build_knn_graph, D, k=10)
    G = perf['result']
    results['graph_construction'] = {
        'runtime': perf['runtime_seconds'],
        'memory': perf['peak_memory_mb'],
        'n_edges': G.number_of_edges()
    }
    
    # 4. Louvain clustering
    perf = measure_runtime_memory(louvain_clustering, G, resolution=1.0)
    results['louvain'] = {
        'runtime': perf['runtime_seconds'],
        'memory': perf['peak_memory_mb']
    }
    
    # 5. Spectral clustering (supervised, use n=3 as example)
    perf = measure_runtime_memory(apply_spectral_distance, D, 3)
    results['spectral'] = {
        'runtime': perf['runtime_seconds'],
        'memory': perf['peak_memory_mb']
    }
    
    # 6. Hierarchical clustering
    perf = measure_runtime_memory(apply_hierarchical, D, 3)
    results['hierarchical'] = {
        'runtime': perf['runtime_seconds'],
        'memory': perf['peak_memory_mb']
    }
    
    # 7. Two-stage spectral
    perf = measure_runtime_memory(apply_two_stage_spectral, D, G)
    results['two_stage'] = {
        'runtime': perf['runtime_seconds'],
        'memory': perf['peak_memory_mb']
    }
    
    results['n_residues'] = len(coords)
    
    return results


def benchmark_scalability(df, parser, size_bins):
    """
    Test scalability across protein sizes
    
    size_bins: [(min_size, max_size, label), ...]
    """
    print("="*70)
    print("SCALABILITY ANALYSIS")
    print("="*70)
    
    all_benchmarks = []
    
    for min_size, max_size, label in size_bins:
        # Get proteins in this size range
        size_df = df[(df['length'] >= min_size) & (df['length'] < max_size)]
        
        if len(size_df) == 0:
            continue
        
        # Sample up to 5 proteins
        sample_df = size_df.sample(n=min(5, len(size_df)), random_state=42)
        
        print(f"\n{label} ({min_size}-{max_size} residues): {len(sample_df)} proteins")
        
        for _, row in sample_df.iterrows():
            try:
                benchmark = benchmark_single_protein(
                    row['pdb_id'], 
                    row['chain'], 
                    parser
                )
                benchmark['size_bin'] = label
                benchmark['pdb_chain'] = f"{row['pdb_id']}_{row['chain']}"
                all_benchmarks.append(benchmark)
                
                print(f"  {row['pdb_id']}_{row['chain']} ({benchmark['n_residues']} residues): "
                      f"Total {sum(b['runtime'] for b in benchmark.values() if isinstance(b, dict)):.2f}s")
                
            except Exception as e:
                print(f"  {row['pdb_id']}_{row['chain']}: FAILED - {str(e)}")
                continue
    
    return all_benchmarks


def estimate_disk_usage():
    """
    Estimate disk space requirements
    """
    print("\n" + "="*70)
    print("DISK SPACE REQUIREMENTS")
    print("="*70)
    
    data_dir = Path('data')
    
    sizes = {}
    
    for subdir in data_dir.rglob('*'):
        if subdir.is_file():
            category = subdir.parent.name
            size_mb = subdir.stat().st_size / 1024 / 1024
            
            if category not in sizes:
                sizes[category] = 0
            sizes[category] += size_mb
    
    print("\nDisk usage by category:")
    total = 0
    for category, size in sorted(sizes.items(), key=lambda x: -x[1]):
        print(f"  {category:30s}: {size:8.1f} MB")
        total += size
    
    print(f"  {'TOTAL':30s}: {total:8.1f} MB ({total/1024:.2f} GB)")
    
    return sizes


def main():
    # Load data
    df = pd.read_csv('data/pdb_clustering/selected_proteins_mmseqs2.csv')
    parser = ProteinStructureParser()
    
    print("="*70)
    print("PERFORMANCE BENCHMARKING")
    print("="*70)
    print(f"Dataset: {len(df)} proteins")
    print(f"Size range: {df['length'].min()}-{df['length'].max()} residues\n")
    
    # Define size bins
    size_bins = [
        (50, 200, 'Small'),
        (200, 500, 'Medium'),
        (500, 1000, 'Large'),
        (1000, 2500, 'Very Large')
    ]
    
    # Run scalability benchmarks
    benchmarks = benchmark_scalability(df, parser, size_bins)
    
    # Convert to DataFrame for analysis
    benchmark_data = []
    for bench in benchmarks:
        base = {
            'pdb_chain': bench['pdb_chain'],
            'n_residues': bench['n_residues'],
            'size_bin': bench['size_bin']
        }
        
        for step, metrics in bench.items():
            if isinstance(metrics, dict) and 'runtime' in metrics:
                row = base.copy()
                row['step'] = step
                row['runtime_seconds'] = metrics['runtime']
                row['memory_mb'] = metrics['memory']
                benchmark_data.append(row)
    
    benchmark_df = pd.DataFrame(benchmark_data)
    
    # Save detailed results
    output_dir = Path('data/results')
    output_dir.mkdir(parents=True, exist_ok=True)
    benchmark_df.to_csv(output_dir / 'performance_benchmarks.csv', index=False)
    
    # Print summary by size bin
    print("\n" + "="*70)
    print("AVERAGE RUNTIME BY SIZE BIN")
    print("="*70)
    
    runtime_summary = benchmark_df.groupby(['size_bin', 'step'])['runtime_seconds'].mean().unstack(fill_value=0)
    print("\n", runtime_summary.round(3))
    
    # Print summary by method
    print("\n" + "="*70)
    print("AVERAGE RUNTIME BY METHOD")
    print("="*70)
    
    method_summary = benchmark_df.groupby('step').agg({
        'runtime_seconds': ['mean', 'std', 'min', 'max'],
        'memory_mb': ['mean', 'max']
    }).round(3)
    print("\n", method_summary)
    
    # Disk usage
    disk_usage = estimate_disk_usage()
    
    # Create summary report
    summary = {
        'dataset_size': len(df),
        'size_range': (int(df['length'].min()), int(df['length'].max())),
        'total_disk_gb': sum(disk_usage.values()) / 1024,
        'avg_runtime_per_protein': {
            'parsing': benchmark_df[benchmark_df['step'] == 'parsing']['runtime_seconds'].mean(),
            'distance_matrix': benchmark_df[benchmark_df['step'] == 'distance_matrix']['runtime_seconds'].mean(),
            'graph_construction': benchmark_df[benchmark_df['step'] == 'graph_construction']['runtime_seconds'].mean(),
            'louvain': benchmark_df[benchmark_df['step'] == 'louvain']['runtime_seconds'].mean(),
            'spectral': benchmark_df[benchmark_df['step'] == 'spectral']['runtime_seconds'].mean(),
            'two_stage': benchmark_df[benchmark_df['step'] == 'two_stage']['runtime_seconds'].mean()
        },
        'peak_memory_mb': benchmark_df['memory_mb'].max()
    }
    
    import json
    with open(output_dir / 'performance_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print("\n" + "="*70)
    print("ESTIMATED TOTAL RUNTIME (full dataset)")
    print("="*70)
    
    avg_total = benchmark_df.groupby('pdb_chain')['runtime_seconds'].sum().mean()
    estimated_full = avg_total * len(df)
    
    print(f"Average per protein: {avg_total:.2f} seconds")
    print(f"Estimated for {len(df)} proteins: {estimated_full:.1f} seconds ({estimated_full/60:.1f} minutes)")
    
    print("\n" + "="*70)
    print("FILES SAVED")
    print("="*70)
    print(f"Detailed benchmarks: {output_dir / 'performance_benchmarks.csv'}")
    print(f"Summary: {output_dir / 'performance_summary.json'}")


if __name__ == "__main__":
    main()
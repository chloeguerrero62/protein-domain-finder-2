## Overview

This project implements a **comprehensive comparative study** of clustering algorithms for automated protein domain identification from 3D structure:

### Methods Compared:
1. **Louvain Community Detection** (unsupervised) - Graph-based clustering
2. **Spectral Clustering** (supervised) - On distance matrices and graphs
3. **Hierarchical Clustering** (supervised) - Agglomerative approach
4. **Two-Stage Spectral** (unsupervised) - Automatic domain count estimation + spectral clustering

### Key Research Questions:
- Can unsupervised methods accurately predict domain count?
- How do graph-based vs. distance-based methods compare?
- What is the performance gap between supervised and unsupervised approaches?

## Quick Start
```bash
# Compare all methods
python scripts/compare_all_methods.py

# Results saved to:
# - data/results/method_comparison.csv (detailed)
# - data/results/method_summary.csv (summary table)
```

## Methods

### Clustering Algorithms

**1. Louvain (Unsupervised)**
- Graph community detection via modularity optimization
- No prior knowledge of domain count required
- Tends to over-cluster (finds structural motifs)

**2. Spectral Clustering (Supervised)**
- Eigendecomposition of similarity/adjacency matrices
- Requires true domain count (oracle knowledge)
- Two variants: distance-based and graph-based

**3. Hierarchical Clustering (Supervised)**
- Agglomerative clustering with average linkage
- Requires true domain count
- Works directly on distance matrices

**4. Two-Stage Spectral (Unsupervised)**
- Stage 1: Estimate domain count via silhouette score
- Stage 2: Apply spectral clustering with estimated count
- Best of both worlds approach

### Evaluation Metrics

- **Domain Count Accuracy**: Exact matches, Mean Absolute Error
- **Cluster Quality**: Adjusted Rand Index (ARI), Normalized Mutual Information (NMI)
- **Boundary Detection**: F1 score for domain boundaries

## Results

[This section will be populated after running experiments]

Expected findings:
- Supervised methods achieve near-perfect domain count accuracy
- Unsupervised Louvain over-clusters by 8-15 domains on average
- Two-stage spectral achieves moderate performance (MAE ~2-3 domains)
```

## **Project Structure (Updated)**
```
structure-domain-finder/
│
├── src/
│   ├── features/
│   │   ├── structure_parser.py
│   │   ├── distance_matrix.py
│   │   ├── graph_builder.py
│   │   └── sequential_graph_builder.py
│   ├── models/
│   │   ├── louvain_clustering.py
│   │   └── domain_count_estimator.py      # NEW
│   └── evaluation/                         # NEW
│       ├── clustering_comparison.py
│       └── metrics.py
│
├── scripts/
│   ├── compare_all_methods.py             # NEW - Main comparison
│   ├── compare_clustering_methods.py       # Keep for quick tests
│   ├── two_stage_pipeline.py              # Keep
│   └── ...
│
└── data/results/
    ├── method_comparison.csv               # NEW - Detailed results
    └── method_summary.csv                  # NEW - Summary table
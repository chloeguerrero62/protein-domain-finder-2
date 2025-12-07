# Project Organization - Structure-Domain Finder

## Directory Structure

```
structure-domain-finder/
├── data/                           # Data directory (mostly .gitignored)
│   ├── pdb_clustering/            # Dataset preparation
│   │   ├── pdb_seqres.txt         # Raw PDB sequences (2 GB, not in git)
│   │   ├── pdb_seqres_filtered.fasta
│   │   ├── pdb_representatives_30.fasta
│   │   ├── selected_proteins_mmseqs2.csv  # Final dataset (in git)
│   │   └── pdb_metadata.csv
│   ├── selected_structures/       # PDB files (~50 MB, not in git)
│   ├── distance_matrices/         # Distance matrices (~50 MB, not in git)
│   ├── graphs/                    # Protein graphs (~50 MB, not in git)
│   ├── clusters/                  # Clustering results (~10 MB, not in git)
│   └── results/                   # All evaluation results
│       ├── best_parameters.json
│       ├── hyperparameter_grid_search.csv
│       ├── test_set_performance.csv
│       ├── random_controls.csv
│       ├── ablation_graph_construction.csv
│       ├── ablation_similarity_kernel.csv
│       ├── performance_benchmarks.csv
│       ├── performance_summary.json
│       ├── method_comparison.csv
│       └── method_summary.csv
│
├── docs/                          # Documentation
│   ├── EVALUATION.md             # Complete evaluation framework
│   └── optimization_plan.md      # Hyperparameter optimization strategy
│
├── scripts/                       # Executable scripts
│   ├── prefilter.py              # Data preparation (step 1)
│   ├── cluster_protein_selection.py  # Stratified sampling (step 2)
│   ├── download_selected_pdbs.py     # Download structures (step 3)
│   ├── distance_matrix.py        # Feature extraction (step 4)
│   ├── graph_builder.py          # Graph construction (step 5)
│   ├── louvain_clustering.py     # Louvain method
│   ├── two_stage_pipeline.py     # Two-stage spectral
│   ├── compare_clustering_methods.py  # Quick method comparison
│   ├── compare_all_methods.py    # Full comparison
│   ├── visualize_domains.py      # Visualization
│   ├── hyperparameter_optimization.py  # Grid search (OPTIMIZATION)
│   ├── random_controls.py        # Baselines and ablations (CONTROLS)
│   └── performance_benchmarking.py    # Runtime/memory profiling (BENCHMARKS)
│
├── src/                           # Source modules
│   ├── features/
│   │   ├── structure_parser.py   # PDB parsing
│   │   ├── distance_matrix.py    # Distance computation
│   │   ├── graph_builder.py      # Graph construction
│   │   └── sequential_graph_builder.py  # Sequential bias (deprecated)
│   ├── models/
│   │   ├── louvain_clustering.py      # Louvain implementation
│   │   └── domain_count_estimator.py  # k estimation methods
│   └── evaluation/
│       ├── clustering_comparison.py   # Method wrappers
│       └── metrics.py                 # Evaluation metrics
│
├── notebooks/                     # Jupyter notebooks (TO CREATE)
│   └── full_workflow.ipynb       # Complete analysis notebook
│
├── .gitignore                     # Git ignore rules
├── README.md                      # Project overview and usage
├── environment.yml                # Conda environment
└── requirements.txt               # Pip requirements
```

## Workflow Execution Order

### Phase 1: Data Preparation (One-Time Setup)
```bash
# 1. Download PDB sequences
cd data/pdb_clustering
wget https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
gunzip pdb_seqres.txt.gz

# 2. Filter sequences
cd ../..
python scripts/prefilter.py

# 3. MMseqs2 clustering (in data/pdb_clustering)
mmseqs createdb pdb_seqres_filtered.fasta pdb_filtered_DB
mmseqs cluster pdb_filtered_DB pdb_clusters_30 tmp \
    --min-seq-id 0.3 -c 0.8 --cov-mode 1 --cluster-mode 2 -s 7.5
mmseqs createsubdb pdb_clusters_30 pdb_filtered_DB pdb_representatives_DB
mmseqs convert2fasta pdb_representatives_DB pdb_representatives_30.fasta

# 4. Select diverse protein set
python scripts/cluster_protein_selection.py

# 5. Download PDB structures
python scripts/download_selected_pdbs.py
```

### Phase 2: Feature Extraction
```bash
# 6. Compute distance matrices
python scripts/distance_matrix.py

# 7. Build graphs
python scripts/graph_builder.py
```

### Phase 3: Evaluation (Main Experiments)
```bash
# 8. Hyperparameter optimization (Two-Stage Spectral)
python scripts/hyperparameter_optimization.py
# Output: data/results/best_parameters.json
# Runtime: ~20 minutes

# 9. Random controls and ablation studies
python scripts/random_controls.py
# Output: data/results/random_controls.csv
#         data/results/ablation_*.csv
# Runtime: ~10 minutes

# 10. Performance benchmarking
python scripts/performance_benchmarking.py
# Output: data/results/performance_benchmarks.csv
# Runtime: ~5 minutes

# 11. Full method comparison (with optimized parameters)
python scripts/compare_all_methods.py
# Output: data/results/method_comparison.csv
# Runtime: ~15 minutes
```

### Phase 4: Analysis and Reporting
```bash
# 12. Generate visualizations
python scripts/visualize_domains.py

# 13. Create Jupyter notebook (interactive analysis)
jupyter notebook notebooks/full_workflow.ipynb
```

## Key Files by Purpose

### Assignment Requirements Checklist

| Requirement | Files | Notes |
|-------------|-------|-------|
| **Task** | README.md, docs/EVALUATION.md | Protein domain detection |
| **Data** | data/pdb_clustering/selected_proteins_mmseqs2.csv | 144 proteins, stratified |
| **Features** | src/features/*.py | Distance matrices, graphs |
| **Model Families** | src/models/*.py, src/evaluation/*.py | Louvain, Spectral, Hierarchical |
| **Hyperparameter Optimization** | scripts/hyperparameter_optimization.py | Two-Stage Spectral grid search |
| **Evaluation** | docs/EVALUATION.md, scripts/random_controls.py | Comprehensive framework |
| **Random Controls** | scripts/random_controls.py | 4 baselines + ablations |
| **Performance Estimation** | scripts/performance_benchmarking.py | Runtime, memory, scalability |
| **45-min Presentation** | TO CREATE | Slides from results |
| **Jupyter Notebook** | TO CREATE | notebooks/full_workflow.ipynb |
| **7-page Paper** | TO CREATE | Based on EVALUATION.md + README.md |

## Deliverables Status

- [x] Code repository (GitHub)
- [x] Data preparation pipeline
- [x] Feature extraction
- [x] Multiple clustering methods
- [x] Hyperparameter optimization script
- [x] Random controls script
- [x] Ablation studies
- [x] Performance benchmarking
- [x] Evaluation framework documentation
- [ ] Jupyter notebook (TO CREATE)
- [ ] Presentation slides (TO CREATE)
- [ ] 7-page paper (TO CREATE)

## Time Estimates

| Task | Time | Notes |
|------|------|-------|
| Data preparation (one-time) | 4 hours | MMseqs2 is slow |
| Feature extraction | 5 minutes | 144 proteins |
| Hyperparameter optimization | 20 minutes | 30 proteins, grid search |
| Random controls | 10 minutes | All proteins |
| Performance benchmarking | 5 minutes | Subset |
| Full method comparison | 15 minutes | All proteins, all methods |
| **Total evaluation time** | **~1 hour** | After setup |

## Next Steps

1. **Run all evaluation scripts** (in order from Phase 3)
2. **Create Jupyter notebook** showing:
   - Data exploration
   - Method comparison
   - Hyperparameter optimization results
   - Random control analysis
   - Performance benchmarks
   - Visualizations
3. **Create presentation** (~15 slides):
   - Problem statement (2 slides)
   - Methods overview (3 slides)
   - Results (5 slides)
   - Hyperparameter optimization (2 slides)
   - Performance benchmarks (1 slide)
   - Conclusions and limitations (2 slides)
4. **Write 7-page paper**:
   - Abstract
   - Introduction
   - Methods
   - Results
   - Discussion
   - Conclusion
   - References

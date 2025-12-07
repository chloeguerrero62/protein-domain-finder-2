# Structure Domain Finder

A computational pipeline for automated protein domain identification using graph-based clustering on 3D structural data.

## Overview

This project implements an unsupervised method for identifying protein domains from 3D structure using:
- **Graph construction** from C-alpha distance matrices
- **Louvain community detection** for domain clustering
- **Evaluation** against Pfam domain annotations

The pipeline processes diverse protein structures selected through sequence clustering to ensure robust generalization.

---

## Quick Start

### 1. Setup Environment
```bash
# Create conda environment
conda env create -f environment.yml
conda activate domain-finder
```

### 2. Run Complete Pipeline
```bash
# Automated workflow (data acquisition)
chmod +x complete_mmseqs2_workflow.sh
./complete_mmseqs2_workflow.sh

# Domain detection pipeline
python scripts/distance_matrix.py      # Compute distance matrices
python scripts/graph_builder.py        # Build KNN graphs
python scripts/louvain_clustering.py   # Detect domains
```

---

## Pipeline Architecture
```
┌─────────────────────────────────────────────────────────────┐
│                    DATA ACQUISITION                         │
└─────────────────────────────────────────────────────────────┘
                            │
    PDB Sequences (952K+) ──┘
                            │
                            ▼
                  ┌──────────────────┐
                  │  Pre-filtering   │  ← Quality control
                  │  (864K sequences)│     (length, resolution, completeness)
                  └──────────────────┘
                            │
                            ▼
                  ┌──────────────────┐
                  │ MMseqs2 Cluster  │  ← 30% sequence identity
                  │  (29K clusters)  │
                  └──────────────────┘
                            │
                            ▼
                  ┌──────────────────┐
                  │ Stratified Sample│  ← By domain count & size
                  │  (~100 proteins) │
                  └──────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────┐
│                  DOMAIN DETECTION                           │
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
              ┌──────────────────────────┐
              │ Parse PDB Structures     │
              │ Extract C-alpha coords   │
              └──────────────────────────┘
                            │
                            ▼
              ┌──────────────────────────┐
              │ Distance Matrices        │  ← Pairwise Euclidean distances
              │ (N × N per protein)      │
              └──────────────────────────┘
                            │
                            ▼
              ┌──────────────────────────┐
              │ KNN Graphs               │  ← k=10, inverse weighting
              │ (residue networks)       │
              └──────────────────────────┘
                            │
                            ▼
              ┌──────────────────────────┐
              │ Louvain Clustering       │  ← Community detection
              │ (domain assignments)     │
              └──────────────────────────┘
                            │
                            ▼
              ┌──────────────────────────┐
              │ Evaluation vs Pfam       │  ← ARI, NMI, F1 scores
              └──────────────────────────┘
```

---

## Installation

### Prerequisites
- Python 3.10+
- MMseqs2
- Conda (recommended)

### Environment Setup
```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/structure-domain-finder.git
cd structure-domain-finder

# Create environment
conda env create -f environment.yml
conda activate domain-finder

# Verify installation
python -c "import Bio, networkx, community; print('✓ All packages installed')"
```

### Key Dependencies
- **Biopython** - PDB structure parsing
- **NetworkX** - Graph construction and manipulation
- **python-louvain** - Community detection
- **scikit-learn** - Evaluation metrics
- **MMseqs2** - Sequence clustering

---

## Usage

### Data Acquisition (Automated)

Select diverse protein structures from PDB:
```bash
./complete_mmseqs2_workflow.sh
```

This script:
1. Downloads PDB sequences (~1.4M sequences)
2. Filters by quality (length, resolution, completeness)
3. Clusters at 30% sequence identity using MMseqs2
4. Selects ~100 diverse representatives
5. Downloads PDB structure files

**Manual steps:** See [Data Acquisition](#data-acquisition-manual) below.

---

### Domain Detection Pipeline

#### 1. Compute Distance Matrices
```bash
python scripts/distance_matrix.py
```

**Input:** PDB structures in `data/selected_structures/`  
**Output:** Distance matrices in `data/distance_matrices/*.npy`  
**Details:** Computes pairwise C-alpha distances for each protein

---

#### 2. Build KNN Graphs
```bash
python scripts/graph_builder.py
```

**Input:** Distance matrices  
**Output:** Graphs in `data/graphs/*.pkl`  
**Details:** Creates k-nearest neighbor graphs (k=10) with inverse distance weighting

---

#### 3. Detect Domains (Louvain Clustering)
```bash
python scripts/louvain_clustering.py
```

**Input:** Protein graphs  
**Output:** Domain assignments in `data/clusters/*.npy`  
**Details:** Applies Louvain algorithm to detect communities (domains)

---

#### 4. Visualize Results (Optional)
```bash
python scripts/visualize_domains.py
```

**Output:** Visualizations in `data/visualizations/`  
**Details:** Creates plots showing domain assignments and sizes

---

## Project Structure
```
structure-domain-finder/
│
├── src/                           # Source code (tracked in Git)
│   ├── features/
│   │   ├── structure_parser.py   # PDB parsing
│   │   ├── distance_matrix.py    # Distance computation
│   │   └── graph_builder.py      # Graph construction
│   └── models/
│       └── louvain_clustering.py # Domain detection
│
├── scripts/                       # Executable scripts (tracked in Git)
│   ├── prefilter.py              # Quality filtering
│   ├── cluster_protein_selection.py  # Stratified sampling
│   ├── download_selected_pdbs.py     # Structure download
│   ├── distance_matrix.py        # Batch distance computation
│   ├── graph_builder.py          # Batch graph building
│   ├── louvain_clustering.py     # Batch clustering
│   └── visualize_domains.py      # Result visualization
│
├── data/                          # Data files (NOT in Git, except CSVs)
│   ├── pdb_clustering/
│   │   └── selected_proteins_mmseqs2.csv  # Selected proteins (tracked)
│   ├── selected_structures/      # PDB files (NOT tracked - 50MB)
│   ├── distance_matrices/        # .npy files (NOT tracked - 50MB)
│   ├── graphs/                   # .pkl files (NOT tracked - 50MB)
│   ├── clusters/                 # Results (NOT tracked)
│   └── results/                  # Summary CSVs (tracked)
│
├── complete_mmseqs2_workflow.sh  # Automated data acquisition
├── environment.yml               # Conda environment
├── .gitignore                    # Excludes large data files
└── README.md                     # This file
```

**Note:** Large binary files (PDB structures, distance matrices, graphs) are excluded from Git via `.gitignore`. They can be regenerated using the scripts.

---

## Data Acquisition (Manual)

If you prefer to run steps individually:

### 1. Download PDB Sequences
```bash
cd data/pdb_clustering
wget https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
gunzip pdb_seqres.txt.gz
```

### 2. Pre-filter Sequences
```bash
python scripts/prefilter.py
```

**Filters applied:**
- Length: 50-2000 residues
- Unknown residues: ≤5% X's
- Resolution: ≤3.0 Å (X-ray only)

### 3. Cluster with MMseqs2
```bash
cd data/pdb_clustering

# Create database
mmseqs createdb pdb_seqres_filtered.fasta pdb_filtered_DB

# Cluster at 30% identity
mmseqs cluster \
    pdb_filtered_DB \
    pdb_clusters_30 \
    tmp \
    --min-seq-id 0.3 \
    -c 0.8 \
    --cov-mode 1 \
    --cluster-mode 2 \
    -s 7.5 \
    --threads 8

# Extract representatives
mmseqs createsubdb pdb_clusters_30 pdb_filtered_DB pdb_representatives_DB
mmseqs convert2fasta pdb_representatives_DB pdb_representatives_30.fasta
```

### 4. Select Diverse Proteins
```bash
python scripts/cluster_protein_selection.py
```

**Stratified by:**
- Domain count (1/2/3/4+ domains)
- Size (small/medium/large)

### 5. Download PDB Structures
```bash
python scripts/download_selected_pdbs.py
```

---

## Methods

### Protein Selection

**Diversity criteria:**
- No two proteins share >30% sequence identity
- Stratified by Pfam domain count (1/2/3/4+ domains)
- Balanced by protein size (50-2000 residues)

**Dataset statistics:**
| Stage | Count |
|-------|-------|
| Raw PDB sequences | 952,113 |
| After quality filtering | 864,635 |
| After clustering (30% ID) | 29,289 |
| With Pfam annotations | 22,109 |
| **Final dataset** | **~100** |

### Domain Detection Algorithm

1. **Extract C-alpha coordinates** from PDB structures
2. **Compute distance matrix** - Pairwise Euclidean distances between all residues
3. **Build KNN graph** - Connect each residue to k=10 nearest neighbors
4. **Edge weighting** - Inverse distance: `weight = 1 / (distance + ε)`
5. **Louvain clustering** - Detect communities in the graph
6. **Domain assignment** - Each community = one domain

### Evaluation

Predicted domains are evaluated against Pfam annotations using:
- **Adjusted Rand Index (ARI)** - Cluster agreement
- **Normalized Mutual Information (NMI)** - Information overlap
- **Boundary F1 Score** - Domain boundary accuracy

---

## Output Files

### Tracked in Git (Small Files)
✅ `data/pdb_clustering/selected_proteins_mmseqs2.csv` - Selected proteins  
✅ `data/results/*.csv` - Summary statistics  
✅ All source code in `src/` and `scripts/`

### NOT Tracked (Large Files - Regenerate Locally)
❌ `data/selected_structures/*.ent` - PDB files (~50 MB total)  
❌ `data/distance_matrices/*.npy` - Distance matrices (~50 MB)  
❌ `data/graphs/*.pkl` - Graph files (~50 MB)  
❌ `data/clusters/*.npy` - Clustering results  
❌ `data/pdb_clustering/*_DB*` - MMseqs2 databases (hundreds of MB)

**To regenerate data:** Run the pipeline scripts as described above.

---

## Troubleshooting

### Import Errors
```bash
# Make sure you're in project root
cd structure-domain-finder

# Scripts should add project root to Python path automatically
python scripts/distance_matrix.py
```

### Missing Data Files
```bash
# Regenerate distance matrices
python scripts/distance_matrix.py

# Regenerate graphs
python scripts/graph_builder.py

# Regenerate clusters
python scripts/louvain_clustering.py
```

### Large File Push Errors

Don't commit large files to Git! They're excluded in `.gitignore`:
```bash
# Verify .gitignore is working
git status  # Should NOT show .npy, .pkl, .ent files

# If large files are staged, remove them:
git rm --cached data/distance_matrices/*.npy
git rm --cached data/graphs/*.pkl
```

---

## Citation

If you use this pipeline, please cite:
```
[Your Name]. (2025). Structure Domain Finder: Graph-based protein domain 
identification using Louvain clustering. GitHub repository. 
https://github.com/YOUR_USERNAME/structure-domain-finder
```

---

## License

MIT License - See LICENSE file for details

---

## Contact

**Author:** Your Name  
**Email:** your.email@example.com  
**Course:** CS-584 Machine Learning (Fall 2025)  
**Institution:** Your University

---

## Acknowledgments

- **RCSB PDB** - Protein structure data
- **Pfam** - Domain annotations
- **MMseqs2** - Sequence clustering
- **NetworkX** - Graph algorithms
- **python-louvain** - Community detection

---

## Future Work

- [ ] Hyperparameter optimization (k, resolution, weighting)
- [ ] Alternative clustering algorithms (Spectral, DBSCAN)
- [ ] Deep learning approaches (GCN, GAT)
- [ ] Extension to multi-chain complexes
- [ ] Web interface for interactive domain prediction
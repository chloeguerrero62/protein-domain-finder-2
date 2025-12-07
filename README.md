# Structure Domain Finder


A pipeline for selecting diverse protein structures from the Protein Data Bank (PDB) for structural domain analysis.

## Quick Start: Run the Complete Workflow

You can run the entire pipeline automatically using the provided shell script:

```bash
chmod +x complete_mmseqs2_workflow.sh
./complete_mmseqs2_workflow.sh
```

This script will:
- Download PDB sequences
- Pre-filter sequences
- Create MMseqs2 database
- Cluster sequences at 30% identity
- Extract cluster representatives
- Create cluster membership table
- Select ~100 representative proteins
- Download PDB structure files

See below for details on each step and manual usage if you want to run steps individually.

## Overview


This project implements a workflow to:
0. Download the full PDB sequence file from RCSB
1. Download and filter PDB sequences by quality criteria
2. Cluster sequences at 30% identity using MMseqs2
3. Select a diverse subset of proteins stratified by domain count and chain length

## Environment Setup

```bash
# Create conda environment
conda env create -f environment.yml
conda activate domain-finder

# Or manually install dependencies
conda create -n mmseqs2 python=3.10 -y
conda activate mmseqs2
conda install -c conda-forge -c bioconda mmseqs2 biopython requests pandas -y
```

## Pipeline Overview


### Automated Workflow
All steps below are run by `complete_mmseqs2_workflow.sh`.

### Manual Workflow Diagram
```
PDB Sequences (952K+)
        │
        ▼
┌───────────────────┐
│  1. Pre-filtering │  ← prefilter.py
│  (Quality control)│
└───────────────────┘
        │
        ▼
  Filtered Sequences (~864K)
        │
        ▼
┌───────────────────┐
│  2. MMseqs2       │  ← mmseqs2 CLI
│  Clustering (30%) │
└───────────────────┘
        │
        ▼
  Cluster Representatives (~29K)
        │
        ▼
┌───────────────────┐
│  3. Stratified    │  ← cluster_protein_selection.py
│  Selection        │
└───────────────────┘
        │
        ▼
  Final Dataset (~136 proteins)
```

### 0. Downloading the Full PDB Sequence File


**Purpose**: Download the complete set of protein and nucleic acid sequences from all PDB entries.

**Automated:**
Handled by `complete_mmseqs2_workflow.sh` (Step 1)

**Manual:**
```bash
cd data/pdb_clustering
wget https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
gunzip pdb_seqres.txt.gz
```

**Input:**
- `pdb_seqres.txt` — Contains all sequences from the PDB in FASTA format

This file is used as the starting point for all subsequent filtering and clustering steps.

---
### Downloading Structure Files for Selected Proteins

**Purpose**: After selection, download the actual 3D structure files (PDB format) for the chosen proteins.

**Automated:**
Handled by `complete_mmseqs2_workflow.sh` (Final step)

**Manual:**
Script: `scripts/download_selected_pdbs.py`

**Input**:
- `selected_proteins_mmseqs2.csv` — List of selected proteins (output from the selection script)

**Output**:
- PDB files for each selected protein, saved in `data/selected_structures/` (default)

**Usage**:
```bash
cd scripts
python download_selected_pdbs.py
```

**How it works**:
1. Reads the CSV file of selected proteins
2. For each protein, downloads the corresponding PDB file from the RCSB PDB using Biopython's `PDBList`
3. Saves all downloaded files in the specified output directory

**Note**: You can change the output directory or file format by editing the script arguments.

---
### Parsing Downloaded Structures

#### `src/features/structure_parser.py`

**Purpose**: Parse downloaded PDB structure files and extract C-alpha (CA) atom coordinates for each residue in a specified chain.

**Key Features:**
- Uses Biopython's `PDBParser` to read PDB files in `.ent` format.
- Extracts CA coordinates and residue numbers for a given PDB ID and chain.
- Handles missing files and missing CA atoms gracefully.

**Usage Example:**
```python
from src.features.structure_parser import ProteinStructureParser
parser = ProteinStructureParser(pdb_dir='data/selected_structures')
coords, residue_ids = parser.parse_structure('1a0b', 'A')
```

**Returns:**
- `coords`: NumPy array of shape (n_residues, 3) with CA coordinates
- `residue_ids`: List of residue numbers

**Typical Use:**
Used as a module in downstream scripts to extract structural features for selected proteins.

---

#### `scripts/parse_all_pdbs.py`

**Purpose**: Batch-parse all selected PDB structures and summarize parsing success/failure for each protein chain.

**How it works:**
1. Loads the list of selected proteins from `selected_proteins_mmseqs2.csv`.
2. Iterates over each protein, using `ProteinStructureParser` to extract CA coordinates.
3. Records the number of residues, success/failure, and any errors for each chain.
4. Saves a summary table to `data/results/parsing_results.csv`.

**Usage:**
```bash
python scripts/parse_all_pdbs.py
```

**Output:**
- `data/results/parsing_results.csv`: Table with columns for PDB chain, number of residues, success status, and error messages (if any).

**Typical Use:**
Run after downloading structures to verify parsing and prepare for downstream structural analysis.

---
## Scripts

---
### New and Updated Scripts (2025)

#### `scripts/louvain_clustering.py`
**Purpose:** Perform community detection (clustering) on protein structure graphs using the Louvain algorithm.
**How it works:**
1. Loads graph representations of protein structures (e.g., from distance matrices or graph files).
2. Applies the Louvain method to detect communities or clusters within each graph.
3. Outputs cluster assignments for each residue or node.
**Usage:**
```bash
python scripts/louvain_clustering.py
```
**Output:**
- Cluster/community assignments for each protein (format depends on implementation, e.g., `.csv`, `.json`).

#### `scripts/visualize_domains.py`
**Purpose:** Visualize protein domains, clusters, or structural features for selected proteins.
**How it works:**
1. Loads structure, domain, or clustering data for selected proteins.
2. Generates visualizations (e.g., domain architectures, cluster colorings, or 3D structure plots).
3. Saves figures or interactive plots for analysis and presentation.
**Usage:**
```bash
python scripts/visualize_domains.py
```
**Output:**
- Figures or plots (e.g., `.png`, `.pdf`, or interactive HTML) in a results or figures directory.

#### `scripts/distance_matrix.py`
**Purpose:** Compute and save C-alpha (CA) distance matrices for each selected protein structure.
**How it works:**
1. Loads the list of selected proteins from `selected_proteins_mmseqs2.csv`.
2. For each protein, loads the corresponding PDB structure using `ProteinStructureParser`.
3. Computes the pairwise distance matrix between all CA atoms (using NumPy).
4. Saves each distance matrix as a `.npy` file (NumPy binary) in `data/distance_matrices/`.
**Usage:**
```bash
python scripts/distance_matrix.py
```
**Output:**
- `data/distance_matrices/<pdb_id>_<chain_id>.npy` for each protein chain.

#### `scripts/graph_builder.py`
**Purpose:** (If implemented) Build graph representations of protein structures, typically from distance matrices or CA coordinates.
**How it works:**
1. Loads distance matrices or coordinates for selected proteins.
2. Constructs graphs (e.g., nodes = residues, edges = spatial proximity or other criteria).
3. Saves graph data for downstream analysis or machine learning.
**Usage:**
```bash
python scripts/graph_builder.py
```
**Output:**
- Graph files (format depends on implementation, e.g., `.gpickle`, `.json`, or `.csv`).

---
### Data and Results Folders

- `data/selected_structures/`: Contains downloaded PDB `.ent` files for selected proteins (used by structure parsing and distance matrix scripts).
- `data/results/`: Stores results such as `parsing_results.csv` and any generated figures.
- `data/processed/`: For processed or intermediate CSVs, e.g., `selected_proteins_mmseqs2.csv`.
- `data/distance_matrices/`: (If used) Contains distance matrices for each protein chain.

---
### Workflow Extension

After downloading and parsing structures, you can now:
1. Compute distance matrices for all selected proteins using `scripts/distance_matrix.py`.
2. (Optionally) Build graph representations for further structural or network analysis using `scripts/graph_builder.py`.

These steps enable downstream structural bioinformatics, clustering, or machine learning tasks on your curated protein set.

---
### Step-by-Step Script Reference

#### 1. Download PDB sequences
**Manual:**
```bash
cd data/pdb_clustering
wget https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
gunzip pdb_seqres.txt.gz
```
**Automated:** Handled by `complete_mmseqs2_workflow.sh` (Step 1)

#### 2. Pre-filter sequences (`scripts/prefilter.py`)
**Purpose:** Pre-filter PDB sequences prior to clustering.
**Usage:**
```bash
cd scripts
python prefilter.py
```
**Outputs:**
- `pdb_seqres_filtered.fasta` — Sequences passing all filters
- `pdb_metadata.csv` — Metadata for passing sequences

#### 3. Create MMseqs2 database
**Manual:**
```bash
cd data/pdb_clustering
mmseqs createdb pdb_seqres_filtered.fasta pdb_filtered_DB
```
**Automated:** Handled by `complete_mmseqs2_workflow.sh` (Step 3)

#### 4. Cluster at 30% identity
**Manual:**
```bash
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
```
**Automated:** Handled by `complete_mmseqs2_workflow.sh` (Step 4)

#### 5. Extract cluster representatives
**Manual:**
```bash
mmseqs createsubdb pdb_clusters_30 pdb_filtered_DB pdb_representatives_DB
mmseqs convert2fasta pdb_representatives_DB pdb_representatives_30.fasta
```
**Automated:** Handled by `complete_mmseqs2_workflow.sh` (Step 5)

#### 6. Create cluster membership table
**Manual:**
```bash
mmseqs createtsv \
        pdb_filtered_DB \
        pdb_filtered_DB \
        pdb_clusters_30 \
        pdb_clusters_30.tsv
```
**Automated:** Handled by `complete_mmseqs2_workflow.sh` (Step 6)

#### 7. Select ~100 representative proteins (`scripts/cluster_protein_selection.py`)
**Purpose:** Select a diverse subset of proteins from cluster representatives.
**Usage:**
```bash
cd data/pdb_clustering
python ../../scripts/cluster_protein_selection.py
```
**Output:**
- `selected_proteins_mmseqs2.csv` — Final selected proteins with annotations

#### 8. Download PDB structures (`scripts/download_selected_pdbs.py`)
**Purpose:** Download the actual 3D structure files (PDB format) for the chosen proteins.
**Usage:**
```bash
cd scripts
python download_selected_pdbs.py
```
**Output:**
- PDB files for each selected protein, saved in `data/selected_structures/` (default)

---
### Parsing Downloaded Structures

#### `src/features/structure_parser.py`
**Purpose:** Parse downloaded PDB structure files and extract C-alpha (CA) atom coordinates for each residue in a specified chain.
**Usage Example:**
```python
from src.features.structure_parser import ProteinStructureParser
parser = ProteinStructureParser(pdb_dir='data/selected_structures')
coords, residue_ids = parser.parse_structure('1a0b', 'A')
```
**Returns:**
- `coords`: NumPy array of shape (n_residues, 3) with CA coordinates
- `residue_ids`: List of residue numbers

#### `scripts/parse_all_pdbs.py`
**Purpose:** Batch-parse all selected PDB structures and summarize parsing success/failure for each protein chain.
**Usage:**
```bash
python scripts/parse_all_pdbs.py
```
**Output:**
- `data/results/parsing_results.csv`: Table with columns for PDB chain, number of residues, success status, and error messages (if any).

### complete_mmseqs2_workflow.sh

**Purpose**: Run the entire pipeline from sequence download to structure retrieval with a single command.

**Usage:**
```bash
chmod +x complete_mmseqs2_workflow.sh
./complete_mmseqs2_workflow.sh
```

**Steps performed:**
1. Download PDB sequences
2. Pre-filter sequences
3. Create MMseqs2 database
4. Cluster at 30% identity
5. Extract cluster representatives
6. Create cluster membership table
7. Select 144 representative proteins
8. Download PDB structures

**Outputs:**
- Cluster representatives: `pdb_representatives_30.fasta`
- Cluster membership: `pdb_clusters_30.tsv`
- Selected proteins: `selected_proteins_mmseqs2.csv`
- PDB structures: `selected_structures/`

### 1. `scripts/prefilter.py`

**Purpose**: Pre-filter PDB sequences prior to clustering.

**Input**: 
- `pdb_seqres.txt` — Full PDB sequence file from RCSB

**Outputs**:
- `pdb_seqres_filtered.fasta` — Sequences passing all filters
- `pdb_metadata.csv` — Metadata for passing sequences

**Filters Applied**:
| Filter | Criteria |
|--------|----------|
| Length | 50–2000 residues |
| Unknown residues | ≤5% "X" characters |
| Resolution | ≤3.0 Å (X-ray structures only) |

**Usage**:
```bash
cd scripts
python prefilter.py
```

**How it works**:
1. Parses input FASTA file sequence by sequence
2. Applies length and unknown residue filters
3. Fetches metadata (resolution, experimental method) from RCSB API using parallelized requests
4. Applies resolution filter for X-ray structures
5. Outputs filtered FASTA and metadata CSV

---

### 2. MMseqs2 Clustering (CLI)

**Purpose**: Cluster sequences at 30% sequence identity to reduce redundancy.

**Commands**:
```bash
cd data/pdb_clustering

# Download PDB sequences
wget https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
gunzip pdb_seqres.txt.gz

# Create MMseqs2 database
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

# Extract representative sequences
mmseqs createsubdb pdb_clusters_30 pdb_filtered_DB pdb_representatives_DB
mmseqs convert2fasta pdb_representatives_DB pdb_representatives_30.fasta

# Extract full cluster membership
mmseqs createtsv pdb_filtered_DB pdb_filtered_DB pdb_clusters_30 pdb_clusters_30.tsv
```

**Clustering Parameters**:
| Parameter | Value | Meaning |
|-----------|-------|---------|
| `--min-seq-id` | 0.3 | 30% sequence identity threshold |
| `-c` | 0.8 | 80% coverage required |
| `--cov-mode` | 1 | Coverage of target sequence |
| `--cluster-mode` | 2 | Greedy set cover algorithm |
| `-s` | 7.5 | Sensitivity |

**Results**:
- 952,113 sequences → 29,289 clusters

---

### 3. `scripts/cluster_protein_selection.py`

**Purpose**: Select a diverse subset of proteins from cluster representatives.

**Input**:
- `pdb_representatives_30.fasta` — Cluster representative sequences
- `pdb_metadata.csv` — Metadata from pre-filtering

**Output**:
- `selected_proteins_mmseqs2.csv` — Final selected proteins with annotations

**Selection Strategy**:

Proteins are selected using **stratified sampling** to ensure diversity:

**By Domain Count (Pfam annotations)**:
| Domains | Target % |
|---------|----------|
| 1 | ~20% |
| 2 | ~30% |
| 3 | ~25% |
| 4+ | ~25% |

**By Chain Length (within each domain bin)**:
| Size | Residues | Target % |
|------|----------|----------|
| Small | 50–300 | 25% |
| Medium | 300–800 | 50% |
| Large | 800–2000 | 25% |

**Usage**:
```bash
cd data/pdb_clustering
python ../../scripts/cluster_protein_selection.py
```

**How it works**:
1. Loads cluster representative sequences
2. Queries PDBe API for Pfam domain annotations (~29K API calls)
3. Bins proteins by domain count and length
4. Performs stratified random sampling from each bin
5. Outputs selected proteins to CSV

**Output Columns**:
| Column | Description |
|--------|-------------|
| pdb_chain | PDB ID + chain (e.g., "101M_A") |
| pdb_id | 4-character PDB ID |
| chain | Chain identifier |
| length | Sequence length |
| n_domains | Number of Pfam domains |
| domain_ids | List of Pfam domain IDs |
| domain_bin | Domain category (1, 2, 3, 4+) |
| length_bin | Size category (small, medium, large) |

---

## Directory Structure

```
structure-domain-finder/
├── README.md
├── environment.yml
├── requirements.txt
├── complete_mmseqs2_workflow.sh
├── data/
│   ├── pdb_clustering/
│   │   ├── pdb_seqres.txt
│   │   ├── pdb_seqres.txt.gz
│   │   ├── pdb_seqres_filtered.fasta
│   │   ├── pdb_seqres_filtered.fasta.gz
│   │   ├── pdb_metadata.csv
│   │   ├── pdb_metadata.csv.gz
│   │   ├── pdb_filtered_DB*
│   │   ├── pdb_clusters_30*
│   │   ├── pdb_representatives_DB*
│   │   ├── pdb_representatives_30.fasta
│   │   ├── selected_proteins_mmseqs2.csv
│   │   └── tmp/
│   ├── processed/
│   │   └── selected_proteins_mmseqs2.csv
│   ├── results/
│   │   ├── parsing_results.csv
│   │   └── figures/
│   └── selected_structures/
│       ├── pdb1eg2.ent
│       ├── pdb1kyq.ent
│       ├── ... (many .ent files)
├── docs/
│   └── .gitkeep
├── notebooks/
│   └── .gitkeep
├── scripts/
│   ├── prefilter.py
│   ├── cluster_protein_selection.py
│   ├── download_selected_pdbs.py
│   ├── parse_all_pdbs.py
│   ├── pfam_annotations.json
│   └── selected_structures/ (empty)
├── src/
│   ├── __init__.py
│   ├── data/
│   │   └── prefilter.py
│   ├── evaluation/
│   │   └── __init__.py
│   ├── features/
│   │   ├── __init__.py
│   │   └── structure_parser.py
│   └── models/
│       └── __init__.py
└── .gitignore
```

---

## Results Summary

| Stage | Count |
|-------|-------|
| Raw PDB sequences | ~1.4M |
| After quality filtering | ~864,635 |
| After clustering (30% identity) | 29,289 clusters |
| With Pfam annotations | 22,109 |
| Final selected proteins | 136 |

---

## Dependencies

- Python 3.10+
- MMseqs2
- Biopython
- pandas
- requests

See `environment.yml` for full dependency list.

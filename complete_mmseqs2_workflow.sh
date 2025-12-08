#!/bin/bash
# complete_mmseqs2_workflow.sh

set -e  # Exit on error

echo "=== MMseqs2 Clustering Workflow ==="
echo "Date: $(date)"
echo

# Step 1: Download sequences
echo "Step 1: Downloading PDB sequences..."
if [ ! -f "pdb_seqres.txt" ]; then
    wget https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
    gunzip pdb_seqres.txt.gz
fi

# Step 2: Pre-filter
echo "Step 2: Pre-filtering sequences..."
python prefilter_pdb_sequences.py

# Step 3: Create database
echo "Step 3: Creating MMseqs2 database..."
mmseqs createdb pdb_seqres_filtered.fasta pdb_filtered_DB

# Step 4: Cluster
echo "Step 4: Clustering at 30% identity..."
mkdir -p tmp
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

# Step 5: Extract representatives
echo "Step 5: Extracting cluster representatives..."
mmseqs createsubdb pdb_clusters_30 pdb_filtered_DB pdb_representatives_DB
mmseqs convert2fasta pdb_representatives_DB pdb_representatives_30.fasta

# Step 6: Create cluster TSV
echo "Step 6: Creating cluster membership table..."
mmseqs createtsv \
    pdb_filtered_DB \
    pdb_filtered_DB \
    pdb_clusters_30 \
    pdb_clusters_30.tsv

# Step 7: Select 144 proteins
echo "Step 7: Selecting 144 representative proteins..."
python select_proteins_from_clusters.py

# Step 8: Download structures
echo "Step 8: Downloading PDB structures..."
python download_selected_pdbs.py

echo
echo "=== Workflow Complete ==="
echo "Results:"
echo "  - Cluster representatives: pdb_representatives_30.fasta"
echo "  - Cluster membership: pdb_clusters_30.tsv"
echo "  - Selected 144 proteins: selected_proteins_mmseqs2.csv"
echo "  - PDB structures: selected_structures/"

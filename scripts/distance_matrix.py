# scripts/distance_matrix.py
"""
Compute distance matrices for all selected proteins
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from src.features.structure_parser import ProteinStructureParser
from src.features.distance_matrix import compute_distance_matrix

# Setup
output_dir = Path('data/distance_matrices')
output_dir.mkdir(parents=True, exist_ok=True)

# Load dataset
df = pd.read_csv('data/pdb_clustering/selected_proteins_mmseqs2.csv')
print(f"Processing {len(df)} proteins...\n")

# Create parser
parser = ProteinStructureParser(pdb_dir='data/selected_structures')

# Track results
results = []
successful = 0
failed = 0

for i, row in df.iterrows():
    pdb_id = row['pdb_id']
    chain_id = row['chain']
    pdb_chain = f"{pdb_id}_{chain_id}"
    
    print(f"[{i+1}/{len(df)}] {pdb_chain}...", end=' ')
    
    try:
        # Parse structure (gets coords fresh each time)
        coords, res_ids = parser.parse_structure(pdb_id, chain_id)
        
        # Compute distance matrix
        D = compute_distance_matrix(coords)
        
        # Save distance matrix
        output_file = output_dir / f'{pdb_id}_{chain_id}_distmat.npy'
        np.save(output_file, D)
        
        successful += 1
        results.append({
            'pdb_chain': pdb_chain,
            'success': True,
            'n_residues': len(coords),
            'matrix_shape': D.shape
        })
        
        print(f"✓ {D.shape}")
        
    except FileNotFoundError as e:
        failed += 1
        results.append({
            'pdb_chain': pdb_chain,
            'success': False,
            'error': 'File not found'
        })
        print(f"✗ Missing PDB file")
        
    except Exception as e:
        failed += 1
        results.append({
            'pdb_chain': pdb_chain,
            'success': False,
            'error': str(e)
        })
        print(f"✗ {e}")

# Summary
print(f"\n{'='*60}")
print(f"SUMMARY")
print(f"{'='*60}")
print(f"Successful: {successful}/{len(df)}")
print(f"Failed: {failed}/{len(df)}")

if failed > 0:
    print(f"\nFailed proteins:")
    for r in results:
        if not r['success']:
            print(f"  {r['pdb_chain']}: {r.get('error', 'Unknown error')}")

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv('data/results/distance_matrix_results.csv', index=False)
print(f"\nResults saved to: data/results/distance_matrix_results.csv")
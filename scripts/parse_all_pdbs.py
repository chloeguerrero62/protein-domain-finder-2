import sys
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from src.features.structure_parser import ProteinStructureParser
import pandas as pd

# Load dataset
df = pd.read_csv('data/pdb_clustering/selected_proteins_mmseqs2.csv')
print(f"Processing {len(df)} proteins...\n")

# Create parser
parser = ProteinStructureParser(pdb_dir='data/selected_structures')

# Store results
results = []

# Loop through each protein
for i, row in df.iterrows():
    pdb_id = row['pdb_id']
    chain_id = row['chain']
    
    print(f"[{i+1}/{len(df)}] Processing {pdb_id}_{chain_id}...", end=' ')
    
    try:
        # Parse structure
        coords, res_ids = parser.parse_structure(pdb_id, chain_id)
        
        # Store results
        results.append({
            'pdb_chain': f"{pdb_id}_{chain_id}",
            'pdb_id': pdb_id,
            'chain': chain_id,
            'n_residues': len(coords),
            'success': True,
            'error': None
        })
        
        print(f"✓ ({len(coords)} residues)")
        
    except Exception as e:
        # Handle errors gracefully
        results.append({
            'pdb_chain': f"{pdb_id}_{chain_id}",
            'pdb_id': pdb_id,
            'chain': chain_id,
            'n_residues': 0,
            'success': False,
            'error': str(e)
        })
        
        print(f"✗ ERROR: {e}")

# Summary
print(f"\n{'='*60}")
print(f"SUMMARY")
print(f"{'='*60}")
successful = sum(1 for r in results if r['success'])
failed = sum(1 for r in results if not r['success'])
print(f"Successful: {successful}/{len(df)}")
print(f"Failed: {failed}/{len(df)}")

if failed > 0:
    print(f"\nFailed proteins:")
    for r in results:
        if not r['success']:
            print(f"  {r['pdb_chain']}: {r['error']}")

# Save results to CSV
results_df = pd.DataFrame(results)
results_df.to_csv('data/results/parsing_results.csv', index=False)
print(f"\nResults saved to: results/parsing_results.csv")
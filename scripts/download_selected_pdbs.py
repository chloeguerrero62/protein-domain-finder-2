"""
Download PDB files for selected proteins
"""

import pandas as pd
from Bio.PDB import PDBList

def download_selected_structures(csv_file, output_dir="data/selected_structures"):
    """
    Download PDB structures for selected proteins
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    df = pd.read_csv(csv_file)
    
    pdbl = PDBList()
        
    for i, row in df.iterrows():
        pdb_id = row['pdb_id'].lower()
        
        print(f"  [{i+1}/{len(df)}] {pdb_id}")
        
        try:
            pdbl.retrieve_pdb_file(
                pdb_id,
                pdir=output_dir,
                file_format='pdb'
            )
        except Exception as e:
            print(f"    ERROR: {e}")


if __name__ == "__main__":
    download_selected_structures("data/pdb_clustering/selected_proteins_mmseqs2.csv")
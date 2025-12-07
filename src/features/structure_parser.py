"""
Parse PDB structures and extract C-alpha coordinates
"""

from Bio.PDB import PDBParser
import numpy as np
from pathlib import Path

class ProteinStructureParser:
    """Parse PDB files and extract structural features"""
    
    def __init__(self, pdb_dir="data/selected_structures"):
        self.pdb_dir = Path(pdb_dir)
        self.parser = PDBParser(QUIET=True, PERMISSIVE=True)
    
    def parse_structure(self, pdb_id, chain_id):
        """
        Extract C-alpha coordinates from PDB structure
        
        Parameters:
        -----------
        pdb_id : str
            4-letter PDB ID (e.g., '1a0b')
        chain_id : str
            Chain identifier (e.g., 'A')
        
        Returns:
        --------
        coords : np.ndarray, shape (n_residues, 3)
            C-alpha coordinates
        residue_ids : list
            Residue numbers
        """
        # Find PDB file (.ent format from your downloads)
        pdb_file = self.pdb_dir / f"pdb{pdb_id.lower()}.ent"
        
        if not pdb_file.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        # Parse structure
        structure = self.parser.get_structure(pdb_id, pdb_file)
        
        coords = []
        residue_ids = []
        
        for model in structure:
            for chain in model:
                if chain.id != chain_id:
                    continue
                
                for residue in chain:
                    # Skip heteroatoms
                    if residue.id[0] != ' ':
                        continue
                    
                    if 'CA' in residue:
                        coords.append(residue['CA'].coord)
                        residue_ids.append(residue.id[1])
            
            break  # Only first model
        
        if len(coords) == 0:
            raise ValueError(f"No C-alpha atoms found")
        
        return np.array(coords), residue_ids
"""
Select 100 representative proteins from MMseqs2 clusters
Stratified by:
1. Pfam domain count (1, 2, 3, 4+)
2. Chain length (small, medium, large)
"""

import pandas as pd
import numpy as np
from collections import defaultdict

def load_cluster_representatives(fasta_file):
    """
    Load representative sequences from FASTA
    Returns dict: {pdb_chain: sequence}
    """
    from Bio import SeqIO
    
    representatives = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        pdb_chain = record.id
        representatives[pdb_chain] = str(record.seq)
    
    return representatives

def get_pfam_domains(pdb_id, chain_id):
    """
    Query Pfam annotations for a PDB chain
    Returns: number of domains, list of domain IDs
    (aka counts the number of unique pfam domains are present in a given chain)
    """
    import requests
    
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{pdb_id}"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            
            if pdb_id.lower() in data:
                pfam_mappings = data[pdb_id.lower()]['Pfam']
                
                # Get domains for specific chain
                chain_domains = []
                for domain_id, domain_info in pfam_mappings.items():
                    for mapping in domain_info['mappings']:
                        if mapping['chain_id'] == chain_id:
                            chain_domains.append(domain_id)
                
                return len(set(chain_domains)), list(set(chain_domains))
    except:
        pass
    
    return None, None

def stratified_selection(representatives, metadata_df, target_n=100):
    """
    Selects 100 proteins stratified by domain count and length
    
    Target distribution
    - 1 domain: 20% (20 proteins)
    - 2 domains: 30% (30 proteins)
    - 3 domains: 25% (25 proteins)
    - 4+ domains: 25% (25 proteins)
    
    Within each domain bin, stratify by length:
    - Small (50-300): 25%
    - Medium (300-800): 50%
    - Large (800-2000): 25%
    """
    print("Fetching Pfam domain annotations (29,289 API calls): ")
    
    # Annotate with Pfam domains
    annotated = []
    
    for i, (pdb_chain, seq) in enumerate(representatives.items()):
        if i % 1000 == 0:
            print(f"  Progress: {i}/{len(representatives)}")
        
        pdb_id = pdb_chain.split('_')[0]
        chain_id = pdb_chain.split('_')[1] if '_' in pdb_chain else 'A'
        
        n_domains, domain_ids = get_pfam_domains(pdb_id, chain_id)
        
        if n_domains is None or n_domains == 0:
            # No Pfam annotations - skip
            continue
        
        length = len(seq)
        
        # Bin domain count
        if n_domains == 1:
            domain_bin = '1'
        elif n_domains == 2:
            domain_bin = '2'
        elif n_domains == 3:
            domain_bin = '3'
        else:
            domain_bin = '4+'
        
        # Bin length
        if length < 300:
            length_bin = 'small'
        elif length < 800:
            length_bin = 'medium'
        else:
            length_bin = 'large'
        
        annotated.append({
            'pdb_chain': pdb_chain,
            'pdb_id': pdb_id,
            'chain': chain_id,
            'length': length,
            'n_domains': n_domains,
            'domain_ids': domain_ids,
            'domain_bin': domain_bin,
            'length_bin': length_bin
        })
    
    df = pd.DataFrame(annotated)
    
    print(f"\nAnnotated {len(df)} proteins with Pfam domains")
    print("\nDomain count distribution:")
    print(df['domain_bin'].value_counts().sort_index())
    print("\nLength distribution:")
    print(df['length_bin'].value_counts())
    
    # Define target numbers
    target_dist = {
        '1': 29,
        '2': 43,
        '3': 36,
        '4+': 36
    }
    
    # Stratified sampling
    selected = []
    
    for domain_bin, target_count in target_dist.items():
        # Get proteins in this domain bin
        domain_df = df[df['domain_bin'] == domain_bin]
        
        if len(domain_df) < target_count:
            print(f"WARNING: Only {len(domain_df)} proteins with "
                  f"{domain_bin} domains (need {target_count})")
            selected.extend(domain_df.to_dict('records'))
            continue
        
        # Stratify by length within domain bin
        length_targets = {
            'small': int(target_count * 0.25),
            'medium': int(target_count * 0.50),
            'large': int(target_count * 0.25)
        }
        
        for length_bin, length_count in length_targets.items():
            bin_df = domain_df[domain_df['length_bin'] == length_bin]
            
            if len(bin_df) < length_count:
                print(f"WARNING: Only {len(bin_df)} proteins with "
                      f"{domain_bin} domains, {length_bin} size")
                selected.extend(bin_df.to_dict('records'))
            else:
                # Random sample
                sample = bin_df.sample(n=length_count, random_state=42)
                selected.extend(sample.to_dict('records'))
    
    selected_df = pd.DataFrame(selected)
    
    print(f"\n=== FINAL SELECTION ===")
    print(f"Total proteins: {len(selected_df)}")
    print("\nDomain distribution:")
    print(selected_df['domain_bin'].value_counts().sort_index())
    print("\nLength distribution:")
    print(selected_df.groupby('domain_bin')['length_bin'].value_counts())
    
    return selected_df

if __name__ == "__main__":
    # Load representatives
    print("Loading cluster representatives...")
    reps = load_cluster_representatives("pdb_representatives_30.fasta")
    print(f"Loaded {len(reps)} representatives")
    
    # Load metadata
    metadata = pd.read_csv("pdb_metadata.csv")
    
    # Select 100 proteins
    selected = stratified_selection(reps, metadata, target_n=100)
    
    # Save
    selected.to_csv("selected_proteins_mmseqs2.csv", index=False)
    print("\nSaved to: selected_proteins_mmseqs2.csv")
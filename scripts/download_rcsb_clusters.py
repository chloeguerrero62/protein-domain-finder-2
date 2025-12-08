"""
Download RCSB Pre-Clustered Sequences at 30% Identity

This script replaces the manual MMseqs2 clustering workflow by directly
using RCSB's pre-computed sequence clusters at 30% identity.

RCSB provides clustering data at: 
- 30%, 40%, 50%, 70%, 90%, 95%, 100% sequence identity levels
- Updates weekly with new PDB releases

Workflow:
1. Download cluster representatives from RCSB
2. Apply quality filters (length, resolution)
3. Sample diverse proteins by Pfam domain annotations
4. Download PDB structures for selected proteins
"""

import requests
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
from io import StringIO
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
import json
from tqdm import tqdm


# =============================================================================
# STEP 1: DOWNLOAD RCSB PRE-CLUSTERED DATA
# =============================================================================

def download_rcsb_clusters(identity_level=30, output_dir='data/pdb_clustering'):
    """
    Download RCSB pre-clustered sequences
    
    RCSB provides cluster representatives at various identity levels.
    We use the GraphQL API to get cluster information.
    
    Parameters:
    -----------
    identity_level : int
        Sequence identity level (30, 40, 50, 70, 90, 95, 100)
    output_dir : str
        Output directory
    
    Returns:
    --------
    clusters_df : pd.DataFrame
        DataFrame with cluster information
    """
    
    print("="*70)
    print(f"DOWNLOADING RCSB PRE-CLUSTERED DATA ({identity_level}% identity)")
    print("="*70)
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # RCSB provides clustering data via their data API
    # We'll query for all PDB entries and their cluster memberships
    
    print("\nFetching cluster data from RCSB...")
    print("Note: This may take several minutes for the full PDB database")
    
    # Query RCSB Search API for all entries with clustering info
    query = {
        "query": {
            "type": "group",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.ncbi_taxonomy_id",
                        "operator": "exists"
                    }
                }
            ],
            "logical_operator": "and"
        },
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "entry"
    }
    
    # Get all PDB IDs first
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    try:
        response = requests.post(search_url, json=query, timeout=60)
        if response.status_code == 200:
            data = response.json()
            all_pdb_ids = [hit['identifier'] for hit in data.get('result_set', [])]
            print(f"✓ Found {len(all_pdb_ids)} PDB entries")
        else:
            print(f"⚠ Search API failed (status {response.status_code})")
            print("  Falling back to manual sequence download...")
            return download_sequences_fallback(output_dir)
    except Exception as e:
        print(f"⚠ Search API error: {e}")
        print("  Falling back to manual sequence download...")
        return download_sequences_fallback(output_dir)
    
    # Get entity-level clustering information
    # RCSB clusters at the entity (chain) level
    print("\nFetching cluster memberships...")
    
    # We'll use the data API to get cluster info in batches
    cluster_data = []
    batch_size = 100
    
    for i in tqdm(range(0, len(all_pdb_ids), batch_size), desc="Batches"):
        batch = all_pdb_ids[i:i+batch_size]
        
        for pdb_id in batch:
            # Query for entity sequences and cluster membership
            entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}"
            
            try:
                resp = requests.get(entity_url, timeout=5)
                if resp.status_code == 200:
                    entity_data = resp.json()
                    
                    # Extract sequence clusters (if available)
                    if 'rcsb_cluster_membership' in entity_data:
                        clusters = entity_data['rcsb_cluster_membership']
                        
                        for cluster in clusters:
                            if cluster.get('identity') == identity_level:
                                cluster_data.append({
                                    'pdb_id': pdb_id.lower(),
                                    'entity_id': entity_data.get('rcsb_id'),
                                    'cluster_id': cluster.get('cluster_id'),
                                    'identity': cluster.get('identity'),
                                    'sequence': entity_data.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can')
                                })
            except:
                continue
    
    print(f"✓ Retrieved clustering data for {len(cluster_data)} entities")
    
    # Convert to DataFrame
    clusters_df = pd.DataFrame(cluster_data)
    
    # Save cluster data
    output_file = output_path / f'rcsb_clusters_{identity_level}.csv'
    clusters_df.to_csv(output_file, index=False)
    
    print(f"✓ Saved cluster data to: {output_file}")
    
    return clusters_df


def download_sequences_fallback(output_dir='data/pdb_clustering'):
    """
    Fallback: Download PDB sequence database and parse cluster information
    
    RCSB provides a sequences file with cluster annotations in the headers.
    This is more reliable than the API for large-scale downloads.
    """
    
    print("\nUsing fallback method: Downloading PDB sequence file...")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Download the full PDB sequence file (includes cluster info in headers)
    seq_url = "https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
    
    import urllib.request
    import gzip
    
    seq_file = output_path / 'pdb_seqres.txt.gz'
    
    if not seq_file.exists():
        print(f"Downloading: {seq_url}")
        print("This is a ~2 GB file and may take several minutes...")
        
        urllib.request.urlretrieve(seq_url, seq_file)
        print("✓ Download complete")
    else:
        print(f"✓ Using existing file: {seq_file}")
    
    # Parse sequences (cluster info is in headers)
    print("\nParsing sequences...")
    
    sequences = []
    
    with gzip.open(seq_file, 'rt') as f:
        for record in SeqIO.parse(f, 'fasta'):
            # Header format: >101M_A mol:protein length:154  MYOGLOBIN
            header_parts = record.description.split()
            
            pdb_chain = record.id
            pdb_id = pdb_chain.split('_')[0]
            chain = pdb_chain.split('_')[1] if '_' in pdb_chain else 'A'
            
            length = len(record.seq)
            
            sequences.append({
                'pdb_id': pdb_id.lower(),
                'chain': chain,
                'pdb_chain': pdb_chain,
                'length': length,
                'sequence': str(record.seq),
                'description': record.description
            })
    
    print(f"✓ Parsed {len(sequences)} sequences")
    
    # Note: The pdb_seqres file doesn't include explicit cluster IDs
    # We'll need to use RCSB's clustering API or compute clusters ourselves
    # For now, save sequences and move to next step
    
    seq_df = pd.DataFrame(sequences)
    seq_file_out = output_path / 'pdb_sequences.csv'
    seq_df.to_csv(seq_file_out, index=False)
    
    print(f"✓ Saved sequences to: {seq_file_out}")
    
    return seq_df


# =============================================================================
# STEP 2: GET CLUSTER REPRESENTATIVES
# =============================================================================

def get_cluster_representatives_api(identity_level=30):
    """
    Use RCSB Cluster API to get representatives directly
    
    This is the most reliable method - RCSB provides cluster representatives
    via their sequence cluster API.
    """
    
    print("\n" + "="*70)
    print("FETCHING CLUSTER REPRESENTATIVES FROM RCSB API")
    print("="*70)
    
    # RCSB cluster API endpoint
    # Format: clusters-by-entity-{identity}-identity.txt
    cluster_url = "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-30.txt"

    print(f"\nDownloading from: {cluster_url}")
    
    try:
        response = requests.get(cluster_url, timeout=30)
        
        if response.status_code != 200:
            print(f"⚠ Failed to download (status {response.status_code})")
            return None
        
        print("✓ Download successful")
        
        # Parse cluster file
        # Format: Each line is a cluster with representative first
        # Example: 1A0A_1 1A0B_1 1A0C_1 ...
        
        clusters = []
        representatives = []
        
        for line in response.text.strip().split('\n'):
            entities = line.split()
            
            if entities:
                # First entity is the representative
                rep = entities[0]
                cluster_members = entities
                
                clusters.append({
                    'representative': rep,
                    'n_members': len(cluster_members),
                    'members': cluster_members
                })
                
                representatives.append(rep)
        
        print(f"\n✓ Parsed {len(clusters)} clusters")
        print(f"✓ {len(representatives)} cluster representatives")
        
        return clusters, representatives
        
    except Exception as e:
        print(f"⚠ Error: {e}")
        return None


# =============================================================================
# STEP 3: FETCH METADATA FOR REPRESENTATIVES
# =============================================================================

def get_pdb_metadata_batch(pdb_ids):
    """
    Fetch metadata for multiple PDB entries
    Returns resolution, experimental method, release date
    """
    
    metadata = {}
    
    def fetch_single(pdb_id):
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        
        try:
            response = requests.get(url, timeout=5)
            if response.status_code == 200:
                data = response.json()
                
                exp_method = data.get('exptl', [{}])[0].get('method', 'Unknown')
                
                resolution = None
                if 'rcsb_entry_info' in data:
                    resolution = data['rcsb_entry_info'].get('resolution_combined', [None])[0]
                
                release_date = data.get('rcsb_accession_info', {}).get('initial_release_date')
                
                return pdb_id, {
                    'exp_method': exp_method,
                    'resolution': resolution,
                    'release_date': release_date
                }
        except:
            pass
        
        return pdb_id, None
    
    print(f"\nFetching metadata for {len(pdb_ids)} PDB entries...")
    
    with ThreadPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(fetch_single, pdb_id): pdb_id for pdb_id in pdb_ids}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Metadata"):
            pdb_id, meta = future.result()
            if meta:
                metadata[pdb_id] = meta
    
    print(f"✓ Retrieved metadata for {len(metadata)} entries")
    
    return metadata


def get_entity_sequences(entity_ids):
    """
    Fetch sequences for specific entity IDs
    Entity ID format: PDB_ID_ENTITY_NUMBER (e.g., 1A0A_1, 1A0A_2)
    """
    
    sequences = {}
    
    def fetch_single(entity_id):
        # Parse entity ID: PDB_ENTITY (e.g., 1A0A_1)
        parts = entity_id.split('_')
        if len(parts) != 2:
            return entity_id, None
        
        pdb_id = parts[0]
        entity_num = parts[1]
        
        url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_num}"
        
        try:
            response = requests.get(url, timeout=5)
            if response.status_code == 200:
                data = response.json()
                
                sequence = data.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', '')
                chains = data.get('entity_poly', {}).get('pdbx_strand_id', '').split(',')
                
                return entity_id, {
                    'sequence': sequence,
                    'length': len(sequence),
                    'chains': chains,
                    'pdb_id': pdb_id
                }
        except:
            pass
        
        return entity_id, None
    
    print(f"\nFetching sequences for {len(entity_ids)} entities...")
    
    with ThreadPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(fetch_single, eid): eid for eid in entity_ids}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Sequences"):
            entity_id, seq_data = future.result()
            if seq_data:
                sequences[entity_id] = seq_data
    
    print(f"✓ Retrieved sequences for {len(sequences)} entities")
    
    return sequences


# =============================================================================
# STEP 4: FILTER AND SELECT PROTEINS
# =============================================================================

def filter_representatives(representatives, sequences, metadata):
    """
    Apply quality filters to cluster representatives
    
    Filters:
    1. Length: 50-2000 residues
    2. Resolution: ≤3.0 Å (for X-ray structures)
    3. Completeness: ≤5% unknown residues
    """
    
    print("\n" + "="*70)
    print("APPLYING QUALITY FILTERS")
    print("="*70)
    
    filtered = []
    
    for entity_id in tqdm(representatives, desc="Filtering"):
        if entity_id not in sequences:
            continue
        
        seq_data = sequences[entity_id]
        pdb_id = seq_data['pdb_id']
        sequence = seq_data['sequence']
        length = seq_data['length']
        
        # Filter 1: Length
        if length < 50 or length > 2000:
            continue
        
        # Filter 2: Unknown residues
        x_count = sequence.count('X')
        if x_count / length > 0.05:
            continue
        
        # Filter 3: Resolution (if X-ray)
        if pdb_id in metadata:
            meta = metadata[pdb_id]
            exp_method = meta.get('exp_method', '')
            resolution = meta.get('resolution')
            
            if 'X-RAY' in exp_method.upper():
                if resolution is None or resolution > 3.0:
                    continue
        else:
            # No metadata - skip
            continue
        
        # Passed all filters
        filtered.append({
            'entity_id': entity_id,
            'pdb_id': pdb_id,
            'chains': seq_data['chains'],
            'length': length,
            'sequence': sequence,
            'exp_method': meta.get('exp_method'),
            'resolution': meta.get('resolution'),
            'release_date': meta.get('release_date')
        })
    
    print(f"\n✓ Filtered to {len(filtered)} high-quality proteins")
    
    return pd.DataFrame(filtered)


# =============================================================================
# STEP 5: STRATIFIED SAMPLING WITH PFAM ANNOTATIONS
# =============================================================================

def get_pfam_domains(pdb_id, chain_id):
    """Query Pfam annotations for a PDB chain"""
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{pdb_id}"
    
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            
            if pdb_id.lower() in data:
                pfam_mappings = data[pdb_id.lower()]['Pfam']
                
                chain_domains = []
                for domain_id, domain_info in pfam_mappings.items():
                    for mapping in domain_info['mappings']:
                        if mapping['chain_id'] == chain_id:
                            chain_domains.append(domain_id)
                
                return len(set(chain_domains)), list(set(chain_domains))
    except:
        pass
    
    return None, None


def stratified_selection_with_pfam(filtered_df, target_n=144):
    """
    Select diverse proteins stratified by domain count and length
    
    Target distribution:
    - 1 domain: 20%
    - 2 domains: 30%
    - 3 domains: 25%
    - 4+ domains: 25%
    
    Within each bin, stratify by size:
    - Small (50-300): 25%
    - Medium (300-800): 50%
    - Large (800-2000): 25%
    """
    
    print("\n" + "="*70)
    print("FETCHING PFAM DOMAIN ANNOTATIONS")
    print("="*70)
    
    # Annotate with Pfam domains
    annotated = []
    
    for i, row in tqdm(filtered_df.iterrows(), total=len(filtered_df), desc="Pfam"):
        pdb_id = row['pdb_id']
        chains = row['chains']
        
        # Use first chain
        if isinstance(chains, list) and len(chains) > 0:
            chain_id = chains[0]
        elif isinstance(chains, str):
            chain_id = chains.split(',')[0] if ',' in chains else chains
        else:
            chain_id = 'A'
        
        n_domains, domain_ids = get_pfam_domains(pdb_id, chain_id)
        
        if n_domains is None or n_domains == 0:
            continue
        
        length = row['length']
        
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
            'pdb_chain': f"{pdb_id}_{chain_id}",
            'pdb_id': pdb_id,
            'chain': chain_id,
            'length': length,
            'n_domains': n_domains,
            'domain_ids': str(domain_ids),
            'domain_bin': domain_bin,
            'length_bin': length_bin,
            'exp_method': row['exp_method'],
            'resolution': row['resolution'],
            'release_date': row['release_date']
        })
    
    df = pd.DataFrame(annotated)
    
    print(f"\n✓ Annotated {len(df)} proteins with Pfam domains")
    print("\nDomain count distribution:")
    print(df['domain_bin'].value_counts().sort_index())
    print("\nLength distribution:")
    print(df['length_bin'].value_counts())
    
    # Stratified sampling
    target_dist = {
        '1': int(target_n * 0.20),
        '2': int(target_n * 0.30),
        '3': int(target_n * 0.25),
        '4+': int(target_n * 0.25)
    }
    
    selected = []
    
    for domain_bin, target_count in target_dist.items():
        domain_df = df[df['domain_bin'] == domain_bin]
        
        if len(domain_df) < target_count:
            print(f"⚠ Only {len(domain_df)} proteins with {domain_bin} domains (need {target_count})")
            selected.extend(domain_df.to_dict('records'))
            continue
        
        # Stratify by length
        length_targets = {
            'small': int(target_count * 0.25),
            'medium': int(target_count * 0.50),
            'large': int(target_count * 0.25)
        }
        
        for length_bin, length_count in length_targets.items():
            bin_df = domain_df[domain_df['length_bin'] == length_bin]
            
            if len(bin_df) < length_count:
                print(f"⚠ Only {len(bin_df)} proteins with {domain_bin} domains, {length_bin} size")
                selected.extend(bin_df.to_dict('records'))
            else:
                sample = bin_df.sample(n=length_count, random_state=42)
                selected.extend(sample.to_dict('records'))
    
    selected_df = pd.DataFrame(selected)
    
    print(f"\n{'='*70}")
    print("FINAL SELECTION")
    print("="*70)
    print(f"Total proteins: {len(selected_df)}")
    print("\nDomain distribution:")
    print(selected_df['domain_bin'].value_counts().sort_index())
    print("\nLength distribution:")
    print(selected_df.groupby('domain_bin')['length_bin'].value_counts())
    
    return selected_df


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """
    Complete workflow using RCSB pre-clustered data
    """
    
    print("="*70)
    print("PROTEIN SELECTION FROM RCSB PRE-CLUSTERED DATA")
    print("="*70)
    
    output_dir = Path('data/pdb_clustering')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Get cluster representatives from RCSB
    result = get_cluster_representatives_api(identity_level=30)
    
    if result is None:
        print("\n⚠ Failed to get cluster data from API")
        print("Please check your internet connection or try again later.")
        return
    
    clusters, representatives = result
    
    # Save cluster info
    cluster_file = output_dir / 'rcsb_clusters_30.json'
    with open(cluster_file, 'w') as f:
        json.dump(clusters, f, indent=2)
    print(f"\n✓ Saved cluster data to: {cluster_file}")
    
    # Step 2: Get sequences for representatives
    sequences = get_entity_sequences(representatives[:1000])  # Limit for testing
    
    # Step 3: Get metadata for PDB IDs
    unique_pdb_ids = list(set([seq['pdb_id'] for seq in sequences.values()]))
    metadata = get_pdb_metadata_batch(unique_pdb_ids)
    
    # Step 4: Filter by quality criteria
    filtered_df = filter_representatives(representatives, sequences, metadata)
    
    # Save filtered representatives
    filtered_file = output_dir / 'filtered_representatives.csv'
    filtered_df.to_csv(filtered_file, index=False)
    print(f"✓ Saved filtered representatives to: {filtered_file}")
    
    # Step 5: Stratified selection with Pfam annotations
    selected_df = stratified_selection_with_pfam(filtered_df, target_n=144)
    
    # Save final selection
    output_file = output_dir / 'selected_proteins_mmseqs2.csv'
    selected_df.to_csv(output_file, index=False)
    
    print(f"\n{'='*70}")
    print("COMPLETE")
    print("="*70)
    print(f"✓ Final selection saved to: {output_file}")
    print(f"✓ Total proteins selected: {len(selected_df)}")
    print(f"\nNext step: Run download_selected_pdbs.py to get structures")


if __name__ == "__main__":
    main()
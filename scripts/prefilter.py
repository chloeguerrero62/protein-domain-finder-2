'''
pre-filtration of pdb files prior to clustering with mmseqs2
filters pertain to chain length, resolutin, and completeness
'''

import gzip
from Bio import SeqIO
from Bio.PDB import PDBList, PDBParser, MMCIFParser
import requests
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

################################################################################
#extracts the experimental approach used to determine the structure, the resolution of the structure in 
# angstroms, and the release date; returns dictionary with these categories as keys
def get_pdb_metadata(pdb_id):
    """
    Fetch PDB metadata from RCSB API
    Returns: resolution, exp_method, release_date
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            
            # Extract metadata
            exp_method = data.get('exptl', [{}])[0].get('method', 'Unknown')
            
            # Resolution (only for X-ray)
            resolution = None
            if 'rcsb_entry_info' in data:
                resolution = data['rcsb_entry_info'].get('resolution_combined', [None])[0]
            
            release_date = data.get('rcsb_accession_info', {}).get('initial_release_date')
            
            return {
                'exp_method': exp_method,
                'resolution': resolution,
                'release_date': release_date
            }
    except:
        return None
    return None

################################################################################
#checks sequence quality based on length and the proportion of unknown residues
def check_sequence_quality(seq_record):
    """
    Check if sequence passes quality filters
    Criteria:
    - Length: 50-2000 residues
    - No excessive X's (unknown residues)
    """
    seq_str = str(seq_record.seq)
    seq_len = len(seq_str)
    
    # Length filter
    if seq_len < 50 or seq_len > 2000:
        return False, "length"
    
    # Check for excessive unknown residues (X)
    x_count = seq_str.count('X')
    if x_count / seq_len > 0.05:  # > 5% unknown
        return False, "unknown_residues"
    return True, "pass"



################################################################################
def filter_pdb_sequences(input_fasta, output_fasta, metadata_csv):
    """
    Main filtering function
    """
    passed = []
    failed = {}
    
    # Track unique PDB IDs for metadata lookup; sequences that fail quality 
    # check are appended to failed dict with reason while passing seqs are
    # appended to passed list of dicts with basic info about them 
    unique_pdbs = set()
    
    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Parse header: >101M_A mol:protein length:154  MYOGLOBIN
            pdb_chain = record.id  # e.g., "101M_A"
            pdb_id = pdb_chain.split('_')[0].upper()  # e.g., "101M"
            chain_id = pdb_chain.split('_')[1] if '_' in pdb_chain else 'A'
            
            # Quality filter 1: Sequence quality
            is_valid, reason = check_sequence_quality(record)
            if not is_valid:
                failed[pdb_chain] = reason
                continue
            
            # Add to unique PDBs for metadata lookup
            unique_pdbs.add(pdb_id)
            
            # Write to output
            SeqIO.write(record, out_f, "fasta")
            passed.append({
                'pdb_id': pdb_id,
                'chain': chain_id,
                'length': len(record.seq),
                'description': record.description
            })
    
    print(f"Passed sequence filters: {len(passed)}")
    print(f"Failed filters: {len(failed)}")
    print()
    
    
    # Now fetch metadata for unique PDB IDs; metadata includes resolution, 
    # release date, experimental method
    print(f"Fetching metadata for {len(unique_pdbs)} unique PDB entries (parallelized)...")
    pdb_metadata = {}
    pdb_list = list(unique_pdbs)
    completed = 0
    
    # Use ThreadPoolExecutor for parallel API requests
    with ThreadPoolExecutor(max_workers=20) as executor:
        # Submit all tasks
        future_to_pdb = {executor.submit(get_pdb_metadata, pdb_id): pdb_id 
                         for pdb_id in pdb_list}
        
        # Process results as they complete
        for future in as_completed(future_to_pdb):
            pdb_id = future_to_pdb[future]
            completed += 1
            
            if completed % 500 == 0:
                print(f"  Progress: {completed}/{len(pdb_list)}")
            
            try:
                metadata = future.result()
                if metadata:
                    pdb_metadata[pdb_id] = metadata
            except Exception as e:
                pass  # Skip failed requests
    
    print(f"Retrieved metadata for {len(pdb_metadata)} entries")
    print()
    
    # Apply resolution filter    
    final_passed = []
    for entry in passed:
        pdb_id = entry['pdb_id']
        
        if pdb_id not in pdb_metadata:
            # No metadata available - exclude to be safe
            continue
        
        metadata = pdb_metadata[pdb_id]
        exp_method = metadata['exp_method']
        resolution = metadata['resolution']
        
        # Resolution filter (only for X-ray)
        if 'X-RAY' in exp_method.upper():
            if resolution is None or resolution > 3.0:
                failed[f"{pdb_id}_{entry['chain']}"] = "resolution"
                continue
        
        # Passed all filters!
        entry['exp_method'] = exp_method
        entry['resolution'] = resolution
        entry['release_date'] = metadata['release_date']
        final_passed.append(entry)
    
    print(f"Final dataset after all filters: {len(final_passed)}")
    print()
    
    # Save metadata to CSV
    import csv
    with open(metadata_csv, 'w', newline='') as csv_f:
        fieldnames = ['pdb_id', 'chain', 'length', 'exp_method', 
                      'resolution', 'release_date', 'description']
        writer = csv.DictWriter(csv_f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(final_passed)
    
    print(f"Saved metadata to: {metadata_csv}")
    
    # Print summary statistics
    print("\nSummary Stats:")
    print(f"Total sequences processed: {len(passed) + len(failed)}")
    print(f"Passed all filters: {len(final_passed)}")
    print(f"Failed filters: {len(failed)}")
    print("\nFailure reasons:")
    from collections import Counter
    failure_counts = Counter(failed.values())
    for reason, count in failure_counts.most_common():
        print(f"  {reason}: {count}")

if __name__ == "__main__":
    input_file = "../data/pdb_clustering/pdb_seqres.txt"
    output_file = "../data/pdb_clustering/pdb_seqres_filtered.fasta"
    metadata_file = "../data/pdb_clustering/pdb_metadata.csv"
    
    filter_pdb_sequences(input_file, output_file, metadata_file)



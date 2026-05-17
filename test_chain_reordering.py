#!/usr/bin/env python3
"""
Quick test script for chain reordering functionality on sub1_P system.
Tests the reorder_trajectory_by_chains function.
"""
import os
import sys
import configparser
import MDAnalysis as mda
from utils import reorder_trajectory_by_chains, resolve_trajectory_file, resolve_topology_file, resolve_reference_pdb_file

def main():
    # Load nucleolin.ini
    config_file = 'inis/nucleolin.ini'
    if not os.path.isfile(config_file):
        print(f"ERROR: Config file not found: {config_file}")
        sys.exit(1)
    
    cp = configparser.ConfigParser()
    cp.read(config_file)
    
    # Extract settings
    input_dir = cp.get('paths', 'input_dir')
    system = 'sub1_P'
    variation = 'pictus'
    rep = 1
    traj_format = cp.get('systems', 'traj_format', fallback='xtc')
    top_format = cp.get('systems', 'top_format', fallback='tpr')
    
    print(f"\n{'='*60}")
    print(f"Test: Chain Reordering for {system}/{variation} rep{rep}")
    print(f"{'='*60}")
    
    # Resolve file paths
    traj_file, _ = resolve_trajectory_file(system, variation, rep, traj_format, base_dir=input_dir)
    top_file, _ = resolve_topology_file(system, variation, top_format, base_dir=input_dir)
    ref_pdb, _ = resolve_reference_pdb_file(system, variation, base_dir=input_dir)
    
    print(f"\nInput files:")
    print(f"  Topology:  {top_file}")
    print(f"  Trajectory: {traj_file}")
    print(f"  Reference PDB: {ref_pdb}")
    
    # Check files exist
    for fname, fpath in [('topology', top_file), ('trajectory', traj_file), ('PDB', ref_pdb)]:
        if not os.path.isfile(fpath):
            print(f"ERROR: {fname} not found: {fpath}")
            sys.exit(1)
    
    # Load original trajectory
    print(f"\nLoading original trajectory...")
    u_orig = mda.Universe(top_file, traj_file)
    print(f"  Frames: {len(u_orig.trajectory)}")
    print(f"  Atoms: {len(u_orig.atoms)}")
    
    # Show original chain info
    traj_chain_attr = getattr(u_orig.atoms, 'chainIDs', None)
    if traj_chain_attr is None:
        traj_chain_attr = getattr(u_orig.atoms, 'chainids', None)
    if traj_chain_attr is not None:
        traj_chains = sorted(set(str(c).strip() for c in traj_chain_attr if str(c).strip()))
        print(f"  Trajectory chains: {traj_chains}")
        for ch in traj_chains:
            chain_atoms = u_orig.select_atoms(f'chainid {ch}')
            chain_resids = chain_atoms.residues.resids
            n_res = len(set(chain_resids))
            print(f"    Chain {ch}: {len(chain_atoms)} atoms, {n_res} residues (resids {min(chain_resids)}-{max(chain_resids)})")
    
    # Show PDB chain info
    print(f"\nLoading reference PDB...")
    u_pdb = mda.Universe(ref_pdb)
    ref_chain_attr = getattr(u_pdb.atoms, 'chainIDs', None)
    if ref_chain_attr is None:
        ref_chain_attr = getattr(u_pdb.atoms, 'chainids', None)
    if ref_chain_attr is not None:
        ref_chains = sorted(set(str(c).strip() for c in ref_chain_attr if str(c).strip()))
        print(f"  PDB chains: {ref_chains}")
        for ch in ref_chains:
            chain_atoms = u_pdb.select_atoms(f'chainid {ch}')
            chain_resids = chain_atoms.residues.resids
            n_res = len(set(chain_resids))
            print(f"    Chain {ch}: {len(chain_atoms)} atoms, {n_res} residues (resids {min(chain_resids)}-{max(chain_resids)})")
    
    # Get chain mapping from config
    import json
    chain_mapping_str = cp.get('rmsf', 'chain_mapping', fallback='none')
    if chain_mapping_str.lower() == 'none':
        print(f"\nERROR: No chain_mapping found in config for this test")
        sys.exit(1)
    
    chain_mapping_raw = json.loads(chain_mapping_str)
    chain_mapping = chain_mapping_raw.get(system)
    if not chain_mapping:
        print(f"ERROR: No chain_mapping for system {system} in config")
        sys.exit(1)
    
    print(f"\nChain mapping (PDB → trajectory): {chain_mapping}")
    
    # Test the reordering function
    print(f"\nReordering trajectory...")
    try:
        reordered_traj_path = reorder_trajectory_by_chains(
            u_orig,
            ref_pdb,
            chain_mapping=chain_mapping,
            verbose=True
        )
        print(f"✓ Reordered trajectory written to: {reordered_traj_path}")
    except Exception as e:
        print(f"ERROR during reordering: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Verify reordered trajectory
    print(f"\nVerifying reordered trajectory...")
    u_reordered = mda.Universe(top_file, reordered_traj_path)
    print(f"  Frames: {len(u_reordered.trajectory)}")
    print(f"  Atoms: {len(u_reordered.atoms)}")
    
    # Check atom ordering matches PDB
    reord_chain_attr = getattr(u_reordered.atoms, 'chainIDs', None)
    if reord_chain_attr is None:
        reord_chain_attr = getattr(u_reordered.atoms, 'chainids', None)
    if reord_chain_attr is not None:
        reord_chains = sorted(set(str(c).strip() for c in reord_chain_attr if str(c).strip()))
        print(f"  Reordered chains: {reord_chains}")
        for ch in reord_chains:
            chain_atoms = u_reordered.select_atoms(f'chainid {ch}')
            chain_resids = chain_atoms.residues.resids
            n_res = len(set(chain_resids))
            print(f"    Chain {ch}: {len(chain_atoms)} atoms, {n_res} residues (resids {min(chain_resids)}-{max(chain_resids)})")
    
    print(f"\n{'='*60}")
    print(f"✓ Test passed! Reordered trajectory is ready.")
    print(f"{'='*60}\n")

if __name__ == '__main__':
    main()

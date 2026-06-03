#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "mdanalysis>=2.8.0",
#     "parmed",
#     "numpy",
# ]
# ///

import os
import sys
import argparse
import re

try:
    import MDAnalysis as mda
    from MDAnalysis import transformations
except ImportError:
    print("ERROR: MDAnalysis is not installed.", file=sys.stderr)
    print("Please run this script using 'uv run' or install it via pip:", file=sys.stderr)
    print("  pip install MDAnalysis", file=sys.stderr)
    sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert VMD/TCL trajectory preparation scripts (unwrap, center, wrap, strip, renumber) to Python using MDAnalysis."
    )
    parser.add_argument(
        "-t", "--topology", required=True,
        help="Input topology file (e.g., .tpr, .gro, .pdb)"
    )
    parser.add_argument(
        "-i", "--trajectory", required=True,
        help="Input trajectory file (e.g., .xtc, .trr, .dcd)"
    )
    parser.add_argument(
        "--out-pdb",
        help="Output reference PDB file path (first frame)"
    )
    parser.add_argument(
        "--out-psf",
        help="Output reference PSF file path (first frame)"
    )
    parser.add_argument(
        "--out-traj",
        help="Output stripped trajectory file path (e.g., .dcd, .xtc)"
    )
    parser.add_argument(
        "--anchor", default="protein",
        help="MDAnalysis selection string for the centering anchor (default: 'protein')"
    )
    parser.add_argument(
        "--strip-sel", default="not water and not (type CLA SOD POT)",
        help="MDAnalysis selection string for atoms to keep (default: 'not water and not (type CLA SOD POT)')"
    )
    parser.add_argument(
        "--guess-bonds", action="store_true",
        help="Force guessing of bonds instead of using topology bonds"
    )
    parser.add_argument(
        "--protein-sel", default="protein",
        help="MDAnalysis selection string to define protein residues for renumbering (default: 'protein')"
    )
    return parser.parse_args()


def renumber_residues_and_set_chains(u, strip_sel, protein_sel_str="protein"):
    """
    Renumber protein residues to start from 1 for each chain, and assign chain IDs.
    Matches TCL behavior:
      For protein atoms with segname (segid) matching 'PRO[A-Z]', extracts the letter
      as chain ID, and renumbers residues for that chain starting from 1.
    """
    if 'chainID' not in u.atoms.topology_attrs:
        u.add_TopologyAttr('chainID')
        
    protein = u.select_atoms(protein_sel_str)
    
    resid_map = {}
    next_resid_by_chain = {}
    renumbered_count = 0
    unique_chains = set()
    
    # We iterate over residues of the selected atoms
    for res in strip_sel.residues:
        # Check if the residue has atoms that are part of the protein
        is_protein = any(atom in protein for atom in res.atoms)
        segid = res.segment.segid
        match = re.match(r"^PRO([A-Z])$", segid)
        
        if is_protein and match:
            chain_id = match.group(1)
            original_resid = res.resid
            
            key = (chain_id, original_resid)
            if key not in resid_map:
                if chain_id not in next_resid_by_chain:
                    next_resid_by_chain[chain_id] = 1
                resid_map[key] = next_resid_by_chain[chain_id]
                next_resid_by_chain[chain_id] += 1
                
            new_resid = resid_map[key]
            res.resid = new_resid
            res.atoms.chainIDs = chain_id
            renumbered_count += 1
            unique_chains.add(chain_id)
            
    if renumbered_count > 0:
        print(f"Successfully renumbered {renumbered_count} residues across {len(unique_chains)} chains: {sorted(list(unique_chains))}")
    else:
        print("WARNING: No protein residues matched the segment name pattern 'PRO[A-Z]'.")
        print("Residues were not renumbered, and chain IDs were not set.")


def generate_conect_records(strip_sel):
    """
    Generate CONECT records for a selection of atoms.
    Returns list of formatted CONECT lines.
    """
    # Create mapping from Atom object to 1-based index in the selection (serial number)
    atom_to_serial = {atom: i + 1 for i, atom in enumerate(strip_sel)}
    
    conect_lines = []
    for atom in strip_sel:
        atom_serial = atom_to_serial[atom]
        neighbor_serials = []
        
        # Safely access bonded_atoms
        try:
            bonded_atoms = atom.bonded_atoms
        except mda.exceptions.NoDataError:
            bonded_atoms = []
            
        for neighbor in bonded_atoms:
            if neighbor in atom_to_serial:
                bonded_serial = atom_to_serial[neighbor]
                # Only write each bond once: when atom_serial < bonded_serial
                if bonded_serial > atom_serial:
                    neighbor_serials.append(bonded_serial)
        
        if not neighbor_serials:
            continue
            
        neighbor_serials.sort()
        # Write up to 4 bonded partners per line
        for k in range(0, len(neighbor_serials), 4):
            chunk = neighbor_serials[k:k+4]
            # Format: CONECT followed by up to 5 right-aligned 5-character integers
            line = f"CONECT{atom_serial:5d}" + "".join(f"{val:5d}" for val in chunk)
            conect_lines.append(line)
            
    return conect_lines


def finalize_pdb_with_conect(pdb_path, conect_lines):
    """
    Read PDB file, strip any existing END line, append CONECT records, and add END line.
    """
    if not os.path.exists(pdb_path):
        return
        
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
        
    new_lines = []
    for line in lines:
        stripped = line.strip()
        if stripped == 'END':
            continue
        new_lines.append(line)
        
    # Append CONECT records
    for line in conect_lines:
        new_lines.append(line + '\n')
        
    # Append END record
    new_lines.append('END\n')
    
    with open(pdb_path, 'w') as f:
        f.writelines(new_lines)


def save_psf_via_parmed(strip_sel, psf_path):
    """
    Convert MDAnalysis AtomGroup to ParmEd structure and save as PSF.
    """
    try:
        import parmed as pmd
    except ImportError:
        print("ERROR: The 'parmed' package is required to write PSF files.", file=sys.stderr)
        print("Please run this script using 'uv run' or install it via pip:", file=sys.stderr)
        print("  pip install parmed", file=sys.stderr)
        sys.exit(1)
        
    print("Converting AtomGroup to ParmEd structure...")
    pmd_struct = strip_sel.convert_to("PARMED")
    
    print(f"Writing PSF to {psf_path}...")
    pmd_struct.save(psf_path, overwrite=True)
    print(f"Successfully wrote {psf_path}")


def main():
    args = parse_args()
    
    if not (args.out_pdb or args.out_psf or args.out_traj):
        print("ERROR: At least one of --out-pdb, --out-psf, or --out-traj must be specified.", file=sys.stderr)
        sys.exit(1)
        
    # Create output directories if they don't exist
    for out_file in [args.out_pdb, args.out_psf, args.out_traj]:
        if out_file:
            out_dir = os.path.dirname(out_file)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
                
    print(f"Loading topology: {args.topology}")
    print(f"Loading trajectory: {args.trajectory}")
    
    # Load Universe
    u = mda.Universe(args.topology, args.trajectory)
    print(f"System contains {len(u.atoms)} atoms and {len(u.trajectory)} frames.")
    
    # Check box dimensions
    if u.atoms.dimensions is None or any(d <= 0 for d in u.atoms.dimensions[:3]):
        print("WARNING: Box dimensions are missing or zero. PBC wrapping/unwrapping might not function correctly.")
    else:
        print(f"Periodic box dimensions: {u.atoms.dimensions}")
        
    # Guess bonds if requested or if no bonds are present in the topology
    has_bonds = False
    try:
        if hasattr(u, 'bonds') and len(u.bonds) > 0:
            has_bonds = True
    except mda.exceptions.NoDataError:
        pass
        
    if args.guess_bonds or not has_bonds:
        print("Guessing bonds for the universe (this might take a few seconds)...")
        u.atoms.guess_bonds()
        print(f"Guessed {len(u.bonds)} bonds.")
    else:
        print(f"Using {len(u.bonds)} bonds from topology file.")
        
    # Select anchor atoms
    anchor = u.select_atoms(args.anchor)
    if len(anchor) == 0:
        print(f"ERROR: Anchor selection '{args.anchor}' matched 0 atoms.", file=sys.stderr)
        sys.exit(1)
    print(f"Anchor group '{args.anchor}' matches {len(anchor)} atoms.")
    
    # Set up transformations
    print(f"Adding PBC unwrap and wrap transformations (anchor: {args.anchor})...")
    workflow = [
        transformations.unwrap(u.atoms),
        transformations.center_in_box(anchor, center='mass'),
        transformations.wrap(u.atoms, compound='fragments')
    ]
    u.trajectory.add_transformations(*workflow)
    
    # Select stripped atoms
    try:
        strip_sel = u.select_atoms(args.strip_sel)
    except Exception as e:
        print(f"ERROR: Selection string '{args.strip_sel}' is invalid: {e}", file=sys.stderr)
        print("Note: If your topology does not contain atom types, you might need to use 'resname' or 'name' instead of 'type'.", file=sys.stderr)
        sys.exit(1)
        
    print(f"Stripped selection '{args.strip_sel}' matches {len(strip_sel)} atoms.")
    if len(strip_sel) == 0:
        print("ERROR: Stripped selection matched 0 atoms.", file=sys.stderr)
        sys.exit(1)
        
    # Perform PDB/PSF writing using frame 0
    if args.out_pdb or args.out_psf:
        print("Setting trajectory to frame 0 for reference structure generation...")
        u.trajectory[0]  # This loads frame 0 and applies transformations
        
        # Apply renumbering and chain assignment in-place
        renumber_residues_and_set_chains(u, strip_sel, protein_sel_str=args.protein_sel)
        
        if args.out_pdb:
            print(f"Writing reference PDB to {args.out_pdb}...")
            strip_sel.write(args.out_pdb)
            
            # Generate and append CONECT records
            conect_lines = generate_conect_records(strip_sel)
            finalize_pdb_with_conect(args.out_pdb, conect_lines)
            print(f"Successfully wrote reference PDB: {args.out_pdb}")
            
        if args.out_psf:
            save_psf_via_parmed(strip_sel, args.out_psf)
            
    # Perform trajectory writing
    if args.out_traj:
        print(f"Writing stripped trajectory to {args.out_traj}...")
        with mda.Writer(args.out_traj, n_atoms=strip_sel.n_atoms) as writer:
            for ts in u.trajectory:
                # Note: transformations are applied on-the-fly as we iterate
                writer.write(strip_sel)
        print(f"Successfully wrote trajectory: {args.out_traj}")


if __name__ == "__main__":
    main()

import os
import pickle
import argparse
import json
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import hydrogenbonds
from MDAnalysis.exceptions import NoDataError
from utils import (
    convert_time_from_ps,
    validate_time_unit,
    transform_trajectory,
    build_complex_selection,
    resolve_trajectory_file,
    resolve_topology_file,
    resolve_reference_pdb_file,
)
from parallelization import ParallelConfig, get_run_kwargs, safe_run


VALID_HBOND_PRESETS = {
    'custom',
    'protein_protein',
    'protein_nucleic',
    'protein_ligand',
}


def _hbond_preset_defaults(preset):
    """Return default selection configuration for a known H-bond preset."""
    if preset == 'protein_protein':
        return {
            'acceptors_sel': 'protein and name O OXT OT1 OT2 OD1 OD2 OE1 OE2 OG OG1 OH ND1 NE2 N',
            'hydrogens_sel': 'protein and name H* HN HT* HXT',
            'between_pairs': [['protein', 'protein']],
        }
    if preset == 'protein_nucleic':
        return {
            'acceptors_sel': '(protein and name O OXT OT1 OT2 OD1 OD2 OE1 OE2 OG OG1 OH ND1 NE2 N) or (nucleic and name O* N*)',
            'hydrogens_sel': '(protein and name H* HN HT* HXT) or (nucleic and name H*)',
            'between_pairs': [['protein', 'nucleic']],
        }
    if preset == 'protein_ligand':
        ligand_sel = (
            'not protein and not nucleic and '
            'not (resname WAT or resname TIP3 or resname SOL or resname HOH) and '
            'not (name NA or name CL or name K or name CA or name MG)'
        )
        return {
            'acceptors_sel': f'(protein and name O OXT OT1 OT2 OD1 OD2 OE1 OE2 OG OG1 OH ND1 NE2 N) or ({ligand_sel} and name O* N*)',
            'hydrogens_sel': f'(protein and name H* HN HT* HXT) or ({ligand_sel} and name H*)',
            'between_pairs': [['protein', ligand_sel]],
        }
    return {
        'acceptors_sel': None,
        'hydrogens_sel': None,
        'between_pairs': None,
    }


def _resolve_hbond_selection_config(hbonds_preset, acceptors_sel, hydrogens_sel, between_pairs):
    """Merge preset defaults with explicit overrides from configuration."""
    preset = (hbonds_preset or 'custom').strip().lower()
    if preset not in VALID_HBOND_PRESETS:
        raise ValueError(
            f"Unknown hbonds_preset '{hbonds_preset}'. "
            f"Accepted values: {', '.join(sorted(VALID_HBOND_PRESETS))}"
        )

    resolved = _hbond_preset_defaults(preset)
    if acceptors_sel is not None:
        resolved['acceptors_sel'] = acceptors_sel
    if hydrogens_sel is not None:
        resolved['hydrogens_sel'] = hydrogens_sel
    if between_pairs is not None:
        resolved['between_pairs'] = between_pairs

    return preset, resolved['acceptors_sel'], resolved['hydrogens_sel'], resolved['between_pairs']


def _is_portable_scalar(value):
    """Return True when *value* is a safe scalar for pickle portability."""
    return isinstance(
        value,
        (str, bytes, bool, int, float, complex, np.generic),
    ) or value is None


def _to_portable_array(name, value):
    """Convert analysis outputs to arrays and reject embedded custom objects."""
    arr = np.array(value, copy=True)
    if arr.dtype == object:
        for item in arr.flat:
            if not _is_portable_scalar(item):
                raise TypeError(
                    f"{name} contains non-portable object of type "
                    f"{type(item).__name__}. Re-run H-bonds analysis with "
                    "portable result extraction."
                )
    return arr


def _build_hbonds_payload(hbonds, *, d_a_cutoff, d_h_a_angle_cutoff,
                          start_frame, update_selections, acceptors_sel,
                          hydrogens_sel, between_pairs, parallel_backend,
                          n_workers, hbonds_preset='custom', time_unit='ps',
                          time_corrected=False, atom_labels_by_index=None):
    """Build a strict portable dict payload for H-bonds plotting."""
    return {
        'times': _to_portable_array('times', hbonds.times),
        'count_by_time': _to_portable_array('count_by_time', hbonds.count_by_time()),
        'count_by_type': _to_portable_array('count_by_type', hbonds.count_by_type()),
        'count_by_ids': _to_portable_array('count_by_ids', hbonds.count_by_ids()),
        'time_unit': str(time_unit),
        'time_corrected': bool(time_corrected),
        'd_a_cutoff': float(d_a_cutoff),
        'd_h_a_angle_cutoff': float(d_h_a_angle_cutoff),
        'start_frame': int(start_frame),
        'update_selections': bool(update_selections),
        'acceptors_sel': acceptors_sel,
        'hydrogens_sel': hydrogens_sel,
        'between_pairs': between_pairs,
        'hbonds_preset': hbonds_preset,
        'parallel_backend': parallel_backend,
        'n_workers': int(n_workers),
        'atom_labels_by_index': atom_labels_by_index or {},
    }


def _build_atom_labels_map(universe, count_by_ids):
    """Return residue-only label mapping for donor/hydrogen/acceptor ids used in count_by_ids."""
    labels = {}
    rows = np.asarray(count_by_ids)
    if rows.size == 0:
        return labels

    # Collect referenced atom identifiers (id/index convention depends on MDAnalysis backend).
    referenced = set()
    for row in rows:
        referenced.add(int(row[0]))
        referenced.add(int(row[1]))
        referenced.add(int(row[2]))

    # Build a robust lookup by both atom.index and atom.id.
    try:
        atom_iter = iter(universe.atoms)
    except Exception:
        return labels

    for atom in atom_iter:
        try:
            label = f"{atom.resname}{atom.resid}"
            idx = int(atom.index)
            aid = int(atom.id)
        except Exception:
            continue
        if idx in referenced and idx not in labels:
            labels[idx] = label
        if aid in referenced and aid not in labels:
            labels[aid] = label

    return labels


def _guess_bonds_if_missing(universe):
    """Try to infer bonds so donor-hydrogen pairs can be built from topology."""
    try:
        universe.atoms.guess_bonds()
        return True
    except Exception:
        return False


def _build_hbond_analysis_with_fallback(
    universe,
    *,
    acceptors_sel,
    hydrogens_sel,
    between_pairs,
    d_a_cutoff,
    d_h_a_angle_cutoff,
    update_selections,
):
    """Create HydrogenBondAnalysis with robust fallbacks for sparse topologies."""
    kwargs = {
        'acceptors_sel': acceptors_sel,
        'hydrogens_sel': hydrogens_sel,
        'd_a_cutoff': d_a_cutoff,
        'd_h_a_angle_cutoff': d_h_a_angle_cutoff,
        'update_selections': update_selections,
    }
    if between_pairs is not None:
        kwargs['between'] = between_pairs

    # First try the default path that relies on topology bonds.
    try:
        return hydrogenbonds.HydrogenBondAnalysis(
            universe,
            donors_sel=None,
            **kwargs,
        )
    except (NoDataError, IndexError) as first_error:
        # Retry after bond guessing for topologies lacking explicit bonds
        # (common for nucleic acid base hydrogens in CHARMM PSFs).
        if _guess_bonds_if_missing(universe):
            try:
                print("INFO: Guessed bonds for H-bond analysis fallback.")
                return hydrogenbonds.HydrogenBondAnalysis(
                    universe,
                    donors_sel=None,
                    **kwargs,
                )
            except (NoDataError, IndexError):
                pass

        # Final fallback: avoid topology donor-hydrogen mapping and let
        # MDAnalysis pair hydrogens to explicit donor candidates by distance.
        if hydrogens_sel:
            donors_sel = f"not ({hydrogens_sel})"
            print("INFO: Falling back to distance-based donor matching for H-bonds.")
            return hydrogenbonds.HydrogenBondAnalysis(
                universe,
                donors_sel=donors_sel,
                **kwargs,
            )

        raise ValueError(
            "HydrogenBondAnalysis could not infer donor-hydrogen pairs from this topology. "
            "Provide explicit --atom-selections (acceptors and hydrogens), optionally together with --between-pairs. "
            "Underlying error: "
            f"{first_error}"
        ) from first_error

# Example selection strings for Hydrogen Bonds analysis:
#
# Option 1 - Atom-focused analysis:
#   acceptors_sel = 'protein and name O*'    # Acceptor atoms
#   hydrogens_sel = 'nucleic and name H*'    # Hydrogen atoms
#
# Option 2 - Pair-focused analysis:
#   between_pairs = [['protein and resnum 73', 'nucleic and name NH*'],
#                    ['protein and resnum 73', 'nucleic and name N*']]
#
# Additional parameters:
#   d_a_cutoff = 3.5         # Distance cutoff between donor and acceptor
#   d_h_a_angle_cutoff = 150 # Angle cutoff for hydrogen bond

def run_hbonds_analysis(systems, variations, num_replicates, d_a_cutoff, d_h_a_angle_cutoff, start_frame, traj_format, top_format='top', acceptors_sel=None, hydrogens_sel=None, between_pairs=None, update_selections=True, wrap_selection='auto', output_dir=None, parallel_backend='serial', n_workers=None, strict_wrapping=False, wrapped_manifest=None, input_dir=None, require_reference_pdb=False, hbonds_preset='custom', time_interval_between_frames=None, time_unit='ps'):
    """
    Runs the Hydrogen Bonds analysis for each system and variation and saves results as individual pickle files.

    Parameters
    ----------
    wrap_selection : str or None, optional
        Controls PBC wrapping.  ``'auto'`` (default) unwraps/centres on
        protein and wraps everything else.  ``None`` disables wrapping.
        Any other string is an MDAnalysis selection for atoms to wrap.
    output_dir : str or None, optional
        Directory for pickle output.  When *None* (default), pickles are
        written to the current working directory.
    parallel_backend : str, optional
        MDAnalysis parallel backend: ``'serial'`` (default),
        ``'multiprocessing'``, or ``'dask'``.
    n_workers : int or None, optional
        Number of parallel workers.  ``None`` auto-detects.
    """
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    _pkl = (lambda name: os.path.join(output_dir, name)) if output_dir else (lambda name: name)

    # Build parallelisation config — H-bonds supports parallel backends.
    parallel_cfg = ParallelConfig(backend=parallel_backend, n_workers=n_workers)
    run_kwargs = get_run_kwargs(parallel_cfg, analysis_key='hbonds')
    validate_time_unit(time_unit)

    print("--- Starting Hydrogen Bonds calculation for each system and variation... ---")

    preset_name, resolved_acceptors_sel, resolved_hydrogens_sel, resolved_between_pairs = _resolve_hbond_selection_config(
        hbonds_preset,
        acceptors_sel,
        hydrogens_sel,
        between_pairs,
    )
    print(f"H-bonds preset: {preset_name}")

    reps = range(1, num_replicates + 1)

    for system in systems:
        for variation in variations[system]:
            for rep in reps:
                pickle_file = _pkl(f'hbonds_plot_{system}_{variation}_rep{rep}.pkl')
                if os.path.exists(pickle_file):
                    print(f"Skipping Hydrogen Bonds analysis for {system}, {variation}, replicate {rep} because the data already exists in {pickle_file}.")
                    continue

                print(f"Processing {system}, {variation}, replicate {rep}.")
                if require_reference_pdb:
                    ref_pdb, _ = resolve_reference_pdb_file(system, variation, base_dir=input_dir)
                    if not os.path.isfile(ref_pdb):
                        raise FileNotFoundError(
                            f"Required reference PDB not found for {system}/{variation}: {ref_pdb}"
                        )

                wrapped_key = (system, variation, rep)
                wrapped_traj_file = None
                if wrapped_manifest:
                    wrapped_traj_file = wrapped_manifest.get(wrapped_key)
                    if wrapped_traj_file and not os.path.isfile(wrapped_traj_file):
                        wrapped_traj_file = None

                if wrapped_traj_file:
                    traj_file = wrapped_traj_file
                else:
                    traj_file, _ = resolve_trajectory_file(
                        system, variation, rep, traj_format, base_dir=input_dir
                    )

                top_file, _ = resolve_topology_file(system, variation, top_format, base_dir=input_dir)
                u = mda.Universe(top_file, traj_file)

                # PBC handling: wrap_selection controls which atoms
                # are wrapped back into the primary image.
                complex_ag, ligand_ag, rest_ag = build_complex_selection(u, wrap_selection=wrap_selection)
                if complex_ag is not None and wrapped_traj_file is None:
                    # Keep wrapping serial and run the H-bond analysis on a
                    # clean trajectory without attached transformations so that
                    # parallel backends remain usable.
                    transform_trajectory(u, complex_ag, rest_ag,
                                         ligand_selection=ligand_ag,
                                         strict_wrapping=strict_wrapping)

                    wrapped_traj_local = f'wrapped_{system}_production_{variation}_rep{rep}.{traj_format}'
                    if not os.path.exists(wrapped_traj_local):
                        print(f"Creating wrapped trajectory (serial pre-processing): {wrapped_traj_local}")
                        with mda.Writer(wrapped_traj_local, n_atoms=u.atoms.n_atoms) as writer:
                            for _ts in u.trajectory:
                                writer.write(u.atoms)

                    del u
                    u = mda.Universe(top_file, wrapped_traj_local)

                if resolved_acceptors_sel is None and resolved_hydrogens_sel is None and resolved_between_pairs is None:
                    raise ValueError("You must provide --between-pairs and/or --atom-selections for Hydrogen Bonds analysis.")

                if (resolved_acceptors_sel is None) != (resolved_hydrogens_sel is None):
                    raise ValueError("--atom-selections requires BOTH acceptors and hydrogens selections.")

                hbonds = None
                if resolved_between_pairs is not None or (resolved_acceptors_sel is not None and resolved_hydrogens_sel is not None):
                    hbonds = _build_hbond_analysis_with_fallback(
                        u,
                        acceptors_sel=resolved_acceptors_sel,
                        hydrogens_sel=resolved_hydrogens_sel,
                        between_pairs=resolved_between_pairs,
                        d_a_cutoff=d_a_cutoff,
                        d_h_a_angle_cutoff=d_h_a_angle_cutoff,
                        update_selections=update_selections,
                    )

                if hbonds is None:
                    raise ValueError("Failed to create HydrogenBondAnalysis object. Check your input parameters.")

                safe_run(hbonds, run_kwargs, start=start_frame, stop=None, step=1)

                count_by_ids = hbonds.count_by_ids()
                atom_labels_by_index = _build_atom_labels_map(u, count_by_ids)

                if traj_format.lower() == 'dcd' and time_interval_between_frames is not None:
                    # Rebuild time axis from frame index + user-provided dt (ps).
                    frame_idx = np.arange(start_frame, start_frame + len(hbonds.times), dtype=float)
                    corrected_times_ps = frame_idx * float(time_interval_between_frames)
                    hbonds.times = convert_time_from_ps(corrected_times_ps, time_unit)
                    payload_time_unit = time_unit
                    payload_time_corrected = True
                else:
                    payload_time_unit = 'ps'
                    payload_time_corrected = False

                hbonds_result = _build_hbonds_payload(
                    hbonds,
                    d_a_cutoff=d_a_cutoff,
                    d_h_a_angle_cutoff=d_h_a_angle_cutoff,
                    start_frame=start_frame,
                    update_selections=update_selections,
                    acceptors_sel=resolved_acceptors_sel,
                    hydrogens_sel=resolved_hydrogens_sel,
                    between_pairs=resolved_between_pairs,
                    parallel_backend=parallel_cfg.backend,
                    n_workers=run_kwargs.get('n_workers', 1),
                    hbonds_preset=preset_name,
                    time_unit=payload_time_unit,
                    time_corrected=payload_time_corrected,
                    atom_labels_by_index=atom_labels_by_index,
                )

                with open(pickle_file, 'wb') as f:
                    pickle.dump(hbonds_result, f)

    print("Finished Hydrogen Bonds calculation for all systems and variations.")


def main():
    """Main function to parse arguments and run Hydrogen Bonds analysis."""
    parser = argparse.ArgumentParser(description='Run Hydrogen Bonds analysis on MD trajectories')

    parser.add_argument('--systems', type=str, required=True,
                        help='JSON string or file path containing list of systems (e.g., \'["6q2b", "5h3r"]\')')
    parser.add_argument('--variations', type=str, required=False, default=None,
                        help='JSON string or file path containing variations dict (e.g., \'{"6q2b": ["wild", "k84r"]}\'). Optional: when not provided, no variation is used in filenames.')
    parser.add_argument('--default-variation', type=str, default='',
                        help='Default variation name when --variations not provided (default: empty string, meaning no variation in filename)')
    parser.add_argument('--num-replicates', type=int, default=3,
                        help='Number of replicates per system and variation (default: 3)')
    parser.add_argument('--traj-format', type=str, default='dcd',
                        help='Trajectory file format (default: dcd)')
    parser.add_argument('--top-format', type=str, default='top',
                        help='Topology file format (default: top)')
    parser.add_argument('--start-frame', type=int, default=500,
                        help='Starting frame for analysis, after equilibration (default: 500)')

    parser.add_argument('--d-a-cutoff', type=float, default=3.5,
                        help='Distance cutoff between donor and acceptor (default: 3.5)')
    parser.add_argument('--d-h-a-angle-cutoff', type=float, default=150.0,
                        help='Angle cutoff for hydrogen bond (default: 150.0)')
    parser.add_argument('--update-selections', action='store_true', default=True,
                        help='Update selections over time (default: True)')
    parser.add_argument('--no-update-selections', dest='update_selections', action='store_false',
                        help='Do not update selections over time for better performance')

    # Selection options. --between-pairs and --atom-selections can be combined.
    parser.add_argument('--atom-selections', nargs=2, metavar=('ACCEPTORS', 'HYDROGENS'),
                        help='Explicit atom selections: acceptors_sel and hydrogens_sel')
    parser.add_argument('--between-pairs', type=str,
                        help='JSON string for pair-focused analysis: list of [selection1, selection2] pairs')

    # PBC wrapping control
    parser.add_argument('--wrap-selection', type=str, default='auto',
                        help='PBC wrap selection: "auto" (wrap non-protein), "none" (disable), '
                             'or a custom MDAnalysis selection string (default: auto)')
    parser.add_argument('--strict-wrapping', action='store_true',
                        help='Fail if fragment-aware wrapping cannot be applied')

    # Parallelization
    parser.add_argument('--parallel-backend', type=str, default='serial',
                        choices=['serial', 'multiprocessing', 'dask'],
                        help='MDAnalysis parallel backend (default: serial)')
    parser.add_argument('--n-workers', type=int, default=None,
                        help='Number of parallel workers (default: auto-detect)')
    parser.add_argument('--time-interval-between-frames', type=float, default=None,
                        help='Frame interval (ps) for DCD time-axis correction')
    parser.add_argument('--time-unit', type=str, default='ps',
                        help='Target time unit for corrected DCD axis (default: ps)')

    args = parser.parse_args()

    # Parse systems and variations from JSON
    try:
        if os.path.isfile(args.systems):
            with open(args.systems, 'r', encoding='utf-8') as f:
                systems = json.load(f)
        else:
            systems = json.loads(args.systems)

        if args.variations:
            if os.path.isfile(args.variations):
                with open(args.variations, 'r', encoding='utf-8') as f:
                    variations = json.load(f)
            else:
                variations = json.loads(args.variations)
        else:
            # Auto-generate empty variation for each system (no variation in filename)
            variations = {system: [args.default_variation] for system in systems}
            print(f"  NOTE: No '--variations' provided. Using empty variation (no variation in filename).")
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"Error parsing JSON: {e}")
        return 1

    # Parse between_pairs if provided
    between_pairs = None
    acceptors_sel = None
    hydrogens_sel = None

    if args.between_pairs:
        try:
            between_pairs = json.loads(args.between_pairs)
        except json.JSONDecodeError as e:
            print(f"Error parsing between_pairs JSON: {e}")
            return 1

    if args.atom_selections:
        acceptors_sel, hydrogens_sel = args.atom_selections

    if between_pairs is None and args.atom_selections is None:
        print("Error: You must provide --between-pairs and/or --atom-selections.")
        return 1

    # Run the analysis
    run_hbonds_analysis(
        systems=systems,
        variations=variations,
        num_replicates=args.num_replicates,
        d_a_cutoff=args.d_a_cutoff,
        d_h_a_angle_cutoff=args.d_h_a_angle_cutoff,
        start_frame=args.start_frame,
        traj_format=args.traj_format,
        top_format=args.top_format,
        acceptors_sel=acceptors_sel,
        hydrogens_sel=hydrogens_sel,
        between_pairs=between_pairs,
        update_selections=args.update_selections,
        wrap_selection=None if args.wrap_selection.lower() == 'none' else args.wrap_selection,
        parallel_backend=args.parallel_backend,
        n_workers=args.n_workers,
        strict_wrapping=args.strict_wrapping,
        time_interval_between_frames=args.time_interval_between_frames,
        time_unit=args.time_unit,
    )

    return 0


if __name__ == '__main__':
    exit(main())

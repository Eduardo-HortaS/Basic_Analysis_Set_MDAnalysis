import os
import pickle
import argparse
import json
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import hydrogenbonds
from utils import transform_trajectory, build_complex_selection
from parallelization import ParallelConfig, get_run_kwargs, safe_run


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
                          n_workers):
    """Build a strict portable dict payload for H-bonds plotting."""
    return {
        'times': _to_portable_array('times', hbonds.times),
        'count_by_time': _to_portable_array('count_by_time', hbonds.count_by_time()),
        'count_by_type': _to_portable_array('count_by_type', hbonds.count_by_type()),
        'count_by_ids': _to_portable_array('count_by_ids', hbonds.count_by_ids()),
        'd_a_cutoff': float(d_a_cutoff),
        'd_h_a_angle_cutoff': float(d_h_a_angle_cutoff),
        'start_frame': int(start_frame),
        'update_selections': bool(update_selections),
        'acceptors_sel': acceptors_sel,
        'hydrogens_sel': hydrogens_sel,
        'between_pairs': between_pairs,
        'parallel_backend': parallel_backend,
        'n_workers': int(n_workers),
    }

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

def run_hbonds_analysis(systems, variations, num_replicates, d_a_cutoff, d_h_a_angle_cutoff, start_frame, traj_format, top_format='top', acceptors_sel=None, hydrogens_sel=None, between_pairs=None, update_selections=True, wrap_selection='auto', output_dir=None, parallel_backend='serial', n_workers=None):
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

    print("--- Starting Hydrogen Bonds calculation for each system and variation... ---")

    reps = range(1, num_replicates + 1)

    for system in systems:
        for variation in variations[system]:
            for rep in reps:
                pickle_file = _pkl(f'hbonds_plot_{system}_{variation}_rep{rep}.pkl')
                if os.path.exists(pickle_file):
                    print(f"Skipping Hydrogen Bonds analysis for {system}, {variation}, replicate {rep} because the data already exists in {pickle_file}.")
                    continue

                print(f"Processing {system}, {variation}, replicate {rep}.")
                traj_file = f'{system}/{variation}/{system}_production_{variation}_rep_{rep}.{traj_format}'
                top_file = f'{system}/{variation}/{system}_system_{variation}.{top_format}'
                u = mda.Universe(top_file, traj_file)

                # PBC handling: wrap_selection controls which atoms
                # are wrapped back into the primary image.
                complex_ag, ligand_ag, rest_ag = build_complex_selection(u, wrap_selection=wrap_selection)
                if complex_ag is not None:
                    transform_trajectory(u, complex_ag, rest_ag,
                                         ligand_selection=ligand_ag)

                if acceptors_sel is None and hydrogens_sel is None and between_pairs is None:
                    raise ValueError("You must provide either acceptors_sel and hydrogens_sel, or between_pairs for Hydrogen Bonds analysis.")

                hbonds = None
                if between_pairs is not None:
                    hbonds = hydrogenbonds.HydrogenBondAnalysis(
                        u,
                        donors_sel=None,
                        acceptors_sel=None,
                        hydrogens_sel=None,
                        between=between_pairs,
                        d_a_cutoff=d_a_cutoff,
                        d_h_a_angle_cutoff=d_h_a_angle_cutoff,
                        update_selections=update_selections
                    )
                elif acceptors_sel is not None and hydrogens_sel is not None:
                    hbonds = hydrogenbonds.HydrogenBondAnalysis(
                        u,
                        donors_sel=None,
                        acceptors_sel=acceptors_sel,
                        hydrogens_sel=hydrogens_sel,
                        d_a_cutoff=d_a_cutoff,
                        d_h_a_angle_cutoff=d_h_a_angle_cutoff,
                        update_selections=update_selections
                    )

                if hbonds is None:
                    raise ValueError("Failed to create HydrogenBondAnalysis object. Check your input parameters.")

                safe_run(hbonds, run_kwargs, start=start_frame, stop=None, step=1)

                hbonds_result = _build_hbonds_payload(
                    hbonds,
                    d_a_cutoff=d_a_cutoff,
                    d_h_a_angle_cutoff=d_h_a_angle_cutoff,
                    start_frame=start_frame,
                    update_selections=update_selections,
                    acceptors_sel=acceptors_sel,
                    hydrogens_sel=hydrogens_sel,
                    between_pairs=between_pairs,
                    parallel_backend=parallel_cfg.backend,
                    n_workers=run_kwargs.get('n_workers', 1),
                )

                with open(pickle_file, 'wb') as f:
                    pickle.dump(hbonds_result, f)

    print("Finished Hydrogen Bonds calculation for all systems and variations.")


def main():
    """Main function to parse arguments and run Hydrogen Bonds analysis."""
    parser = argparse.ArgumentParser(description='Run Hydrogen Bonds analysis on MD trajectories')

    parser.add_argument('--systems', type=str, required=True,
                        help='JSON string or file path containing list of systems (e.g., \'["6q2b", "5h3r"]\')')
    parser.add_argument('--variations', type=str, required=True,
                        help='JSON string or file path containing variations dict (e.g., \'{"6q2b": ["wild", "k84r"]}\')')
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

    # Selection options - either atom-focused or pair-focused
    selection_group = parser.add_mutually_exclusive_group(required=True)
    selection_group.add_argument('--atom-selections', nargs=2, metavar=('ACCEPTORS', 'HYDROGENS'),
                                help='Atom-focused analysis: acceptors_sel and hydrogens_sel')
    selection_group.add_argument('--between-pairs', type=str,
                                help='JSON string for pair-focused analysis: list of [donor, acceptor] pairs')

    # PBC wrapping control
    parser.add_argument('--wrap-selection', type=str, default='auto',
                        help='PBC wrap selection: "auto" (wrap non-protein), "none" (disable), '
                             'or a custom MDAnalysis selection string (default: auto)')

    # Parallelization
    parser.add_argument('--parallel-backend', type=str, default='serial',
                        choices=['serial', 'multiprocessing', 'dask'],
                        help='MDAnalysis parallel backend (default: serial)')
    parser.add_argument('--n-workers', type=int, default=None,
                        help='Number of parallel workers (default: auto-detect)')

    args = parser.parse_args()

    # Parse systems and variations from JSON
    try:
        if os.path.isfile(args.systems):
            with open(args.systems, 'r', encoding='utf-8') as f:
                systems = json.load(f)
        else:
            systems = json.loads(args.systems)

        if os.path.isfile(args.variations):
            with open(args.variations, 'r', encoding='utf-8') as f:
                variations = json.load(f)
        else:
            variations = json.loads(args.variations)
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
    else:
        acceptors_sel, hydrogens_sel = args.atom_selections

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
    )

    return 0


if __name__ == '__main__':
    exit(main())

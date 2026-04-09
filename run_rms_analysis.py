import os
import pickle
import argparse
import json
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms, diffusionmap
from utils import transform_trajectory, align_trajectory, build_complex_selection, convert_time_from_ps, validate_time_unit, SUPPORTED_TIME_UNITS, resolve_trajectory_file
from parallelization import ParallelConfig, get_run_kwargs, safe_run

# Example selection strings for different RMS* analyses:

# RMSD Analysis - Compare structural deviation over time:
#   target_selection = 'nucleic'           # Atoms to calculate RMSD for
#   ref_selection = 'protein and backbone' # Atoms used for alignment
#
# RMSF Analysis - Calculate per-residue fluctuations:
#   target_selection = 'protein and backbone' # Atoms to calculate RMSF for
#   ref_selection = 'protein and backbone'    # Atoms used for alignment
#
# 2D-RMSD Analysis - Calculate pairwise frame distances:
#   target_selection = 'protein and backbone' # Atoms to calculate distance matrix for


def validate_chain_intervals(universe, chain_intervals, target_selection):
    """
    Validates that chain_intervals resids (1-based, PDB-style) exist in the trajectory
    and are consistent with the target selection.

    Parameters
    ----------
    universe : mda.Universe
        The loaded MDAnalysis Universe.
    chain_intervals : dict
        Mapping of chain ID to [start_resid, end_resid] (1-based inclusive).
        Example: {"A": [1, 120], "B": [121, 239]}
    target_selection : str
        MDAnalysis selection string for the target atoms.

    Raises
    ------
    ValueError
        If any interval's resids are not found in the trajectory or are invalid.
    """
    target_atoms = universe.select_atoms(target_selection)
    available_resids = set(target_atoms.residues.resids)

    all_chain_resids = []
    for chain_id, (start_resid, end_resid) in chain_intervals.items():
        if start_resid < 1:
            raise ValueError(
                f"Chain '{chain_id}': start_resid={start_resid} is less than 1. "
                f"Chain intervals must use 1-based (PDB-style) residue IDs."
            )
        if start_resid > end_resid:
            raise ValueError(
                f"Chain '{chain_id}': start_resid={start_resid} > end_resid={end_resid}."
            )

        chain_resids = set(range(start_resid, end_resid + 1))
        missing = chain_resids - available_resids
        if missing:
            missing_sorted = sorted(missing)
            raise ValueError(
                f"Chain '{chain_id}': resids {missing_sorted[:10]}{'...' if len(missing_sorted) > 10 else ''} "
                f"not found in target selection '{target_selection}'. "
                f"Available resid range: {min(available_resids)}-{max(available_resids)}."
            )

        all_chain_resids.append((chain_id, chain_resids))

    # Check for overlaps
    for i, (id_a, resids_a) in enumerate(all_chain_resids):
        for j, (id_b, resids_b) in enumerate(all_chain_resids):
            if i >= j:
                continue
            overlap = resids_a & resids_b
            if overlap:
                print(f"WARNING: Chain '{id_a}' and chain '{id_b}' have overlapping resids: {sorted(overlap)[:10]}")

    # Check for gaps
    covered = set()
    for _, resids in all_chain_resids:
        covered |= resids
    uncovered = available_resids - covered
    if uncovered:
        print(f"WARNING: Resids {sorted(uncovered)[:10]}{'...' if len(uncovered) > 10 else ''} "
              f"from target selection are not covered by any chain interval.")


def run_rms_analysis(systems, variations, num_replicates, analysis, target_selection, ref_selection,
                     start_frame, traj_format, top_format='top', group_selections=None, chain_intervals=None,
                     time_interval_between_frames=None, time_unit='ns', ref_suffix='', wrap_selection='auto',
                     output_dir=None, parallel_backend='serial', n_workers=None, strict_wrapping=False):
    """
    After external transforming and aligning of trajectory data for each analysis, system, variation and replicate,
    runs the respective analysis (RMSD (2D and conventional) | RMSF) and saves results as individual pickle files.

    Parameters
    ----------
    systems : list
        List of system names.
    variations : dict
        Mapping of system name to list of variations.
    num_replicates : int
        Number of replicates per system/variation.
    analysis : str
        Type of analysis: 'RMSD', 'RMSF', or '2D-RMSD'.
    target_selection : str
        MDAnalysis selection string for target atoms.
    ref_selection : str
        MDAnalysis selection string for reference atoms (alignment).
    start_frame : int
        Starting frame for analysis (skip equilibration).
    traj_format : str
        Trajectory file format (e.g., 'dcd', 'xtc').
    group_selections : list of str, optional
        Additional group selections for RMSD analysis.
    chain_intervals : dict, optional
        For RMSF chain split. Supports two formats:
        - Per-system: maps system name to chain dict.
          Example: {"Ung_G-C_4": {"A": [1, 239]}, "Mug_G-C_4": {"A": [1, 211]}}
        - Legacy flat: maps chain ID to [start_resid, end_resid] (1-based).
          Example: {"A": [1, 120], "B": [121, 239]}
          Applied identically to every system.
    time_interval_between_frames : float, optional
        Time interval between frames in picoseconds. Used for DCD time-axis correction
        in RMSD analysis. When provided and traj_format is 'dcd', recalculates the
        time column of RMSD results.
    time_unit : str, optional
        Output time unit for the corrected time axis.  Accepted values:
        'fs', 'ps', 'ns', 'us', 'ms', 's'.  Default: 'ns'.
    ref_suffix : str, optional
        Suffix to append to pickle filenames for distinguishing different
        ref_selection iterations (e.g., '_ref0', '_ref1').  Default: ''.
    wrap_selection : str or None, optional
        Controls PBC wrapping.  ``'auto'`` (default) unwraps/centres on
        protein and wraps everything else.  ``None`` disables wrapping.
        Any other string is an MDAnalysis selection for atoms to wrap.
    output_dir : str or None, optional
        Directory for pickle output.  When *None*, pickles are written
        to the current working directory.
    parallel_backend : str, optional
        MDAnalysis parallel backend: ``'serial'`` (default),
        ``'multiprocessing'``, or ``'dask'``.  Only RMSD supports
        non-serial backends; RMSF and 2D-RMSD always run serially.
    n_workers : int or None, optional
        Number of parallel workers.  ``None`` auto-detects.
    """
    validate_time_unit(time_unit)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    _pkl = (lambda name: os.path.join(output_dir, name)) if output_dir else (lambda name: name)

    # Build parallelisation config — only RMSD supports parallel backends.
    # RMSF and 2D-RMSD always run serially, so we only create the config
    # when the analysis actually uses it.
    run_kwargs: dict = {}
    if analysis == 'RMSD':
        parallel_cfg = ParallelConfig(backend=parallel_backend, n_workers=n_workers)
        run_kwargs = get_run_kwargs(parallel_cfg, analysis_key='RMSD')

    print(f"--- Starting {analysis} calculation for each system, condition, and replicate... ---")

    reps = range(1, num_replicates + 1)

    analysis_file_prefix = {
        'RMSD': 'rmsd_plot',
        'RMSF': 'rmsf_plot',
        '2D-RMSD': '2d_rmsd_plot'
    }

    for system in systems:
        for variation in variations[system]:
            for rep in reps:
                # Resolve per-system chain intervals
                system_chain_intervals = None
                if chain_intervals is not None:
                    system_chain_intervals = chain_intervals.get(system)
                    if system_chain_intervals is None:
                        # Legacy flat format: the dict itself is the chain map
                        first_val = next(iter(chain_intervals.values()))
                        if isinstance(first_val, list):
                            system_chain_intervals = chain_intervals

                # Determine if we should skip based on existing files
                # (RMSD handles its own skip-if-exists per selection inside the RMSD block)
                if analysis == 'RMSD':
                    pass  # skip logic is per-selection inside the RMSD block
                elif analysis == 'RMSF' and system_chain_intervals is not None:
                    # Check if all chain pickles exist
                    chain_pickles_exist = all(
                        os.path.exists(_pkl(f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}_chain{chain_id}.pkl'))
                        for chain_id in system_chain_intervals
                    )
                    if chain_pickles_exist:
                        print(f"Skipping {analysis} for {system}, {variation}, replicate {rep} because all chain data files already exist.")
                        continue
                else:
                    pickle_file = _pkl(f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}.pkl')
                    if os.path.exists(pickle_file):
                        print(f"Skipping {analysis} for {system}, {variation}, replicate {rep} because the data already exists in {pickle_file}.")
                        continue

                print(f"Processing {system}, {variation}, replicate {rep}.")
                traj_file, _ = resolve_trajectory_file(
                    system, variation, rep, traj_format
                )
                aligned_traj_file = f'rmsfit_{system}_production_{variation}_reduced_rep{rep}.{traj_format}'
                top_file = f'{system}/{variation}/{system}_system_{variation}.{top_format}'

                if analysis == 'RMSD':
                    u = mda.Universe(top_file, traj_file)

                    # Unwrap, center, and wrap the trajectory to handle
                    # periodic boundary artifacts.  wrap_selection
                    # controls which atoms are wrapped back into the
                    # primary image (default: everything non-protein).
                    complex_ag, ligand_ag, rest_ag = build_complex_selection(
                        u, wrap_selection=wrap_selection)
                    if complex_ag is not None:
                        # Keep wrapping as a serial pre-processing step and
                        # run RMSD on a clean trajectory (no attached
                        # transformations), so MDAnalysis parallel backends can
                        # be used safely.
                        transform_trajectory(u, complex_ag, rest_ag,
                                             ligand_selection=ligand_ag)

                        wrapped_traj_file = f'wrapped_{system}_production_{variation}_rep{rep}.{traj_format}'
                        if not os.path.exists(wrapped_traj_file):
                            print(f"Creating wrapped trajectory (serial pre-processing): {wrapped_traj_file}")
                            with mda.Writer(wrapped_traj_file, n_atoms=u.atoms.n_atoms) as writer:
                                for _ts in u.trajectory:
                                    writer.write(u.atoms)

                        del u
                        u = mda.Universe(top_file, wrapped_traj_file)

                    # Build selection list with target_selection always first.
                    # This guarantees the canonical base pickle
                    # rmsd_plot_*_repX.pkl is always produced for Nextflow,
                    # while optional group selections are emitted as extras.
                    rmsd_selections = [target_selection]
                    if group_selections:
                        for gs in group_selections:
                            if gs != target_selection:
                                rmsd_selections.append(gs)

                    for sel_idx, sel_string in enumerate(rmsd_selections):
                        # Determine pickle filename.
                        # sel_idx == 0 corresponds to target_selection and writes
                        # the canonical base file expected by pipeline processes.
                        if sel_idx == 0:
                            pkl_name = _pkl(f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}{ref_suffix}.pkl')
                        else:
                            pkl_name = _pkl(f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}{ref_suffix}_sel{sel_idx - 1}.pkl')

                        if os.path.exists(pkl_name):
                            print(f"Skipping RMSD (selection {sel_idx}: '{sel_string}') for {system}, {variation}, replicate {rep} — pickle already exists.")
                            continue

                        # Validate that the selection exists in the universe
                        try:
                            sel_atoms = u.select_atoms(sel_string)
                        except Exception as e:
                            print(f"WARNING: Selection '{sel_string}' raised an error for {system}/{variation} rep {rep}: {e} — skipping this selection.")
                            continue
                        if len(sel_atoms) == 0:
                            print(f"WARNING: Selection '{sel_string}' matched 0 atoms in {system}/{variation} rep {rep} — skipping this selection.")
                            continue

                        # Use ref_selection for structural alignment.
                        # When ref_selection differs from sel_string, align on
                        # ref_selection and compute RMSD on sel_string via
                        # groupselections.
                        if ref_selection != sel_string:
                            print(f"  Running RMSD with alignment='{ref_selection}', target='{sel_string}' (sel{sel_idx})...")
                            to_run_rmsd = rms.RMSD(u, u, select=ref_selection,
                                                   groupselections=[sel_string],
                                                   ref_frame=0)
                        else:
                            print(f"  Running RMSD with selection '{sel_string}' (sel{sel_idx})...")
                            to_run_rmsd = rms.RMSD(u, u, select=sel_string, ref_frame=0)
                        safe_run(to_run_rmsd, run_kwargs, start=start_frame, stop=None, step=1)

                        # When groupselections was used, the target RMSD is in
                        # column 3 (col 2 = alignment RMSD).  Normalize so
                        # column 2 always holds the target RMSD for downstream
                        # plotting compatibility.
                        if ref_selection != sel_string:
                            rmsd_raw = to_run_rmsd.results.rmsd
                            # Keep [frame, time, target_rmsd] — drop alignment col
                            to_run_rmsd.results.rmsd = np.column_stack([
                                rmsd_raw[:, 0], rmsd_raw[:, 1], rmsd_raw[:, 3]
                            ])

                        # DCD time-axis correction: recalculate time column when DCD format is used
                        if traj_format.lower() == 'dcd' and time_interval_between_frames is not None:
                            print(f"Applying DCD time-axis correction (dt={time_interval_between_frames} ps, output unit={time_unit})...")
                            rmsd_data = to_run_rmsd.results.rmsd
                            # Column 0 = frame, Column 1 = time, Column 2+ = RMSD values
                            frame_numbers = rmsd_data[:, 0]
                            corrected_time = frame_numbers * time_interval_between_frames  # time in ps
                            corrected_time = convert_time_from_ps(corrected_time, time_unit)
                            rmsd_data[:, 1] = corrected_time

                        # Store portable data (arrays + metadata) to avoid
                        # pickling MDAnalysis runtime objects tied to local
                        # trajectory paths.
                        rmsd_result = {
                            'rmsd_data': np.asarray(to_run_rmsd.results.rmsd),
                            'time_corrected': time_interval_between_frames is not None and traj_format.lower() == 'dcd',
                            'time_unit': time_unit if (time_interval_between_frames is not None and traj_format.lower() == 'dcd') else 'ps',
                            'time_interval_between_frames': time_interval_between_frames,
                            'selection': sel_string,
                            'ref_selection': ref_selection,
                        }

                        with open(pkl_name, 'wb') as f:
                            pickle.dump(rmsd_result, f)
                        print(f"  Saved RMSD pickle: {pkl_name}")

                    # Create aligned trajectory for reuse by subsequent
                    # RMSF / 2D-RMSD runs.  This avoids re-loading the
                    # raw trajectory and re-computing PBC corrections.
                    if not os.path.exists(aligned_traj_file):
                        print(f"Creating aligned trajectory: {aligned_traj_file}")
                        # Reset to frame 0 so AlignTraj uses a deterministic
                        # reference regardless of whether RMSD ran serially
                        # (leaves trajectory at last frame) or in parallel
                        # (leaves trajectory at frame 0).
                        u.trajectory[0]
                        align_trajectory(u, ref_selection, '2D-RMSD', system,
                                         variation, rep, traj_format, start_frame)

                    # Skip the rest of this rep iteration for RMSD
                    continue

                else:  # For RMSF and 2D-RMSD, we need to check if the aligned trajectory already exists.
                    if os.path.exists(aligned_traj_file):
                        print(f"Using pre-aligned trajectory: {aligned_traj_file}")
                        u = mda.Universe(top_file, aligned_traj_file)
                    else:
                        print(f"Aligned trajectory file {aligned_traj_file} does not exist. Creating it now.")
                        u = mda.Universe(top_file, traj_file)

                        # Apply PBC corrections *before* alignment so that
                        # periodic-image artifacts are not baked into the
                        # aligned trajectory.
                        pre_complex_ag, pre_ligand_ag, pre_rest_ag = build_complex_selection(
                            u, wrap_selection=wrap_selection)
                        if pre_complex_ag is not None:
                            transform_trajectory(u, pre_complex_ag, pre_rest_ag,
                                                 ligand_selection=pre_ligand_ag,
                                                 strict_wrapping=strict_wrapping)

                        align_trajectory(u, target_selection, analysis, system, variation, rep, traj_format, start_frame)
                        del u
                        u = mda.Universe(top_file, aligned_traj_file)

                    # PBC handling: wrap_selection controls which atoms
                    # are wrapped back into the primary image.
                    complex_ag, ligand_ag, rest_ag = build_complex_selection(
                        u, wrap_selection=wrap_selection)
                    if complex_ag is not None:
                        transform_trajectory(u, complex_ag, rest_ag,
                                             ligand_selection=ligand_ag,
                                             strict_wrapping=strict_wrapping)

                    if analysis == 'RMSF':
                        print("Calculating RMSF...")

                        # Build selection list with target_selection always first.
                        # This guarantees the canonical base pickle is always produced,
                        # while optional group selections are emitted as extras.
                        rmsf_selections = [target_selection]
                        if group_selections:
                            for gs in group_selections:
                                if gs != target_selection:
                                    rmsf_selections.append(gs)

                        for sel_idx, sel_string in enumerate(rmsf_selections):
                            # Determine pickle filename.
                            # sel_idx == 0 corresponds to target_selection and writes
                            # the canonical base file expected by pipeline processes.
                            if sel_idx == 0:
                                pkl_name_template = f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}.pkl'
                            else:
                                pkl_name_template = f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}_sel{sel_idx - 1}.pkl'

                            # Validate that the selection exists in the universe
                            try:
                                sel_atoms = u.select_atoms(sel_string)
                            except Exception as e:
                                print(f"WARNING: Selection '{sel_string}' raised an error for {system}/{variation} rep {rep}: {e} — skipping this selection.")
                                continue
                            if len(sel_atoms) == 0:
                                print(f"WARNING: Selection '{sel_string}' matched 0 atoms in {system}/{variation} rep {rep} — skipping this selection.")
                                continue

                            # Compute RMSF for this selection
                            rmsf_atoms = sel_atoms
                            to_run_rmsf = rms.RMSF(rmsf_atoms).run()

                            # rms.RMSF computes per-ATOM fluctuations.  Convert to
                            # per-RESIDUE by averaging atom RMSF values within each
                            # residue — this is the standard for publication plots.
                            rmsf_values_atoms = to_run_rmsf.results.rmsf
                            atom_resids = np.array([atom.resid for atom in rmsf_atoms])

                            def _per_residue_rmsf(resid_list):
                                """Average atom-level RMSF for each residue in *resid_list*."""
                                rmsf_list, kept_resids = [], []
                                for resid in resid_list:
                                    mask = atom_resids == resid
                                    vals = rmsf_values_atoms[mask]
                                    if vals.size == 0:
                                        continue
                                    rmsf_list.append(np.mean(vals))
                                    kept_resids.append(resid)
                                return np.array(rmsf_list), np.array(kept_resids)

                            # Apply chain-interval splitting only to the canonical
                            # base selection (sel_idx == 0). For explicit group
                            # selections (e.g., chain-specific selections), writing
                            # per-selection RMSF directly avoids interval/selection
                            # mismatches.
                            if system_chain_intervals is not None and sel_idx == 0:
                                # Validate chain intervals against the trajectory
                                validate_chain_intervals(u, system_chain_intervals, sel_string)

                                # Split RMSF results by chain
                                for chain_id, (start_resid, end_resid) in system_chain_intervals.items():
                                    present_resids = [r for r in range(start_resid, end_resid + 1)
                                                      if r in set(atom_resids)]
                                    chain_rmsf, chain_original_resids = _per_residue_rmsf(present_resids)
                                    chain_renumbered_resids = np.arange(1, len(chain_rmsf) + 1)

                                    chain_result = {
                                        'rmsf': chain_rmsf,
                                        'resids': chain_renumbered_resids,
                                        'chain_id': chain_id,
                                        'original_resids': chain_original_resids,
                                        'selection': sel_string,
                                    }

                                    if sel_idx == 0:
                                        chain_pkl = _pkl(f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}_chain{chain_id}.pkl')
                                    else:
                                        chain_pkl = _pkl(f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}_sel{sel_idx - 1}_chain{chain_id}.pkl')

                                    with open(chain_pkl, 'wb') as f:
                                        pickle.dump(chain_result, f)
                                    print(f"  Saved chain {chain_id} RMSF ({len(chain_rmsf)} residues, sel_idx={sel_idx}) to {chain_pkl}")
                            else:
                                # No chain split — still convert to per-residue RMSF
                                unique_resids = sorted(set(atom_resids))
                                full_rmsf, full_resids = _per_residue_rmsf(unique_resids)

                                rmsf_result = {
                                    'rmsf': full_rmsf,
                                    'resids': full_resids,
                                    'chain_id': '',
                                    'original_resids': full_resids,
                                    'selection': sel_string,
                                }

                                pkl_path = _pkl(pkl_name_template)
                                with open(pkl_path, 'wb') as f:
                                    pickle.dump(rmsf_result, f)
                                print(f"  Saved RMSF ({len(full_rmsf)} residues, sel_idx={sel_idx}) to {pkl_path}")
                    elif analysis == '2D-RMSD':
                        print("Calculating 2D-RMSD...")
                        matrix_2d_rmsd = diffusionmap.DistanceMatrix(u, select=target_selection).run()

                        matrix_result = {
                            'dist_matrix': np.asarray(matrix_2d_rmsd.results.dist_matrix),
                            'selection': target_selection,
                        }

                        with open(_pkl(f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}.pkl'), 'wb') as f:
                            pickle.dump(matrix_result, f)

    print(f"Finished {analysis} calculation for all systems, conditions, and replicates.")


def main():
    """Main function to parse arguments and run RMS analysis."""
    parser = argparse.ArgumentParser(description='Run RMS analysis (RMSD, RMSF, or 2D-RMSD) on MD trajectories')

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

    parser.add_argument('--analysis', type=str, required=True, choices=['RMSD', 'RMSF', '2D-RMSD'],
                        help='Type of RMS analysis to perform')
    parser.add_argument('--target-selection', type=str, required=True,
                        help='MDAnalysis selection string for target atoms')
    parser.add_argument('--ref-selection', type=str, required=True,
                        help='MDAnalysis selection string for reference atoms (alignment)')
    parser.add_argument('--group-selections', type=str, nargs='*',
                        help='Additional group selections for RMSD analysis')

    # Chain intervals for RMSF chain split
    parser.add_argument('--chain-intervals', type=str, default=None,
                        help='JSON string or file path with chain intervals dict for RMSF chain split. '
                             'Maps chain ID to [start_resid, end_resid] (1-based PDB-style). '
                             'Example: \'{"A": [1, 120], "B": [121, 239]}\'')

    # DCD time-axis correction parameters
    parser.add_argument('--time-interval-between-frames', type=float, default=None,
                        help='Time interval between frames in picoseconds (for DCD time-axis correction)')
    parser.add_argument('--time-unit', type=str, default='ns', choices=list(SUPPORTED_TIME_UNITS),
                        help='Output time unit for corrected time axis (default: ns)')

    # PBC wrapping control
    parser.add_argument('--wrap-selection', type=str, default='auto',
                        help='PBC wrap selection: "auto" (wrap non-protein), "none" (disable), '
                             'or a custom MDAnalysis selection string (default: auto)')
    parser.add_argument('--strict-wrapping', action='store_true',
                        help='Fail if fragment-aware wrapping cannot be applied')

    # Parallelization
    parser.add_argument('--parallel-backend', type=str, default='serial',
                        choices=['serial', 'multiprocessing', 'dask'],
                        help='MDAnalysis parallel backend (default: serial). '
                             'Only RMSD supports non-serial backends.')
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

    # Parse chain_intervals if provided
    chain_intervals = None
    if args.chain_intervals is not None:
        try:
            if os.path.isfile(args.chain_intervals):
                with open(args.chain_intervals, 'r', encoding='utf-8') as f:
                    chain_intervals = json.load(f)
            else:
                chain_intervals = json.loads(args.chain_intervals)
        except (json.JSONDecodeError, FileNotFoundError) as e:
            print(f"Error parsing chain_intervals JSON: {e}")
            return 1

    # Run the analysis
    run_rms_analysis(
        systems=systems,
        variations=variations,
        num_replicates=args.num_replicates,
        analysis=args.analysis,
        target_selection=args.target_selection,
        ref_selection=args.ref_selection,
        start_frame=args.start_frame,
        traj_format=args.traj_format,
        top_format=args.top_format,
        group_selections=args.group_selections,
        chain_intervals=chain_intervals,
        time_interval_between_frames=args.time_interval_between_frames,
        time_unit=args.time_unit,
        wrap_selection=None if args.wrap_selection.lower() == 'none' else args.wrap_selection,
        parallel_backend=args.parallel_backend,
        n_workers=args.n_workers,
        strict_wrapping=args.strict_wrapping,
    )

    return 0


if __name__ == '__main__':
    exit(main())

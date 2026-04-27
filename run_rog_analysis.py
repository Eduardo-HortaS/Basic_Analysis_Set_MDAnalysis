import os
import sys
import pickle
import argparse
import json
import MDAnalysis as mda
import numpy as np
from utils import (
    transform_trajectory,
    build_complex_selection,
    convert_time_from_ps,
    validate_time_unit,
    SUPPORTED_TIME_UNITS,
    resolve_trajectory_file,
    resolve_topology_file,
    resolve_reference_pdb_file,
)

# Make the module importable under its own name even when run as ``__main__``.
# This is required so that ``pickle`` can resolve ``RoGResults`` when the script
# is executed directly (e.g. by Nextflow).
sys.modules.setdefault('run_rog_analysis', sys.modules[__name__])


class RoGResults:
    """Container for Radius of Gyration results, compatible with plotting functions."""

    def __init__(self, frames, times, rog_values):
        self.rog_data = np.column_stack((frames, times, rog_values))


# Force the canonical module name so that pickle always stores
# ``run_rog_analysis.RoGResults`` regardless of whether this file is
# executed as ``__main__`` (Nextflow) or imported as a library (executor.py).
RoGResults.__module__ = 'run_rog_analysis'


# Example selection strings for Radius of Gyration analysis:
#
# RoG Analysis - Calculate radius of gyration over time:
#   selection = 'protein and backbone' # Atoms to calculate RoG for (default)
#   selection = 'protein'              # All protein atoms
#   selection = 'nucleic'              # All nucleic acid atoms

def run_rog_analysis(systems, variations, num_replicates, start_frame, traj_format, top_format='top', selection='protein and backbone', time_unit='ns', wrap_selection='auto', output_dir=None, strict_wrapping=False, wrapped_manifest=None, input_dir=None, require_reference_pdb=False):
    """
    Runs the Radius of Gyration analysis for each system and variation and saves results as individual pickle files.

    RoG is always serial — it uses a manual per-frame loop rather than an
    ``AnalysisBase`` subclass, so MDAnalysis parallel backends do not apply.

    Parameters
    ----------
    time_unit : str, optional
        Output time unit: 'fs', 'ps', 'ns', 'us', 'ms', or 's'.  Default: 'ns'.
    wrap_selection : str or None, optional
        Controls PBC wrapping.  ``'auto'`` (default) unwraps/centres on
        protein and wraps everything else.  ``None`` disables wrapping.
        Any other string is an MDAnalysis selection for atoms to wrap.
    output_dir : str or None, optional
        Directory for pickle output.  When *None* (default), pickles are
        written to the current working directory.
    """
    validate_time_unit(time_unit)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    _pkl = (lambda name: os.path.join(output_dir, name)) if output_dir else (lambda name: name)
    print("--- Starting Radius of Gyration calculation for each system and variation... ---")
    analysis = 'RoG'

    reps = range(1, num_replicates + 1)
    analysis_file_prefix = 'rog_plot'

    for system in systems:
        for variation in variations[system]:
            for rep in reps:
                pickle_file = _pkl(f'{analysis_file_prefix}_{system}_{variation}_rep{rep}.pkl')
                if os.path.exists(pickle_file):
                    print(f"Skipping {analysis} for {system}, {variation}, replicate {rep} because the data already exists in {pickle_file}.")
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
                        system, variation, rep, traj_format
                    )

                top_file, _ = resolve_topology_file(system, variation, top_format)
                u = mda.Universe(top_file, traj_file)

                # PBC handling: wrap_selection controls which atoms
                # are wrapped back into the primary image.
                complex_ag, ligand_ag, rest_ag = build_complex_selection(u, wrap_selection=wrap_selection)
                if complex_ag is not None and wrapped_traj_file is None:
                    transform_trajectory(u, complex_ag, rest_ag,
                                         ligand_selection=ligand_ag,
                                         strict_wrapping=strict_wrapping)

                # Select atoms for RoG calculation
                selected_atoms = u.select_atoms(selection)

                print("Calculating Radius of Gyration...")

                # Store results
                rog_results = []
                frames = []
                times = []

                # Iterate through trajectory starting from start_frame
                for ts in u.trajectory[start_frame:]:
                    rog_value = selected_atoms.radius_of_gyration()
                    rog_results.append(rog_value)
                    frames.append(ts.frame)
                    times.append(convert_time_from_ps(ts.time, time_unit))

                rog_analysis_results = RoGResults(frames, times, rog_results)

                rog_result = {
                    'rog_obj': rog_analysis_results,
                    'time_unit': time_unit,
                    'selection': selection,
                }

                with open(pickle_file, 'wb') as f:
                    pickle.dump(rog_result, f)

    print("Finished Radius of Gyration calculation for all systems and variations.")


def main():
    """Main function to parse arguments and run RoG analysis."""
    parser = argparse.ArgumentParser(description='Run Radius of Gyration (RoG) analysis on MD trajectories')

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

    parser.add_argument('--selection', type=str, default='protein and backbone',
                        help='MDAnalysis selection string for atoms to calculate RoG (default: "protein and backbone")')
    parser.add_argument('--time-unit', type=str, default='ns', choices=list(SUPPORTED_TIME_UNITS),
                        help='Output time unit (default: ns)')

    # PBC wrapping control
    parser.add_argument('--wrap-selection', type=str, default='auto',
                        help='PBC wrap selection: "auto" (wrap non-protein), "none" (disable), '
                             'or a custom MDAnalysis selection string (default: auto)')
    parser.add_argument('--strict-wrapping', action='store_true',
                        help='Fail if fragment-aware wrapping cannot be applied')

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

    # Run the analysis
    run_rog_analysis(
        systems=systems,
        variations=variations,
        num_replicates=args.num_replicates,
        start_frame=args.start_frame,
        traj_format=args.traj_format,
        top_format=args.top_format,
        selection=args.selection,
        time_unit=args.time_unit,
        wrap_selection=None if args.wrap_selection.lower() == 'none' else args.wrap_selection,
        strict_wrapping=args.strict_wrapping,
    )

    return 0


if __name__ == '__main__':
    exit(main())

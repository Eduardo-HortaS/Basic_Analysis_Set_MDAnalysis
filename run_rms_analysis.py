import os
import pickle
import argparse
import json
import MDAnalysis as mda
from MDAnalysis.analysis import rms, diffusionmap
from utils import transform_trajectory, align_trajectory

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


def run_rms_analysis(systems, variations, num_replicates, analysis, target_selection, ref_selection, start_frame, traj_format, group_selections=None):
    """
    After external transforming and aligning of trajectory data for each analysis, system, variation and replicate,
    runs the respective analysis (RMSD (2D and conventional) | RMSF) and saves results as individual pickle files.
    """
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
                pickle_file = f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}.pkl'
                if os.path.exists(pickle_file):
                    print(f"Skipping {analysis} for {system}, {variation}, replicate {rep} because the data already exists in {pickle_file}.")
                    continue

                print(f"Processing {system}, {variation}, replicate {rep}.")
                traj_file = f'{system}/{variation}/{system}_production_{variation}_rep_{rep}.{traj_format}'
                aligned_traj_file = f'rmsfit_{system}_production_{variation}_reduced_rep{rep}.{traj_format}'
                top_file = f'{system}/{variation}/{system}_system_{variation}.top'

                if analysis == 'RMSD':
                    u = mda.Universe(top_file, traj_file)
                    to_run_rmsd = rms.RMSD(u, u, select=target_selection, groupselections=group_selections, ref_frame=0)
                    to_run_rmsd.run(start=start_frame, stop=None, step=1)

                    with open(f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}.pkl', 'wb') as f:
                        pickle.dump(to_run_rmsd, f)

                else:  # For RMSF and 2D-RMSD, we need to check if the aligned trajectory already exists.
                    if os.path.exists(aligned_traj_file):
                        print(f"Using pre-aligned trajectory: {aligned_traj_file}")
                        u = mda.Universe(top_file, aligned_traj_file)
                        target_selection = u.select_atoms(target_selection)
                        ref_selection = u.select_atoms(ref_selection)
                        transform_trajectory(u, target_selection, ref_selection)
                    else:
                        print(f"Aligned trajectory file {aligned_traj_file} does not exist. Creating it now.")
                        u = mda.Universe(top_file, traj_file)
                        align_trajectory(u, target_selection, analysis, system, variation, rep, traj_format, start_frame)
                        del u
                        u = mda.Universe(top_file, aligned_traj_file)
                        target_selection = u.select_atoms(target_selection)
                        ref_selection = u.select_atoms(ref_selection)
                        transform_trajectory(u, target_selection, ref_selection)

                    if analysis == 'RMSF':
                        # TODO: Implement protein chain split to avoid considering the movements of all chains as a single entity.
                        print("Calculating RMSF...")
                        u_rmsf_backbone = u.select_atoms(target_selection)
                        to_run_rmsf = rms.RMSF(u_rmsf_backbone).run()

                        with open(f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}.pkl', 'wb') as f:
                            pickle.dump(to_run_rmsf, f)

                    elif analysis == '2D-RMSD':
                        print("Calculating 2D-RMSD...")
                        matrix_2d_rmsd = diffusionmap.DistanceMatrix(u, select=target_selection).run()

                        with open(f'{analysis_file_prefix[analysis]}_{system}_{variation}_rep{rep}.pkl', 'wb') as f:
                            pickle.dump(matrix_2d_rmsd, f)

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
    run_rms_analysis(
        systems=systems,
        variations=variations,
        num_replicates=args.num_replicates,
        analysis=args.analysis,
        target_selection=args.target_selection,
        ref_selection=args.ref_selection,
        start_frame=args.start_frame,
        traj_format=args.traj_format,
        group_selections=args.group_selections
    )

    return 0


if __name__ == '__main__':
    exit(main())

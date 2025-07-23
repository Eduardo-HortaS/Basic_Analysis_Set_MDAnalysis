import os
import pickle
import argparse
import json
import MDAnalysis as mda
import numpy as np
from utils import transform_trajectory

# Example selection strings for Radius of Gyration analysis:
#
# RoG Analysis - Calculate radius of gyration over time:
#   selection = 'protein and backbone' # Atoms to calculate RoG for (default)
#   selection = 'protein'              # All protein atoms
#   selection = 'nucleic'              # All nucleic acid atoms

def run_rog_analysis(systems, variations, num_replicates, start_frame, traj_format, selection='protein and backbone'):
    """
    Runs the Radius of Gyration analysis for each system and variation and saves results as individual pickle files.
    """
    print("--- Starting Radius of Gyration calculation for each system and variation... ---")
    analysis = 'RoG'

    reps = range(1, num_replicates + 1)
    analysis_file_prefix = 'rog_plot'

    for system in systems:
        for variation in variations[system]:
            for rep in reps:
                pickle_file = f'{analysis_file_prefix}_{system}_{variation}_rep{rep}.pkl'
                if os.path.exists(pickle_file):
                    print(f"Skipping {analysis} for {system}, {variation}, replicate {rep} because the data already exists in {pickle_file}.")
                    continue

                print(f"Processing {system}, {variation}, replicate {rep}.")
                traj_file = f'{system}/{variation}/{system}_production_{variation}_rep_{rep}.{traj_format}'
                top_file = f'{system}/{variation}/{system}_system_{variation}.top'

                u = mda.Universe(top_file, traj_file)

                # Transform trajectory
                protein = u.select_atoms('protein')
                not_protein = u.select_atoms('not protein')

                transform_trajectory(u, protein, not_protein)

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
                    times.append(ts.time / 1000.0)  # Convert ps to ns

                # Create results object compatible with plotting functions
                class RoGResults:
                    def __init__(self, frames, times, rog_values):
                        # Create structure similar to RMSD results
                        self.rog_data = np.column_stack((frames, times, rog_values))

                rog_analysis_results = RoGResults(frames, times, rog_results)

                with open(f'{analysis_file_prefix}_{system}_{variation}_rep{rep}.pkl', 'wb') as f:
                    pickle.dump(rog_analysis_results, f)

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
        selection=args.selection
    )

    return 0


if __name__ == '__main__':
    exit(main())

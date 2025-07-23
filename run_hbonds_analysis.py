import os
import pickle
import argparse
import json
import MDAnalysis as mda
from MDAnalysis.analysis import hydrogenbonds

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

def run_hbonds_analysis(systems, variations, num_replicates, d_a_cutoff, d_h_a_angle_cutoff, start_frame, traj_format, acceptors_sel=None, hydrogens_sel=None, between_pairs=None, update_selections=True):
    """
    Runs the Hydrogen Bonds analysis for each system and variation and saves results as individual pickle files.
    """
    print("--- Starting Hydrogen Bonds calculation for each system and variation... ---")

    reps = range(1, num_replicates + 1)

    for system in systems:
        for variation in variations[system]:
            for rep in reps:
                pickle_file = f'hbonds_plot_{system}_{variation}_rep{rep}.pkl'
                if os.path.exists(pickle_file):
                    print(f"Skipping Hydrogen Bonds analysis for {system}, {variation}, replicate {rep} because the data already exists in {pickle_file}.")
                    continue

                print(f"Processing {system}, {variation}, replicate {rep}.")
                traj_file = f'{system}/{variation}/{system}_production_{variation}_rep_{rep}.{traj_format}'
                top_file = f'{system}/{variation}/{system}_system_{variation}.top'
                u = mda.Universe(top_file, traj_file)

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

                hbonds.run(start=start_frame, stop=None, step=1)

                with open(f'hbonds_plot_{system}_{variation}_rep{rep}.pkl', 'wb') as f:
                    pickle.dump(hbonds, f)

                # TODO: Utilize the results from hbonds with proper logic for generating plots such as:
                # hydrogen bonds count by time, count by type, count by ids.
                # Reference: https://userguide.mdanalysis.org/stable/examples/analysis/hydrogen_bonds/hbonds.html

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
        acceptors_sel=acceptors_sel,
        hydrogens_sel=hydrogens_sel,
        between_pairs=between_pairs,
        update_selections=args.update_selections
    )

    return 0


if __name__ == '__main__':
    exit(main())

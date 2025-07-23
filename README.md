# Basic_Analysis_Set_MDAnalysis

A set of scripts using MDAnalysis that works on the following assumptions to calculate **RMSD**, **2D-RMSD**, **RMSF**, **RoG** and **HB**.

## Analysis Requirements

1. **Basic Configuration**
   The analysis will at least consist of x systems, y variations, z replicates, topology format and trajectory format, and start frame (ex. 500).

   ### Hydrogen Bonds Analysis

   Add one of two pairs of atom selections:

   - **Option 1:** Atom-focused analysis
     - `acceptors_sel` (e.g. `'protein and name O*'`)
     - `hydrogens_sel` (e.g. `'nucleic and name H*'`)

   - **Option 2:** Pair-focused analysis
     - `between_pairs` list of lists (e.g. `[['protein and resnum 73', 'nucleic and name NH*'], ['protein and resnum 73', 'nucleic and name N*'], ...]`)

   **Additional parameters:**
   - `d_a_cutoff` (e.g. 3.5)
   - `d_h_a_angle_cutoff` (e.g. 150)
   - `update_selections` may be passed as `False` to improve performance, but this will keep selection from updating over time, which is generally not recommended.

   ### RMS* Analysis (RMSF, RMSD or 2D-RMSD)

   - `reference_selection` (e.g. `'protein and backbone'`)
   - `target_selection` (e.g. `''`)

2. **Analysis Definition**
   The specific analysis must be defined in the `analysis_to_plot_prefix` dictionary.

3. **Plot Generation - All Variations**
   When making `{analysis}_{system}_all_variations_rep*.png` files, the script will take in all y variations for a given replicate number.

4. **Plot Generation - Comparisons**
   When making `{analysis}_{system}_{variation}_vs_{equivalent_system}_{equivalent_variation}_rep*.png` files and the variation is an actual variation (!= 'wild'), the script will take the equivalent variation from the `equivalent_mutantions` dictionary.

5. **Directory Structure**
   The script will be run on the parent folder for both systems, where each has its own folder with the following structure:

   ```text
   ./{system}/{variation}/
   ```

## Configuration

Configuration is delegated to a JSON.



[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16389354.svg)](https://doi.org/10.5281/zenodo.16389354)

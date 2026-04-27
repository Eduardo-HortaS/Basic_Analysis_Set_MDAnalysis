"""
2D-RMSD plotting script.

Generates 2D-RMSD distance matrix heatmaps.
Reads pickle files produced by run_rms_analysis.py (2D-RMSD mode).
"""
import os
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import DEFAULT_DPI, prettify_label
from style import format_selection_context


def plot_2d_rmsd(pickle_file, output_dir='.', dpi=DEFAULT_DPI, cmap='viridis', label=None):
    """
    Generates a 2D-RMSD heatmap from a pickle file.

    Parameters
    ----------
    pickle_file : str
        Path to the 2D-RMSD pickle file.
    output_dir : str
        Directory to save the plot.
    dpi : int
        Output resolution.
    cmap : str
        Matplotlib colormap name.
    label : str, optional
        Custom title suffix.
    """
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)

    selection = ''
    ref_selection = ''
    if isinstance(data, dict):
        if 'dist_matrix' not in data:
            raise ValueError("Invalid 2D-RMSD pickle schema: expected 'dist_matrix' key")
        dist_matrix = np.asarray(data['dist_matrix'])
        selection = data.get('selection', '')
        ref_selection = data.get('ref_selection', '')
    else:
        dist_matrix = data.results.dist_matrix

    if label is None:
        label = os.path.splitext(os.path.basename(pickle_file))[0].replace('2d_rmsd_plot_', '')
    label = prettify_label(label)

    fig, ax = plt.subplots(figsize=(10, 8))

    im = ax.imshow(dist_matrix, cmap=cmap, origin='lower', aspect='equal')
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label(r'RMSD ($\AA$)', fontsize=13)

    ax.set_xlabel('Frame', fontsize=14, fontweight='bold')
    ax.set_ylabel('Frame', fontsize=14, fontweight='bold')
    ax.set_title(f'2D-RMSD Distance Matrix\n{label}', fontsize=16, fontweight='bold', pad=15)
    context_line = format_selection_context(
        target_selection=selection,
        ref_selection=ref_selection,
    )
    if context_line:
        ax.text(
            0.5,
            1.01,
            f'Selections: {context_line}',
            transform=ax.transAxes,
            fontsize=10,
            ha='center',
            va='bottom',
            style='italic',
            color='gray',
        )
    ax.tick_params(axis='both', labelsize=12)

    os.makedirs(output_dir, exist_ok=True)
    output_name = os.path.splitext(os.path.basename(pickle_file))[0] + '.png'
    output_path = os.path.join(output_dir, output_name)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved 2D-RMSD plot to {output_path}")


def main():
    """Main function to parse arguments and generate 2D-RMSD plots."""
    parser = argparse.ArgumentParser(description='Generate 2D-RMSD heatmap plots from pickle files')

    parser.add_argument('--pickle-file', type=str, required=True,
                        help='Path to the 2D-RMSD pickle file')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory to save plots (default: current directory)')
    parser.add_argument('--dpi', type=int, default=DEFAULT_DPI,
                        help=f'Output resolution in DPI (default: {DEFAULT_DPI})')
    parser.add_argument('--cmap', type=str, default='viridis',
                        help='Matplotlib colormap (default: viridis)')
    parser.add_argument('--label', type=str, default=None,
                        help='Custom title suffix')

    args = parser.parse_args()

    if not os.path.exists(args.pickle_file):
        print(f"Error: Pickle file not found: {args.pickle_file}")
        return 1

    plot_2d_rmsd(
        pickle_file=args.pickle_file,
        output_dir=args.output_dir,
        dpi=args.dpi,
        cmap=args.cmap,
        label=args.label
    )

    return 0


if __name__ == '__main__':
    exit(main())

"""
Radius of Gyration (RoG) plotting script.

Generates RoG vs Time plots with an optional KDE density sideplot.
Inspired by paulus_plot_rmsd.py: dual-panel layout with density estimation.

Reads pickle files produced by run_rog_analysis.py.
"""
import os
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.stats import gaussian_kde

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import apply_style, get_color_cycle, format_label_with_stats, prettify_label, DEFAULT_DPI, DEFAULT_FIGSIZE

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import time_unit_label


def plot_rog(pickle_file, output_dir='.', dpi=DEFAULT_DPI, show_kde=True, label=None):
    """
    Generates a Radius of Gyration vs Time plot from a pickle file.

    Parameters
    ----------
    pickle_file : str
        Path to the RoG pickle file.
    output_dir : str
        Directory to save the plot.
    dpi : int
        Output resolution.
    show_kde : bool
        Whether to show the KDE density sideplot.
    label : str, optional
        Custom label for the plot.
    """
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)

    # Handle both old-style (bare RoGResults) and new-style (dict with metadata)
    if isinstance(data, dict):
        rog_data = data['rog_obj'].rog_data
        time_unit = data.get('time_unit', 'ns')
        selection = data.get('selection', '')
    else:
        # Legacy pickle — plain RoGResults object, times already in ns
        rog_data = data.rog_data
        time_unit = 'ns'
        selection = ''

    time_values = rog_data[:, 1]
    rog_values = rog_data[:, 2]

    if label is None:
        label = os.path.splitext(os.path.basename(pickle_file))[0].replace('rog_plot_', '')
    label = prettify_label(label)

    if show_kde:
        fig = plt.figure(figsize=DEFAULT_FIGSIZE)
        gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1], wspace=0.05)
        ax0 = plt.subplot(gs[0])
    else:
        fig, ax0 = plt.subplots(figsize=DEFAULT_FIGSIZE)

    color = get_color_cycle(1)[0]
    mean_rog = np.mean(rog_values)
    plot_label = format_label_with_stats(label, rog_values)

    ax0.plot(time_values, rog_values, color=color, label=plot_label, lw=2.5, alpha=0.9)
    ax0.axhline(mean_rog, color=color, linestyle='--', lw=1.2, alpha=0.5)

    ax0.set_xlabel(f'Time ({time_unit_label(time_unit)})', fontsize=14, fontweight='bold')
    ax0.set_ylabel(r'Radius of Gyration ($\AA$)', fontsize=14, fontweight='bold')
    ax0.set_title('Radius of Gyration vs Time', fontsize=16, fontweight='bold', pad=15)
    if selection:
        ax0.text(0.5, 1.01, f'Selection: {selection}',
                 transform=ax0.transAxes, fontsize=10, ha='center', va='bottom',
                 style='italic', color='gray')
    ax0.legend(loc='lower right', fontsize=11, frameon=True)
    apply_style(ax0)

    if show_kde:
        ax1 = plt.subplot(gs[1])
        rog_range = np.linspace(rog_values.min(), rog_values.max(), 200)
        kde = gaussian_kde(rog_values)
        ax1.plot(kde(rog_range), rog_range, color=color, lw=2.2)
        ax1.set_xlabel('Density', fontsize=13)
        ax1.set_title('RoG Density', fontsize=14, pad=10)
        ax1.set_ylabel(r'RoG ($\AA$)', fontsize=13)
        apply_style(ax1, remove_spines=['top', 'left', 'bottom'])

    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.13, top=0.92, wspace=0.05)

    os.makedirs(output_dir, exist_ok=True)
    output_name = os.path.splitext(os.path.basename(pickle_file))[0] + '.png'
    output_path = os.path.join(output_dir, output_name)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved RoG plot to {output_path}")


def main():
    """Main function to parse arguments and generate RoG plots."""
    parser = argparse.ArgumentParser(description='Generate Radius of Gyration plots from pickle files')

    parser.add_argument('--pickle-file', type=str, required=True,
                        help='Path to the RoG pickle file')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory to save plots (default: current directory)')
    parser.add_argument('--dpi', type=int, default=DEFAULT_DPI,
                        help=f'Output resolution in DPI (default: {DEFAULT_DPI})')
    parser.add_argument('--no-kde', action='store_true',
                        help='Disable the KDE density sideplot')
    parser.add_argument('--label', type=str, default=None,
                        help='Custom label for the plot legend')

    args = parser.parse_args()

    if not os.path.exists(args.pickle_file):
        print(f"Error: Pickle file not found: {args.pickle_file}")
        return 1

    plot_rog(
        pickle_file=args.pickle_file,
        output_dir=args.output_dir,
        dpi=args.dpi,
        show_kde=not args.no_kde,
        label=args.label
    )

    return 0


if __name__ == '__main__':
    exit(main())

"""
RMSD plotting script.

Generates RMSD vs Time plots with an optional KDE density sideplot.
Inspired by paulus_plot_rmsd.py: dual-panel GridSpec layout, mean overlay lines,
legend with statistics, and professional styling.

Reads pickle files produced by run_rms_analysis.py (RMSD mode).
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
from style import (apply_style, get_color_cycle, format_label_with_stats,
                   format_selection_subtitle, prettify_label, DEFAULT_DPI, DEFAULT_FIGSIZE)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import time_unit_label


def plot_rmsd(pickle_file, output_dir='.', dpi=DEFAULT_DPI, show_kde=True, label=None):
    """
    Generates an RMSD vs Time plot from a pickle file.

    Parameters
    ----------
    pickle_file : str
        Path to the RMSD pickle file.
    output_dir : str
        Directory to save the plot.
    dpi : int
        Output resolution.
    show_kde : bool
        Whether to show the KDE density sideplot.
    label : str, optional
        Custom label for the plot. If None, derived from filename.
    """
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)

    # Handle both portable dict schema and legacy object-based schemas.
    if isinstance(data, dict):
        if 'rmsd_data' in data:
            rmsd_data = np.asarray(data['rmsd_data'])
        elif 'rmsd_obj' in data:
            rmsd_data = data['rmsd_obj'].results.rmsd
        else:
            raise ValueError("Invalid RMSD pickle schema: expected 'rmsd_data' key")
        time_corrected = data.get('time_corrected', False)
        time_unit = data.get('time_unit', 'ps')
    else:
        rmsd_obj = data
        rmsd_data = rmsd_obj.results.rmsd
        time_corrected = False
        time_unit = 'ps'

    time_values = rmsd_data[:, 1]
    rmsd_values = rmsd_data[:, 2]

    if label is None:
        label = os.path.splitext(os.path.basename(pickle_file))[0].replace('rmsd_plot_', '')
    label = prettify_label(label)

    time_label = f'Time ({time_unit_label(time_unit)})'

    if show_kde:
        fig = plt.figure(figsize=DEFAULT_FIGSIZE)
        gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1], wspace=0.05)
        ax0 = plt.subplot(gs[0])
    else:
        fig, ax0 = plt.subplots(figsize=DEFAULT_FIGSIZE)

    # Main RMSD vs Time plot
    mean_rmsd = np.mean(rmsd_values)
    color = get_color_cycle(1)[0]
    plot_label = format_label_with_stats(label, rmsd_values)

    ax0.plot(time_values, rmsd_values, color=color, label=plot_label, lw=2.5, alpha=0.9)
    ax0.axhline(mean_rmsd, color=color, linestyle='--', lw=1.2, alpha=0.5)

    ax0.set_xlabel(time_label, fontsize=14, fontweight='bold')
    ax0.set_ylabel(r'RMSD ($\AA$)', fontsize=14, fontweight='bold')
    ax0.set_title('RMSD vs Time', fontsize=16, fontweight='bold', pad=15)
    ax0.legend(loc='lower right', fontsize=11, frameon=True)
    apply_style(ax0)

    if show_kde:
        # KDE density sideplot
        ax1 = plt.subplot(gs[1])
        rmsd_range = np.linspace(rmsd_values.min(), rmsd_values.max(), 200)
        kde = gaussian_kde(rmsd_values)
        ax1.plot(kde(rmsd_range), rmsd_range, color=color, lw=2.2)
        ax1.set_xlabel('Density', fontsize=13)
        ax1.set_title('RMSD Density', fontsize=14, pad=10)
        ax1.set_ylabel(r'RMSD ($\AA$)', fontsize=13)
        apply_style(ax1, remove_spines=['top', 'left', 'bottom'])

    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.13, top=0.92, wspace=0.05)

    os.makedirs(output_dir, exist_ok=True)
    output_name = os.path.splitext(os.path.basename(pickle_file))[0] + '.png'
    output_path = os.path.join(output_dir, output_name)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved RMSD plot to {output_path}")


def _load_rmsd_pickle(pickle_file):
    """Load an RMSD pickle file and return (time_values, rmsd_values, time_unit, selection, ref_selection)."""
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)

    if isinstance(data, dict):
        if 'rmsd_data' in data:
            rmsd_data = np.asarray(data['rmsd_data'])
        elif 'rmsd_obj' in data:
            rmsd_data = data['rmsd_obj'].results.rmsd
        else:
            raise ValueError("Invalid RMSD pickle schema: expected 'rmsd_data' key")
        time_unit = data.get('time_unit', 'ps')
        selection = data.get('selection', '')
        ref_selection = data.get('ref_selection', '')
    else:
        rmsd_obj = data
        rmsd_data = rmsd_obj.results.rmsd
        time_unit = 'ps'
        selection = ''
        ref_selection = ''

    return rmsd_data[:, 1], rmsd_data[:, 2], time_unit, selection, ref_selection


def plot_rmsd_comparison(pickle_files, labels, group_name, output_dir='.',
                         dpi=DEFAULT_DPI, show_kde=True, selection_label=None):
    """
    Overlay multiple RMSD traces on a single plot for comparison.

    Parameters
    ----------
    pickle_files : list of str
        Paths to RMSD pickle files to overlay.
    labels : list of str
        Labels for each trace (e.g. "SystemA / wild").
    group_name : str
        Name of the comparison group (used in title and filename).
    output_dir : str
        Directory to save the plot.
    dpi : int
        Output resolution.
    show_kde : bool
        Whether to show KDE density sideplots for each trace.
    selection_label : str, optional
        The selection string used; shown as subtitle.
    """
    n = len(pickle_files)
    colors = get_color_cycle(n)

    if show_kde:
        fig = plt.figure(figsize=DEFAULT_FIGSIZE)
        gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1], wspace=0.05)
        ax0 = plt.subplot(gs[0])
    else:
        fig, ax0 = plt.subplots(figsize=DEFAULT_FIGSIZE)

    time_unit_str = 'ps'

    for i, (pkl, lbl) in enumerate(zip(pickle_files, labels)):
        time_values, rmsd_values, tu, sel, _ref_sel = _load_rmsd_pickle(pkl)
        time_unit_str = tu
        color = colors[i]
        mean_rmsd = np.mean(rmsd_values)
        plot_label = format_label_with_stats(lbl, rmsd_values)

        ax0.plot(time_values, rmsd_values, color=color, label=plot_label,
                 lw=2.0, alpha=0.85)
        ax0.axhline(mean_rmsd, color=color, linestyle='--', lw=1.0, alpha=0.4)

    time_label = f'Time ({time_unit_label(time_unit_str)})'
    ax0.set_xlabel(time_label, fontsize=14, fontweight='bold')
    ax0.set_ylabel(r'RMSD ($\AA$)', fontsize=14, fontweight='bold')

    title = f'RMSD vs Time \u2014 {prettify_label(group_name)}'
    ax0.set_title(title, fontsize=16, fontweight='bold', pad=15)

    if selection_label:
        ax0.text(0.5, 1.01, f'Selection: {format_selection_subtitle(selection_label)}',
                 transform=ax0.transAxes, fontsize=10, ha='center', va='bottom',
                 style='italic', color='gray')

    ax0.legend(loc='lower right', fontsize=10, frameon=True)
    apply_style(ax0)

    if show_kde:
        ax1 = plt.subplot(gs[1])
        all_rmsd_min = float('inf')
        all_rmsd_max = float('-inf')
        kde_data = []

        for i, (pkl, lbl) in enumerate(zip(pickle_files, labels)):
            _, rmsd_values, _, _, _ = _load_rmsd_pickle(pkl)
            all_rmsd_min = min(all_rmsd_min, rmsd_values.min())
            all_rmsd_max = max(all_rmsd_max, rmsd_values.max())
            kde_data.append((rmsd_values, colors[i]))

        rmsd_range = np.linspace(all_rmsd_min, all_rmsd_max, 200)
        for rmsd_values, color in kde_data:
            kde = gaussian_kde(rmsd_values)
            ax1.plot(kde(rmsd_range), rmsd_range, color=color, lw=1.8, alpha=0.8)

        ax1.set_xlabel('Density', fontsize=13)
        ax1.set_title('RMSD Density', fontsize=14, pad=10)
        ax1.set_ylabel(r'RMSD ($\AA$)', fontsize=13)
        apply_style(ax1, remove_spines=['top', 'left', 'bottom'])

    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.13, top=0.89, wspace=0.05)

    os.makedirs(output_dir, exist_ok=True)
    safe_name = group_name.replace(' ', '_').replace('/', '_')
    output_path = os.path.join(output_dir, f'rmsd_comparison_{safe_name}.png')
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved RMSD comparison plot to {output_path}")


def plot_rmsd_comparison_average(pickle_groups, labels, group_name, output_dir='.',
                                 dpi=DEFAULT_DPI, show_kde=True, selection_label=None):
    """
    Overlay averaged RMSD traces (mean ± std across replicates) on a single plot.

    Parameters
    ----------
    pickle_groups : list of list of str
        For each trace, a list of replicate pickle files to average.
    labels : list of str
        Labels for each averaged trace.
    group_name : str
        Name of the comparison group.
    output_dir : str
        Directory to save the plot.
    dpi : int
        Output resolution.
    show_kde : bool
        Whether to show KDE density sideplot (uses concatenated data).
    selection_label : str, optional
        Selection string shown as subtitle.
    """
    n = len(pickle_groups)
    colors = get_color_cycle(n)

    if show_kde:
        fig = plt.figure(figsize=DEFAULT_FIGSIZE)
        gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1], wspace=0.05)
        ax0 = plt.subplot(gs[0])
    else:
        fig, ax0 = plt.subplots(figsize=DEFAULT_FIGSIZE)

    time_unit_str = 'ps'

    for i, (pkl_list, lbl) in enumerate(zip(pickle_groups, labels)):
        # Load all replicates and interpolate to common time grid
        all_times = []
        all_rmsd = []
        for pkl in pkl_list:
            tv, rv, tu, sel, _ref_sel = _load_rmsd_pickle(pkl)
            time_unit_str = tu
            all_times.append(tv)
            all_rmsd.append(rv)

        # Use the shortest common time range
        t_min = max(t.min() for t in all_times)
        t_max = min(t.max() for t in all_times)
        n_pts = min(len(t) for t in all_times)
        common_time = np.linspace(t_min, t_max, n_pts)

        # Interpolate each replicate to common grid
        interpolated = []
        for tv, rv in zip(all_times, all_rmsd):
            interpolated.append(np.interp(common_time, tv, rv))

        stacked = np.array(interpolated)
        mean_rmsd = np.mean(stacked, axis=0)
        std_rmsd = np.std(stacked, axis=0)

        color = colors[i]
        plot_label = format_label_with_stats(lbl, mean_rmsd)

        ax0.plot(common_time, mean_rmsd, color=color, label=plot_label, lw=2.0, alpha=0.9)
        ax0.fill_between(common_time, mean_rmsd - std_rmsd, mean_rmsd + std_rmsd,
                         color=color, alpha=0.2)
        ax0.axhline(np.mean(mean_rmsd), color=color, linestyle='--', lw=1.0, alpha=0.4)

    time_label = f'Time ({time_unit_label(time_unit_str)})'
    ax0.set_xlabel(time_label, fontsize=14, fontweight='bold')
    ax0.set_ylabel(r'RMSD ($\AA$)', fontsize=14, fontweight='bold')

    title = f'RMSD vs Time \u2014 {prettify_label(group_name)} (averaged)'
    ax0.set_title(title, fontsize=16, fontweight='bold', pad=15)

    if selection_label:
        ax0.text(0.5, 1.01, f'Selection: {format_selection_subtitle(selection_label)}',
                 transform=ax0.transAxes, fontsize=10, ha='center', va='bottom',
                 style='italic', color='gray')

    ax0.legend(loc='lower right', fontsize=10, frameon=True)
    apply_style(ax0)

    if show_kde:
        ax1 = plt.subplot(gs[1])
        all_rmsd_min = float('inf')
        all_rmsd_max = float('-inf')
        kde_data_list = []

        for i, (pkl_list, lbl) in enumerate(zip(pickle_groups, labels)):
            concatenated = np.concatenate([_load_rmsd_pickle(p)[1] for p in pkl_list])
            all_rmsd_min = min(all_rmsd_min, concatenated.min())
            all_rmsd_max = max(all_rmsd_max, concatenated.max())
            kde_data_list.append((concatenated, colors[i]))

        rmsd_range = np.linspace(all_rmsd_min, all_rmsd_max, 200)
        for concat_vals, color in kde_data_list:
            kde = gaussian_kde(concat_vals)
            ax1.plot(kde(rmsd_range), rmsd_range, color=color, lw=1.8, alpha=0.8)

        ax1.set_xlabel('Density', fontsize=13)
        ax1.set_title('RMSD Density', fontsize=14, pad=10)
        ax1.set_ylabel(r'RMSD ($\AA$)', fontsize=13)
        apply_style(ax1, remove_spines=['top', 'left', 'bottom'])

    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.13, top=0.89, wspace=0.05)

    os.makedirs(output_dir, exist_ok=True)
    safe_name = group_name.replace(' ', '_').replace('/', '_')
    output_path = os.path.join(output_dir, f'rmsd_comparison_{safe_name}_avg.png')
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved averaged RMSD comparison plot to {output_path}")


def main():
    """Main function to parse arguments and generate RMSD plots."""
    parser = argparse.ArgumentParser(description='Generate RMSD plots from pickle files')

    parser.add_argument('--pickle-file', type=str, required=True,
                        help='Path to the RMSD pickle file')
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

    plot_rmsd(
        pickle_file=args.pickle_file,
        output_dir=args.output_dir,
        dpi=args.dpi,
        show_kde=not args.no_kde,
        label=args.label
    )

    return 0


if __name__ == '__main__':
    exit(main())

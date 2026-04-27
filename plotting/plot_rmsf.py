"""
RMSF plotting script.

Generates per-residue RMSF bar/line plots.
Supports both whole-protein and per-chain RMSF pickle files.

Reads pickle files produced by run_rms_analysis.py (RMSF mode).
"""
import os
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import (apply_style, get_color_cycle, format_label_with_stats,
                   format_selection_subtitle, format_selection_context,
                   prettify_label, DEFAULT_DPI, DEFAULT_FIGSIZE)


def plot_rmsf(pickle_file, output_dir='.', dpi=DEFAULT_DPI, plot_type='line', label=None):
    """
    Generates a per-residue RMSF plot from a pickle file.

    Parameters
    ----------
    pickle_file : str
        Path to the RMSF pickle file.
    output_dir : str
        Directory to save the plot.
    dpi : int
        Output resolution.
    plot_type : str
        'line' or 'bar' plot type.
    label : str, optional
        Custom label for the plot.
    """
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)

    # Handle chain-split (dict) vs. full RMSF (rms.RMSF object) pickles
    selection = ''
    ref_selection = ''
    if isinstance(data, dict) and 'rmsf' in data:
        # Chain-split pickle
        rmsf_values = data['rmsf']
        resids = data['resids']
        chain_id = data.get('chain_id', '')
        original_resids = data.get('original_resids', resids)
        selection = data.get('selection', '')
        ref_selection = data.get('ref_selection', '')
        if label is None:
            label = f"Chain {chain_id}" if chain_id else os.path.basename(pickle_file)
        title_suffix = f" - Chain {chain_id}" if chain_id else ""
    else:
        # Full RMSF object
        rmsf_values = data.results.rmsf
        resids = np.arange(1, len(rmsf_values) + 1)
        chain_id = ''
        title_suffix = ''
        if label is None:
            label = os.path.splitext(os.path.basename(pickle_file))[0].replace('rmsf_plot_', '')
    label = prettify_label(label)

    color = get_color_cycle(1)[0]

    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    if plot_type == 'bar':
        ax.bar(resids, rmsf_values, color=color, alpha=0.7, width=1.0, label=label)
    else:
        ax.plot(resids, rmsf_values, color=color, lw=1.5, alpha=0.9, label=label)
        ax.fill_between(resids, rmsf_values, alpha=0.15, color=color)

    # Mean line
    mean_rmsf = np.mean(rmsf_values)
    ax.axhline(mean_rmsf, color=color, linestyle='--', lw=1.2, alpha=0.5,
               label=f'Mean = {mean_rmsf:.2f} \u00c5')

    ax.set_xlabel('Residue Number', fontsize=14, fontweight='bold')
    ax.set_ylabel(r'RMSF ($\AA$)', fontsize=14, fontweight='bold')
    ax.set_title(f'RMSF per Residue{title_suffix}', fontsize=16, fontweight='bold', pad=15)
    context_line = format_selection_context(
        target_selection=selection,
        ref_selection=ref_selection,
    )
    if context_line:
        ax.text(0.5, 1.01, f'Selections: {context_line}',
                transform=ax.transAxes, fontsize=10, ha='center', va='bottom',
                style='italic', color='gray')
    ax.legend(loc='upper right', fontsize=11, frameon=True)
    apply_style(ax)

    os.makedirs(output_dir, exist_ok=True)
    output_name = os.path.splitext(os.path.basename(pickle_file))[0] + '.png'
    output_path = os.path.join(output_dir, output_name)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved RMSF plot to {output_path}")


def _load_rmsf_pickle(pickle_file):
    """Load an RMSF pickle file and return data plus selection metadata."""
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)

    if isinstance(data, dict) and 'rmsf' in data:
        rmsf_values = data['rmsf']
        resids = data['resids']
        chain_id = data.get('chain_id', '')
        selection = data.get('selection', '')
        ref_selection = data.get('ref_selection', '')
    else:
        rmsf_values = data.results.rmsf
        resids = np.arange(1, len(rmsf_values) + 1)
        chain_id = ''
        selection = ''
        ref_selection = ''

    return resids, rmsf_values, chain_id, selection, ref_selection


def plot_rmsf_comparison(pickle_files, labels, group_name, output_dir='.',
                         dpi=DEFAULT_DPI, selection_label=None):
    """
    Overlay multiple RMSF traces on a single plot for comparison.

    No fill_between — line only for each trace, color-coded.

    Parameters
    ----------
    pickle_files : list of str
        Paths to RMSF pickle files to overlay.
    labels : list of str
        Labels for each trace (e.g. "SystemA / wild").
    group_name : str
        Name of the comparison group (used in title and filename).
    output_dir : str
        Directory to save the plot.
    dpi : int
        Output resolution.
    selection_label : str, optional
        The selection string used; shown as subtitle.
    """
    n = len(pickle_files)
    colors = get_color_cycle(n)

    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    chain_id_for_title = ''
    plot_selection = ''
    plot_ref_selection = ''
    for i, (pkl, lbl) in enumerate(zip(pickle_files, labels)):
        resids, rmsf_values, chain_id, selection, ref_selection = _load_rmsf_pickle(pkl)
        if chain_id:
            chain_id_for_title = chain_id
        if not plot_selection:
            plot_selection = selection
        if not plot_ref_selection:
            plot_ref_selection = ref_selection
        color = colors[i]
        plot_label = format_label_with_stats(lbl, rmsf_values)

        ax.plot(resids, rmsf_values, color=color, lw=1.5, alpha=0.85, label=plot_label)
        mean_rmsf = np.mean(rmsf_values)
        ax.axhline(mean_rmsf, color=color, linestyle='--', lw=1.0, alpha=0.4)

    ax.set_xlabel('Residue Number', fontsize=14, fontweight='bold')
    ax.set_ylabel(r'RMSF ($\AA$)', fontsize=14, fontweight='bold')

    title_suffix = f' \u2014 Chain {chain_id_for_title}' if chain_id_for_title else ''
    title = f'RMSF per Residue \u2014 {prettify_label(group_name)}{title_suffix}'
    ax.set_title(title, fontsize=16, fontweight='bold', pad=15)

    context_line = format_selection_context(
        target_selection=plot_selection or selection_label,
        ref_selection=plot_ref_selection,
    )
    if context_line:
        ax.text(0.5, 1.01, f'Selections: {format_selection_subtitle(context_line, max_length=120)}',
                transform=ax.transAxes, fontsize=10, ha='center', va='bottom',
                style='italic', color='gray')

    ax.legend(loc='upper right', fontsize=10, frameon=True)
    apply_style(ax)

    os.makedirs(output_dir, exist_ok=True)
    safe_name = group_name.replace(' ', '_').replace('/', '_')
    chain_suffix = f'_chain{chain_id_for_title}' if chain_id_for_title else ''
    output_path = os.path.join(output_dir, f'rmsf_comparison_{safe_name}{chain_suffix}.png')
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved RMSF comparison plot to {output_path}")


def plot_rmsf_comparison_average(pickle_groups, labels, group_name, output_dir='.',
                                  dpi=DEFAULT_DPI, selection_label=None):
    """
    Overlay averaged RMSF traces (mean ± std across replicates) on a single plot.

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
    selection_label : str, optional
        Selection string shown as subtitle.
    """
    n = len(pickle_groups)
    colors = get_color_cycle(n)

    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    chain_id_for_title = ''
    plot_selection = ''
    plot_ref_selection = ''
    for i, (pkl_list, lbl) in enumerate(zip(pickle_groups, labels)):
        all_rmsf = []
        all_resids = None
        for pkl in pkl_list:
            resids, rmsf_values, chain_id, selection, ref_selection = _load_rmsf_pickle(pkl)
            if chain_id:
                chain_id_for_title = chain_id
            if not plot_selection:
                plot_selection = selection
            if not plot_ref_selection:
                plot_ref_selection = ref_selection
            all_rmsf.append(rmsf_values)
            if all_resids is None:
                all_resids = resids

        # Use shortest common length
        min_len = min(len(r) for r in all_rmsf)
        stacked = np.array([r[:min_len] for r in all_rmsf])
        mean_rmsf = np.mean(stacked, axis=0)
        std_rmsf = np.std(stacked, axis=0)
        resids = all_resids[:min_len]

        color = colors[i]
        plot_label = format_label_with_stats(lbl, mean_rmsf)

        ax.plot(resids, mean_rmsf, color=color, lw=1.5, alpha=0.9, label=plot_label)
        ax.fill_between(resids, mean_rmsf - std_rmsf, mean_rmsf + std_rmsf,
                         color=color, alpha=0.2)
        ax.axhline(np.mean(mean_rmsf), color=color, linestyle='--', lw=1.0, alpha=0.4)

    ax.set_xlabel('Residue Number', fontsize=14, fontweight='bold')
    ax.set_ylabel(r'RMSF ($\AA$)', fontsize=14, fontweight='bold')

    title_suffix = f' \u2014 Chain {chain_id_for_title}' if chain_id_for_title else ''
    title = f'RMSF per Residue \u2014 {prettify_label(group_name)} (averaged){title_suffix}'
    ax.set_title(title, fontsize=16, fontweight='bold', pad=15)

    context_line = format_selection_context(
        target_selection=plot_selection or selection_label,
        ref_selection=plot_ref_selection,
    )
    if context_line:
        ax.text(0.5, 1.01, f'Selections: {format_selection_subtitle(context_line, max_length=120)}',
                transform=ax.transAxes, fontsize=10, ha='center', va='bottom',
                style='italic', color='gray')

    ax.legend(loc='upper right', fontsize=10, frameon=True)
    apply_style(ax)

    os.makedirs(output_dir, exist_ok=True)
    safe_name = group_name.replace(' ', '_').replace('/', '_')
    chain_suffix = f'_chain{chain_id_for_title}' if chain_id_for_title else ''
    output_path = os.path.join(output_dir, f'rmsf_comparison_{safe_name}_avg{chain_suffix}.png')
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved averaged RMSF comparison plot to {output_path}")


def main():
    """Main function to parse arguments and generate RMSF plots."""
    parser = argparse.ArgumentParser(description='Generate RMSF plots from pickle files')

    parser.add_argument('--pickle-file', type=str, required=True,
                        help='Path to the RMSF pickle file')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory to save plots (default: current directory)')
    parser.add_argument('--dpi', type=int, default=DEFAULT_DPI,
                        help=f'Output resolution in DPI (default: {DEFAULT_DPI})')
    parser.add_argument('--plot-type', type=str, default='line', choices=['line', 'bar'],
                        help='Plot type: line or bar (default: line)')
    parser.add_argument('--label', type=str, default=None,
                        help='Custom label for the plot legend')

    args = parser.parse_args()

    if not os.path.exists(args.pickle_file):
        print(f"Error: Pickle file not found: {args.pickle_file}")
        return 1

    plot_rmsf(
        pickle_file=args.pickle_file,
        output_dir=args.output_dir,
        dpi=args.dpi,
        plot_type=args.plot_type,
        label=args.label
    )

    return 0


if __name__ == '__main__':
    exit(main())

"""
Hydrogen Bonds plotting script.

Generates three types of H-bond plots from pickle files produced by run_hbonds_analysis.py:
  1. Count by time  — H-bond count vs. time (using hbonds.count_by_time())
  2. Count by type  — bar chart of H-bond types (using hbonds.count_by_type())
  3. Count by ids   — top-N most frequent H-bonds (using hbonds.count_by_ids())

Reference: https://userguide.mdanalysis.org/stable/examples/analysis/hydrogen_bonds/hbonds.html
"""
import os
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import apply_style, get_color_cycle, prettify_label, DEFAULT_DPI, DEFAULT_FIGSIZE

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import convert_time_from_ps, time_unit_label


def _load_hbonds_payload(pickle_file):
    """Load H-bonds payload and provide actionable errors for legacy pickles."""
    try:
        with open(pickle_file, 'rb') as f:
            return pickle.load(f)
    except (OSError, FileNotFoundError, AttributeError) as exc:
        raise ValueError(
            "Failed to load H-bonds pickle. The file appears to contain "
            "non-portable MDAnalysis objects that require original trajectory "
            "paths during unpickling. Regenerate this pickle with the current "
            "run_hbonds_analysis.py, then re-run plotting."
        ) from exc


def plot_hbonds_by_time(times, counts, output_path, dpi=DEFAULT_DPI, time_unit='ps', label=''):
    """
    Plot hydrogen bond count over time.

    Parameters
    ----------
    times : array-like
        Time values in picoseconds.
    counts : array-like
        Hydrogen bond counts per frame.
    output_path : str
        Full path for the output PNG file.
    dpi : int
        Output resolution.
    time_unit : str
        Time unit for the x-axis.  Accepted: 'fs', 'ps', 'ns', 'us', 'ms', 's'.
    label : str
        Label for the plot title.
    """
    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
    color = get_color_cycle(1)[0]

    times = np.asarray(times)
    counts = np.asarray(counts)

    times = convert_time_from_ps(times, time_unit)

    ax.plot(times, counts, color=color, lw=2.0, alpha=0.9)

    mean_count = np.mean(counts)
    ax.axhline(mean_count, color=color, linestyle='--', lw=1.2, alpha=0.5,
               label=f'Mean = {mean_count:.1f}')

    ax.set_xlabel(f'Time ({time_unit_label(time_unit)})', fontsize=14, fontweight='bold')
    ax.set_ylabel(r'$N_{HB}$', fontsize=14, fontweight='bold')
    ax.set_title(f'Hydrogen Bonds over Time\n{label}', fontsize=16, fontweight='bold', pad=15)
    ax.legend(loc='upper right', fontsize=11, frameon=True)
    apply_style(ax)

    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved H-bonds by time plot to {output_path}")


def plot_hbonds_by_type(type_counts, output_path, dpi=DEFAULT_DPI, label=''):
    """
    Plot hydrogen bond counts grouped by donor-acceptor type.

    Parameters
    ----------
    type_counts : array-like
        Output array from ``count_by_type``.
    output_path : str
        Full path for the output PNG file.
    dpi : int
        Output resolution.
    label : str
        Label for the plot title.
    """
    type_counts = np.asarray(type_counts)

    if len(type_counts) == 0:
        print(f"No hydrogen bond types found for {label}. Skipping count-by-type plot.")
        return

    donors = [row[0] for row in type_counts]
    acceptors = [row[1] for row in type_counts]
    counts = [int(row[2]) for row in type_counts]

    type_labels = [f"{d} \u2192 {a}" for d, a in zip(donors, acceptors)]
    colors = get_color_cycle(len(type_labels))

    fig, ax = plt.subplots(figsize=(max(10, len(type_labels) * 1.2), 6))

    bars = ax.bar(range(len(type_labels)), counts, color=colors, alpha=0.8)

    ax.set_xticks(range(len(type_labels)))
    ax.set_xticklabels(type_labels, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Total Count', fontsize=14, fontweight='bold')
    ax.set_title(f'Hydrogen Bonds by Type\n{label}', fontsize=16, fontweight='bold', pad=15)
    apply_style(ax)

    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved H-bonds by type plot to {output_path}")


def plot_hbonds_by_ids(id_counts, output_path, dpi=DEFAULT_DPI, top_n=20, label=''):
    """
    Plot top-N most frequent hydrogen bonds by specific atom IDs.

    Parameters
    ----------
    id_counts : array-like
        Output array from ``count_by_ids``.
    output_path : str
        Full path for the output PNG file.
    dpi : int
        Output resolution.
    top_n : int
        Number of top H-bonds to show.
    label : str
        Label for the plot title.
    """
    id_counts = np.asarray(id_counts)

    if len(id_counts) == 0:
        print(f"No hydrogen bond IDs found for {label}. Skipping count-by-ids plot.")
        return

    # id_counts columns: [donor_ix, hydrogen_ix, acceptor_ix, count]
    # Already sorted by count (descending)
    top_data = id_counts[:top_n]

    hb_labels = [f"D:{int(row[0])}-H:{int(row[1])}-A:{int(row[2])}" for row in top_data]
    counts = [int(row[3]) for row in top_data]

    color = get_color_cycle(1)[0]

    fig, ax = plt.subplots(figsize=(10, max(6, len(hb_labels) * 0.4)))

    y_pos = range(len(hb_labels))
    ax.barh(y_pos, counts, color=color, alpha=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(hb_labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Total Count', fontsize=14, fontweight='bold')
    ax.set_title(f'Top {len(top_data)} Hydrogen Bonds by Atom IDs\n{label}', fontsize=16, fontweight='bold', pad=15)
    apply_style(ax)

    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved H-bonds by IDs plot to {output_path}")


def plot_hbonds(pickle_file, output_dir='.', dpi=DEFAULT_DPI, plot_types=None,
                top_n=20, time_unit='ps', label=None):
    """
    Main function to generate all requested H-bond plots from a pickle file.

    Parameters
    ----------
    pickle_file : str
        Path to the H-bonds pickle file.
    output_dir : str
        Directory to save plots.
    dpi : int
        Output resolution.
    plot_types : list of str
        Which plots to generate: 'time', 'type', 'ids'. Default: all three.
    top_n : int
        Number of top H-bonds for the by-ids plot.
    time_unit : str
        'ps' or 'ns' for the time axis.
    label : str, optional
        Custom label for plot titles.
    """
    if plot_types is None:
        plot_types = ['time', 'type', 'ids']

    hbonds = _load_hbonds_payload(pickle_file)

    if not isinstance(hbonds, dict):
        raise ValueError(
            "Unsupported H-bonds pickle schema. Expected a dict payload "
            "created by run_hbonds_analysis.py with portable arrays."
        )

    required_keys = {'times', 'count_by_time', 'count_by_type', 'count_by_ids'}
    missing = [k for k in required_keys if k not in hbonds]
    if missing:
        raise ValueError(
            f"Invalid H-bonds pickle schema. Missing keys: {', '.join(missing)}"
        )

    times = hbonds['times']
    count_by_time = hbonds['count_by_time']
    count_by_type = hbonds['count_by_type']
    count_by_ids = hbonds['count_by_ids']

    if label is None:
        label = os.path.splitext(os.path.basename(pickle_file))[0].replace('hbonds_plot_', '')
    label = prettify_label(label)

    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(pickle_file))[0]

    if 'time' in plot_types:
        out = os.path.join(output_dir, f'{base_name}_by_time.png')
        plot_hbonds_by_time(times, count_by_time, out, dpi=dpi, time_unit=time_unit, label=label)

    if 'type' in plot_types:
        out = os.path.join(output_dir, f'{base_name}_by_type.png')
        plot_hbonds_by_type(count_by_type, out, dpi=dpi, label=label)

    if 'ids' in plot_types:
        out = os.path.join(output_dir, f'{base_name}_by_ids.png')
        plot_hbonds_by_ids(count_by_ids, out, dpi=dpi, top_n=top_n, label=label)


def main():
    """Main function to parse arguments and generate H-bond plots."""
    parser = argparse.ArgumentParser(description='Generate Hydrogen Bond plots from pickle files')

    parser.add_argument('--pickle-file', type=str, required=True,
                        help='Path to the H-bonds pickle file')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory to save plots (default: current directory)')
    parser.add_argument('--dpi', type=int, default=DEFAULT_DPI,
                        help=f'Output resolution in DPI (default: {DEFAULT_DPI})')
    parser.add_argument('--plot-types', type=str, nargs='+', default=['time', 'type', 'ids'],
                        choices=['time', 'type', 'ids'],
                        help='Which plots to generate (default: all three)')
    parser.add_argument('--top-n', type=int, default=20,
                        help='Number of top H-bonds to show in by-ids plot (default: 20)')
    parser.add_argument('--time-unit', type=str, default='ps', choices=['ps', 'ns'],
                        help='Time unit for the time-axis plot (default: ps)')
    parser.add_argument('--label', type=str, default=None,
                        help='Custom label for plot titles')

    args = parser.parse_args()

    if not os.path.exists(args.pickle_file):
        print(f"Error: Pickle file not found: {args.pickle_file}")
        return 1

    plot_hbonds(
        pickle_file=args.pickle_file,
        output_dir=args.output_dir,
        dpi=args.dpi,
        plot_types=args.plot_types,
        top_n=args.top_n,
        time_unit=args.time_unit,
        label=args.label
    )

    return 0


if __name__ == '__main__':
    exit(main())

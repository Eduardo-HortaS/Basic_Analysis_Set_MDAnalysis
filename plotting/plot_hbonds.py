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
import re
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import (
    apply_style,
    get_color_cycle,
    prettify_label,
    DEFAULT_DPI,
    DEFAULT_FIGSIZE,
)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import convert_time_from_ps, time_unit_label, TIME_UNITS, validate_time_unit


def _convert_time_units(values, source_unit, target_unit):
    """Convert time values between arbitrary supported units."""
    validate_time_unit(source_unit)
    validate_time_unit(target_unit)
    values = np.asarray(values, dtype=float)
    if source_unit == target_unit:
        return values

    if source_unit == 'ps':
        return convert_time_from_ps(values, target_unit)

    values_ps = values * TIME_UNITS[source_unit]
    return convert_time_from_ps(values_ps, target_unit)


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


def _format_hbond_id_label(row, atom_labels_by_index):
    """Format one count_by_ids row as donor-hydrogen-acceptor label."""
    def _resolve_atom_label(atom_id):
        atom_id = int(atom_id)
        return atom_labels_by_index.get(atom_id, str(atom_id))

    donor = _resolve_atom_label(row[0])
    hydrogen = _resolve_atom_label(row[1])
    acceptor = _resolve_atom_label(row[2])
    return f"D:{donor}-H:{hydrogen}-A:{acceptor}"


def _normalize_hbonds_base_name(pickle_path):
    """Return the stable replicate-independent base name for a H-bonds pickle."""
    base_name = os.path.splitext(os.path.basename(pickle_path))[0]
    match = re.match(r'^(?P<prefix>.+)_rep\d+$', base_name)
    return match.group('prefix') if match else base_name


def _group_hbonds_pickles(pickle_files):
    """Group H-bonds pickles by replicate-independent base name."""
    grouped = defaultdict(list)
    for pickle_file in pickle_files:
        grouped[_normalize_hbonds_base_name(pickle_file)].append(pickle_file)
    return dict(sorted(grouped.items()))


def _average_hbonds_times(payloads, target_time_unit):
    """Average time-series counts across replicate payloads."""
    converted_times = []
    count_series = []

    for payload in payloads:
        source_unit = payload.get('time_unit', 'ps')
        converted_times.append(_convert_time_units(payload['times'], source_unit, target_time_unit))
        count_series.append(np.asarray(payload['count_by_time'], dtype=float))

    min_len = min(len(series) for series in count_series)
    if min_len == 0:
        return np.array([]), np.array([]), target_time_unit

    reference_times = converted_times[0][:min_len]
    if any(not np.allclose(reference_times, times[:min_len]) for times in converted_times[1:]):
        print("WARNING: H-bonds replicate time axes differ slightly; averaging by frame index.")

    averaged_counts = np.mean(np.vstack([series[:min_len] for series in count_series]), axis=0)
    return reference_times, averaged_counts, target_time_unit


def _aggregate_hbonds_type_counts(payloads):
    """Aggregate donor/acceptor type counts across replicate payloads."""
    totals = defaultdict(int)
    for payload in payloads:
        rows = np.asarray(payload['count_by_type'])
        for row in rows:
            if len(row) < 3:
                continue
            key = (str(row[0]), str(row[1]))
            totals[key] += int(row[2])

    aggregated = [
        [donor, acceptor, total]
        for (donor, acceptor), total in sorted(totals.items(), key=lambda item: (-item[1], item[0][0], item[0][1]))
    ]
    return np.asarray(aggregated, dtype=object)


def _aggregate_hbonds_id_counts(payloads):
    """Aggregate donor/hydrogen/acceptor ID triplets across replicate payloads."""
    totals = defaultdict(int)
    labels_by_index = {}

    for payload in payloads:
        for atom_id, label in (payload.get('atom_labels_by_index') or {}).items():
            labels_by_index.setdefault(int(atom_id), str(label))

        rows = np.asarray(payload['count_by_ids'])
        for row in rows:
            if len(row) < 4:
                continue
            key = (int(row[0]), int(row[1]), int(row[2]))
            totals[key] += int(row[3])

    aggregated = [
        [donor, hydrogen, acceptor, total]
        for (donor, hydrogen, acceptor), total in sorted(
            totals.items(), key=lambda item: (-item[1], item[0][0], item[0][1], item[0][2])
        )
    ]
    return np.asarray(aggregated, dtype=object), labels_by_index


def plot_hbonds_by_time(times, counts, output_path, dpi=DEFAULT_DPI, time_unit='ps',
                        source_time_unit='ps', label=''):
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

    times = _convert_time_units(times, source_time_unit, time_unit)

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


def plot_hbonds_by_ids(id_counts, output_path, dpi=DEFAULT_DPI, top_n=20, label='', atom_labels_by_index=None):
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

    atom_labels_by_index = atom_labels_by_index or {}

    hb_labels = [_format_hbond_id_label(row, atom_labels_by_index) for row in top_data]
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


def plot_hbonds_average(pickle_files, output_dir='.', dpi=DEFAULT_DPI, plot_types=None,
                        top_n=20, time_unit='ps', label=None):
    """Generate averaged H-bond plots across replicate pickle files."""
    if plot_types is None:
        plot_types = ['time', 'type', 'ids']

    grouped_pickles = _group_hbonds_pickles(pickle_files)
    os.makedirs(output_dir, exist_ok=True)

    for group_name, group_files in grouped_pickles.items():
        payloads = []
        for pickle_file in group_files:
            hbonds = _load_hbonds_payload(pickle_file)
            if not isinstance(hbonds, dict):
                raise ValueError(
                    "Unsupported H-bonds pickle schema. Expected a dict payload "
                    "created by run_hbonds_analysis.py with portable arrays."
                )
            payloads.append(hbonds)

        if not payloads:
            continue

        avg_label = label
        if avg_label is None:
            avg_label = prettify_label(group_name.replace('hbonds_plot_', ''))
        avg_label = prettify_label(f'{avg_label} average')

        base_name = f'{group_name}_avg'

        if 'time' in plot_types:
            times, counts, source_time_unit = _average_hbonds_times(payloads, time_unit)
            out = os.path.join(output_dir, f'{base_name}_by_time.png')
            plot_hbonds_by_time(
                times,
                counts,
                out,
                dpi=dpi,
                time_unit=time_unit,
                source_time_unit=source_time_unit,
                label=avg_label,
            )

        if 'type' in plot_types:
            out = os.path.join(output_dir, f'{base_name}_by_type.png')
            plot_hbonds_by_type(
                _aggregate_hbonds_type_counts(payloads),
                out,
                dpi=dpi,
                label=avg_label,
            )

        if 'ids' in plot_types:
            out = os.path.join(output_dir, f'{base_name}_by_ids.png')
            id_counts, atom_labels_by_index = _aggregate_hbonds_id_counts(payloads)
            plot_hbonds_by_ids(
                id_counts,
                out,
                dpi=dpi,
                top_n=top_n,
                label=avg_label,
                atom_labels_by_index=atom_labels_by_index,
            )


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
    source_time_unit = hbonds.get('time_unit', 'ps')
    atom_labels_by_index = hbonds.get('atom_labels_by_index', {})
    hbonds_preset = hbonds.get('hbonds_preset')

    if label is None:
        label = hbonds_preset or os.path.splitext(os.path.basename(pickle_file))[0].replace('hbonds_plot_', '')
    label = prettify_label(label)

    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(pickle_file))[0]

    if 'time' in plot_types:
        out = os.path.join(output_dir, f'{base_name}_by_time.png')
        plot_hbonds_by_time(
            times,
            count_by_time,
            out,
            dpi=dpi,
            time_unit=time_unit,
            source_time_unit=source_time_unit,
            label=label,
        )

    if 'type' in plot_types:
        out = os.path.join(output_dir, f'{base_name}_by_type.png')
        plot_hbonds_by_type(count_by_type, out, dpi=dpi, label=label)

    if 'ids' in plot_types:
        out = os.path.join(output_dir, f'{base_name}_by_ids.png')
        plot_hbonds_by_ids(
            count_by_ids,
            out,
            dpi=dpi,
            top_n=top_n,
            label=label,
            atom_labels_by_index=atom_labels_by_index,
        )


def main():
    """Main function to parse arguments and generate H-bond plots."""
    parser = argparse.ArgumentParser(description='Generate Hydrogen Bond plots from pickle files')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pickle-file', type=str,
                       help='Path to a single H-bonds pickle file')
    group.add_argument('--pickle-files', type=str, nargs='+',
                       help='Paths to multiple H-bonds pickle files to average by replicate group')
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

    if args.pickle_files:
        missing = [path for path in args.pickle_files if not os.path.exists(path)]
        if missing:
            print(f"Error: Pickle file not found: {missing[0]}")
            return 1

        plot_hbonds_average(
            pickle_files=args.pickle_files,
            output_dir=args.output_dir,
            dpi=args.dpi,
            plot_types=args.plot_types,
            top_n=args.top_n,
            time_unit=args.time_unit,
            label=args.label,
        )
    else:
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

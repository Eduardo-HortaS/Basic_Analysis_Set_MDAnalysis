#!/usr/bin/env python3
"""Generate RMSD/RMSF comparison-group plots from existing pickle files.

This helper is used by Nextflow to reproduce the same plot-group behavior as
executor.py for RMSD and RMSF.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
import sys
from typing import Any


def _import_rmsd_plotters():
    """Import RMSD plot helpers when run as package or standalone script."""
    try:
        from plotting.plot_rmsd import plot_rmsd_comparison, plot_rmsd_comparison_average
        return plot_rmsd_comparison, plot_rmsd_comparison_average
    except ModuleNotFoundError:
        # Nextflow runs this script from a work directory, so make the
        # local plotting directory importable as a fallback.
        script_dir = os.path.dirname(os.path.abspath(__file__))
        if script_dir not in sys.path:
            sys.path.insert(0, script_dir)
        from plot_rmsd import plot_rmsd_comparison, plot_rmsd_comparison_average
        return plot_rmsd_comparison, plot_rmsd_comparison_average


def _import_rmsf_plotters():
    """Import RMSF plot helpers when run as package or standalone script."""
    try:
        from plotting.plot_rmsf import plot_rmsf_comparison, plot_rmsf_comparison_average
        return plot_rmsf_comparison, plot_rmsf_comparison_average
    except ModuleNotFoundError:
        # Nextflow runs this script from a work directory, so make the
        # local plotting directory importable as a fallback.
        script_dir = os.path.dirname(os.path.abspath(__file__))
        if script_dir not in sys.path:
            sys.path.insert(0, script_dir)
        from plot_rmsf import plot_rmsf_comparison, plot_rmsf_comparison_average
        return plot_rmsf_comparison, plot_rmsf_comparison_average


def _parse_bool(value: Any) -> bool:
    """Parse common boolean spellings used by CLI/YAML pipelines."""
    if isinstance(value, bool):
        return value
    val = str(value).strip().lower()
    if val in ('1', 'true', 'yes', 'y', 'on'):
        return True
    if val in ('0', 'false', 'no', 'n', 'off'):
        return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def _parse_plot_groups(raw_value: str) -> dict[str, list[tuple[str, str]]]:
    """Parse plot_groups from JSON text or a JSON file path."""
    value = raw_value.strip()
    if not value:
        return {}

    if os.path.exists(value):
        with open(value, 'r', encoding='utf-8') as f:
            payload = json.load(f)
    else:
        payload = json.loads(value)

    if not isinstance(payload, dict):
        raise ValueError('plot_groups must be a JSON object mapping names to member lists')

    parsed: dict[str, list[tuple[str, str]]] = {}
    for group_name, members in payload.items():
        if not isinstance(members, list):
            print(f"WARNING: plot group '{group_name}' is not a list. Skipping.")
            continue

        validated: list[tuple[str, str]] = []
        for member in members:
            if not isinstance(member, (list, tuple)) or len(member) != 2:
                print(f"WARNING: plot group '{group_name}' has invalid member {member!r}. Skipping member.")
                continue
            system, variation = member
            validated.append((str(system), str(variation)))

        if validated:
            parsed[str(group_name)] = validated

    return parsed


def _normalize_replicate_mode(mode: str) -> str:
    mode_norm = (mode or 'separate').strip().lower()
    if mode_norm not in ('separate', 'average'):
        print(f"WARNING: Unknown replicate_mode '{mode}'. Falling back to 'separate'.")
        return 'separate'
    return mode_norm


def _find_pickle_for_member(work_dir: str, prefix: str, system: str, variation: str,
                            rep: int, sel_idx: int | None = None,
                            ref_idx: int | None = None) -> str | None:
    """Find a pickle matching system/variation/rep with optional ref/sel suffixes."""
    ref_part = f'_ref{ref_idx}' if ref_idx is not None else ''
    if sel_idx is not None:
        pattern = os.path.join(work_dir, f'{prefix}_{system}_{variation}_rep{rep}{ref_part}_sel{sel_idx}.pkl')
    else:
        pattern = os.path.join(work_dir, f'{prefix}_{system}_{variation}_rep{rep}{ref_part}.pkl')
    matches = glob.glob(pattern)
    if matches:
        return matches[0]

    pattern = os.path.join(work_dir, f'{prefix}_{system}_{variation}_rep{rep}{ref_part}*.pkl')
    matches = glob.glob(pattern)
    return matches[0] if matches else None


def _detect_selection_indices(work_dir: str, prefix: str) -> list[int | None]:
    sel_indices: set[int] = set()
    for match in glob.glob(os.path.join(work_dir, f'{prefix}_*_sel*.pkl')):
        base = os.path.basename(match)
        marker = '_sel'
        pos = base.rfind(marker)
        if pos == -1:
            continue
        tail = base[pos + len(marker):]
        num = []
        for ch in tail:
            if ch.isdigit():
                num.append(ch)
            else:
                break
        if num:
            sel_indices.add(int(''.join(num)))
    return sorted(sel_indices) if sel_indices else [None]


def _detect_ref_indices(work_dir: str, prefix: str) -> list[int | None]:
    ref_indices: set[int] = set()
    for match in glob.glob(os.path.join(work_dir, f'{prefix}_*_ref*.pkl')):
        base = os.path.basename(match)
        marker = '_ref'
        pos = base.rfind(marker)
        if pos == -1:
            continue
        tail = base[pos + len(marker):]
        num = []
        for ch in tail:
            if ch.isdigit():
                num.append(ch)
            else:
                break
        if num:
            ref_indices.add(int(''.join(num)))
    return sorted(ref_indices) if ref_indices else [None]


def _get_selection_label_from_pickle(pkl_path: str) -> tuple[str, str]:
    try:
        with open(pkl_path, 'rb') as f:
            data = pickle.load(f)
        if isinstance(data, dict):
            return data.get('selection', ''), data.get('ref_selection', '')
    except Exception:
        pass
    return '', ''


def _plot_rmsd_groups(plot_groups: dict[str, list[tuple[str, str]]], *, work_dir: str,
                      output_dir: str, num_replicates: int, dpi: int,
                      replicate_mode: str, rmsd_show_kde: bool) -> int:
    plot_rmsd_comparison, plot_rmsd_comparison_average = _import_rmsd_plotters()

    produced = 0
    sel_indices = _detect_selection_indices(work_dir, 'rmsd_plot')
    ref_indices = _detect_ref_indices(work_dir, 'rmsd_plot')

    for group_name, members in plot_groups.items():
        for ref_idx in ref_indices:
            ref_suffix = f' (ref{ref_idx})' if ref_idx is not None else ''
            for sel_idx in sel_indices:
                sel_suffix = f' (sel{sel_idx})' if sel_idx is not None else ''

                if replicate_mode == 'average':
                    pickle_groups = []
                    labels = []
                    selection_label = ''
                    ref_selection_label = ''

                    for system, variation in members:
                        rep_pickles = []
                        for rep in range(1, num_replicates + 1):
                            pkl = _find_pickle_for_member(
                                work_dir, 'rmsd_plot', system, variation, rep, sel_idx, ref_idx)
                            if pkl and os.path.exists(pkl):
                                rep_pickles.append(pkl)
                                if not selection_label:
                                    selection_label, ref_selection_label = _get_selection_label_from_pickle(pkl)
                        if rep_pickles:
                            pickle_groups.append(rep_pickles)
                            labels.append(f'{system} / {variation}')

                    if not pickle_groups:
                        print(f"WARNING: No RMSD pickles found for group '{group_name}'{ref_suffix}{sel_suffix}. Skipping.")
                        continue

                    gname = f'{group_name}{ref_suffix}{sel_suffix}'.strip()
                    subtitle = selection_label
                    if ref_selection_label and ref_selection_label != selection_label:
                        subtitle = f'{selection_label} (aligned on: {ref_selection_label})'

                    plot_rmsd_comparison_average(
                        pickle_groups,
                        labels,
                        gname,
                        output_dir=output_dir,
                        dpi=dpi,
                        show_kde=rmsd_show_kde,
                        selection_label=subtitle,
                    )
                    produced += 1
                else:
                    plotted_for_group = False
                    for rep in range(1, num_replicates + 1):
                        pkl_files = []
                        labels = []
                        selection_label = ''
                        ref_selection_label = ''

                        for system, variation in members:
                            pkl = _find_pickle_for_member(
                                work_dir, 'rmsd_plot', system, variation, rep, sel_idx, ref_idx)
                            if pkl and os.path.exists(pkl):
                                pkl_files.append(pkl)
                                labels.append(f'{system} / {variation}')
                                if not selection_label:
                                    selection_label, ref_selection_label = _get_selection_label_from_pickle(pkl)

                        if not pkl_files:
                            continue

                        rep_suffix = f'_rep{rep}' if num_replicates > 1 else ''
                        gname = f'{group_name}{ref_suffix}{sel_suffix}{rep_suffix}'.strip()
                        subtitle = selection_label
                        if ref_selection_label and ref_selection_label != selection_label:
                            subtitle = f'{selection_label} (aligned on: {ref_selection_label})'

                        plot_rmsd_comparison(
                            pkl_files,
                            labels,
                            gname,
                            output_dir=output_dir,
                            dpi=dpi,
                            show_kde=rmsd_show_kde,
                            selection_label=subtitle,
                        )
                        produced += 1
                        plotted_for_group = True

                    if not plotted_for_group:
                        print(f"WARNING: No RMSD pickles found for group '{group_name}'{ref_suffix}{sel_suffix}. Skipping.")

    return produced


def _plot_rmsf_groups(plot_groups: dict[str, list[tuple[str, str]]], *, work_dir: str,
                      output_dir: str, num_replicates: int, dpi: int,
                      replicate_mode: str, target_selection: str) -> int:
    plot_rmsf_comparison, plot_rmsf_comparison_average = _import_rmsf_plotters()

    produced = 0

    # Detect both selection indices and chain IDs in pickles
    sel_indices = _detect_selection_indices(work_dir, 'rmsf_plot')
    
    chain_ids: set[str] = set()
    for pkl in glob.glob(os.path.join(work_dir, 'rmsf_plot_*.pkl')):
        base = os.path.basename(pkl)
        marker = '_chain'
        pos = base.rfind(marker)
        if pos == -1:
            continue
        chain = base[pos + len(marker):].removesuffix('.pkl')
        if chain:
            chain_ids.add(chain)

    chain_id_list = sorted(chain_ids) if chain_ids else [None]

    for group_name, members in plot_groups.items():
        for sel_idx in sel_indices:
            sel_suffix = f' (sel{sel_idx})' if sel_idx is not None else ''
            for chain_id in chain_id_list:
                chain_suffix = f'_chain{chain_id}' if chain_id else ''

                if replicate_mode == 'average':
                    pickle_groups = []
                    labels = []

                    for system, variation in members:
                        rep_pickles = []
                        for rep in range(1, num_replicates + 1):
                            pkl = _find_pickle_for_rmsf_member(
                                work_dir, system, variation, rep, sel_idx, chain_id)
                            if pkl and os.path.exists(pkl):
                                rep_pickles.append(pkl)
                        if rep_pickles:
                            pickle_groups.append(rep_pickles)
                            labels.append(f'{system} / {variation}')

                    if not pickle_groups:
                        continue

                    gname = f'{group_name}{sel_suffix}{chain_suffix}'
                    plot_rmsf_comparison_average(
                        pickle_groups,
                        labels,
                        gname,
                        output_dir=output_dir,
                        dpi=dpi,
                        selection_label=target_selection,
                    )
                    produced += 1
                else:
                    plotted_for_group = False
                    for rep in range(1, num_replicates + 1):
                        pkl_files = []
                        labels = []

                        for system, variation in members:
                            pkl = _find_pickle_for_rmsf_member(
                                work_dir, system, variation, rep, sel_idx, chain_id)
                            if pkl and os.path.exists(pkl):
                                pkl_files.append(pkl)
                                labels.append(f'{system} / {variation}')

                        if not pkl_files:
                            continue

                        rep_suffix = f'_rep{rep}' if num_replicates > 1 else ''
                        gname = f'{group_name}{sel_suffix}{chain_suffix}{rep_suffix}'
                        plot_rmsf_comparison(
                            pkl_files,
                            labels,
                            gname,
                            output_dir=output_dir,
                            dpi=dpi,
                            selection_label=target_selection,
                        )
                        produced += 1
                        plotted_for_group = True

                    if not plotted_for_group:
                        combined_suffix = f'{sel_suffix}{chain_suffix}'.strip()
                        if combined_suffix:
                            print(f"WARNING: No RMSF pickles found for group '{group_name}'{combined_suffix}. Skipping.")

    return produced


def _find_pickle_for_rmsf_member(work_dir: str, system: str, variation: str,
                                 rep: int, sel_idx: int | None = None,
                                 chain_id: str | None = None) -> str | None:
    """Find an RMSF pickle matching system/variation/rep with optional sel/chain suffixes."""
    sel_part = f'_sel{sel_idx}' if sel_idx is not None else ''
    chain_part = f'_chain{chain_id}' if chain_id is not None else ''
    
    pattern = os.path.join(work_dir, f'rmsf_plot_{system}_{variation}_rep{rep}{sel_part}{chain_part}.pkl')
    matches = glob.glob(pattern)
    if matches:
        return matches[0]
    return None


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Generate grouped RMSD/RMSF comparison plots from pickles.')
    parser.add_argument('--analysis', required=True, choices=('rmsd', 'rmsf'))
    parser.add_argument('--plot-groups', required=True,
                        help='JSON object or JSON file path with group definitions')
    parser.add_argument('--work-dir', default='.')
    parser.add_argument('--output-dir', default='.')
    parser.add_argument('--num-replicates', type=int, default=1)
    parser.add_argument('--replicate-mode', default='separate')
    parser.add_argument('--dpi', type=int, default=400)
    parser.add_argument('--target-selection', default='protein and backbone')
    parser.add_argument('--rmsd-show-kde', default='true')
    args = parser.parse_args()

    plot_groups = _parse_plot_groups(args.plot_groups)
    if not plot_groups:
        print('No valid plot_groups provided. Nothing to do.')
        return

    replicate_mode = _normalize_replicate_mode(args.replicate_mode)
    rmsd_show_kde = _parse_bool(args.rmsd_show_kde)
    os.makedirs(args.output_dir, exist_ok=True)

    if args.analysis == 'rmsd':
        count = _plot_rmsd_groups(
            plot_groups,
            work_dir=args.work_dir,
            output_dir=args.output_dir,
            num_replicates=args.num_replicates,
            dpi=args.dpi,
            replicate_mode=replicate_mode,
            rmsd_show_kde=rmsd_show_kde,
        )
    else:
        count = _plot_rmsf_groups(
            plot_groups,
            work_dir=args.work_dir,
            output_dir=args.output_dir,
            num_replicates=args.num_replicates,
            dpi=args.dpi,
            replicate_mode=replicate_mode,
            target_selection=args.target_selection,
        )

    print(f'Generated {count} grouped {args.analysis.upper()} comparison plot(s).')


if __name__ == '__main__':
    main()

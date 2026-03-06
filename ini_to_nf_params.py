#!/usr/bin/env python3
"""Convert an executor.py INI config file into a Nextflow ``-params-file`` YAML.

Usage
-----
    python ini_to_nf_params.py my_config.ini              # prints YAML to stdout
    python ini_to_nf_params.py my_config.ini -o params.yml # writes to file

Then run Nextflow with::

    nextflow run main.nf -params-file params.yml

This ensures that executor.py and the Nextflow pipeline receive identical
parameters, which is a prerequisite for output-parity testing.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

# Allow importing executor.py from the same directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from executor import load_config  # noqa: E402


def _cfg_to_nf_params(cfg: dict) -> dict:
    """Map the executor.py *cfg* dict to Nextflow param names."""
    params: dict = {}

    # Systems  ────────────────────────────────────────────────────────────
    params['systems'] = json.dumps(cfg['systems'])
    params['variations'] = json.dumps(cfg['variations'])
    params['num_replicates'] = cfg['num_replicates']
    params['traj_format'] = cfg['traj_format']
    params['top_format'] = cfg['top_format']
    params['start_frame'] = cfg['start_frame']

    # Paths  ──────────────────────────────────────────────────────────────
    params['input_dir'] = os.path.abspath(cfg['input_dir'])
    params['outdir'] = os.path.abspath(cfg['output_dir'])

    # Analysis toggles  ──────────────────────────────────────────────────
    for key in ('run_rmsd', 'run_rmsf', 'run_2d_rmsd', 'run_rog', 'run_hbonds'):
        params[key] = cfg[key]

    # Plot toggles  ──────────────────────────────────────────────────────
    # Respect the ``run_plots`` master switch: when it is False in the
    # executor config, ALL per-analysis plot flags must be False so that
    # the Nextflow pipeline behaves identically (executor.py skips all
    # plotting when run_plots=false regardless of individual toggles).
    _plots_enabled = cfg.get('run_plots', True)
    for key in ('plot_rmsd', 'plot_rmsf', 'plot_2d_rmsd', 'plot_rog', 'plot_hbonds'):
        params[key] = cfg[key] if _plots_enabled else False

    # Selections  ────────────────────────────────────────────────────────
    params['target_selection'] = cfg['target_selection']
    # NF receives a single ref_selection; for multi-ref the first is used.
    params['ref_selection'] = cfg['ref_selection'][0] if cfg['ref_selection'] else 'protein and backbone'
    params['rog_selection'] = cfg['rog_selection']

    # PBC wrapping
    params['wrap_selection'] = cfg['wrap_selection'] if cfg['wrap_selection'] else 'none'

    # RMSD  ──────────────────────────────────────────────────────────────
    if cfg['time_interval_between_frames'] is not None:
        params['time_interval_between_frames'] = cfg['time_interval_between_frames']
    params['time_unit'] = cfg['time_unit']
    if cfg['group_selections']:
        params['group_selections'] = json.dumps(cfg['group_selections'])

    # RMSF  ──────────────────────────────────────────────────────────────
    if cfg['chain_intervals']:
        params['chain_intervals'] = json.dumps(cfg['chain_intervals'])

    # H-bonds  ───────────────────────────────────────────────────────────
    params['d_a_cutoff'] = cfg['d_a_cutoff']
    params['d_h_a_angle_cutoff'] = cfg['d_h_a_angle_cutoff']
    params['update_selections'] = cfg['update_selections']
    if cfg['acceptors_sel']:
        params['acceptors_sel'] = cfg['acceptors_sel']
    if cfg['hydrogens_sel']:
        params['hydrogens_sel'] = cfg['hydrogens_sel']
    if cfg['between_pairs']:
        params['between_pairs'] = json.dumps(cfg['between_pairs'])

    # Plotting  ──────────────────────────────────────────────────────────
    params['dpi'] = cfg['dpi']
    params['hbonds_top_n'] = cfg.get('hbonds_top_n', 20)

    # Parallelization  ───────────────────────────────────────────────────
    params['parallel_backend'] = cfg.get('parallel_backend', 'serial')
    _nw = cfg.get('n_workers')
    if _nw is not None:
        params['n_workers'] = _nw

    return params


def _dump_yaml(params: dict) -> str:
    """Minimal YAML serialisation (no PyYAML dependency required)."""
    lines: list[str] = []
    for key, value in params.items():
        if isinstance(value, bool):
            lines.append(f"{key}: {'true' if value else 'false'}")
        elif isinstance(value, (int, float)):
            lines.append(f"{key}: {value}")
        elif isinstance(value, str):
            # Quote strings that may confuse YAML parsers
            if any(c in value for c in ':[]{},&*?|>!%@`"\' \t#') or value in ('true', 'false', 'null', 'yes', 'no'):
                escaped = value.replace("'", "''")
                lines.append(f"{key}: '{escaped}'")
            else:
                lines.append(f"{key}: {value}")
        elif value is None:
            lines.append(f"{key}: null")
        else:
            lines.append(f"{key}: '{value}'")
    return '\n'.join(lines) + '\n'


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Convert an executor.py INI config to a Nextflow -params-file YAML.')
    parser.add_argument('config', help='Path to the .INI configuration file')
    parser.add_argument('-o', '--output', default=None,
                        help='Output YAML path (default: stdout)')
    args = parser.parse_args()

    cfg = load_config(args.config)
    params = _cfg_to_nf_params(cfg)
    yaml_str = _dump_yaml(params)

    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(yaml_str)
        print(f"Wrote Nextflow params to {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(yaml_str)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
executor.py — Run the full MD analysis pipeline from a .INI configuration file.

Reads an INI file that specifies systems, variations, analysis toggles,
selections, paths, and plotting options, then executes every enabled
analysis followed by its corresponding plots.

Usage:
    python executor.py config.ini
    python executor.py config.ini --dry-run
    python executor.py config.ini --force
"""

import os
import sys
import glob
import json
import argparse
import configparser
import time as _time

# ─── Ensure repo root is importable ──────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from utils import validate_time_unit, SUPPORTED_TIME_UNITS, resolve_trajectory_file, cleanup_work_directory
from parallelization import VALID_BACKENDS


# ─── INI parsing helpers ─────────────────────────────────────────────────────

def _parse_json_or_path(value):
    """Parse a value as JSON, or read it from a file path if it points to one."""
    value = value.strip()
    if os.path.isfile(value):
        with open(value, 'r', encoding='utf-8') as f:
            return json.load(f)
    return json.loads(value)


def _parse_bool(value):
    return value.strip().lower() in ('true', 'yes', '1', 'on')


def _parse_optional_str(value):
    """Return None for empty / 'none' values, else the stripped string."""
    v = value.strip()
    if not v or v.lower() == 'none':
        return None
    return v


def _parse_optional_float(value):
    v = _parse_optional_str(value)
    return float(v) if v is not None else None


def _parse_optional_int(value):
    v = _parse_optional_str(value)
    return int(v) if v is not None else None


def _parse_optional_json(value):
    v = _parse_optional_str(value)
    if v is None:
        return None
    return _parse_json_or_path(v)


def _normalize_chain_intervals(raw, systems):
    """
    Normalize chain_intervals into a per-system dict.

    Accepts two formats:
      - Per-system (new):  {"Ung_G-C_4": {"A": [1, 239]}, "Mug_G-C_4": {"A": [1, 211]}, ...}
      - Legacy flat:       {"A": [1, 120], "B": [121, 239]}

    Detection: if every value is a dict → per-system; if every value is a list → legacy flat.
    Legacy flat dicts are broadcast to every system in `systems`.

    Returns None when raw is None.
    """
    if raw is None:
        return None

    if not isinstance(raw, dict) or not raw:
        raise ValueError(
            f"chain_intervals must be a JSON dict, got {type(raw).__name__}: {raw!r}"
        )

    first_val = next(iter(raw.values()))

    if isinstance(first_val, dict):
        # Per-system format — validate that all values are dicts
        for key, val in raw.items():
            if not isinstance(val, dict):
                raise ValueError(
                    f"chain_intervals: mixed format detected. Key '{key}' has "
                    f"value type {type(val).__name__}, expected dict (per-system format)."
                )
        return raw

    if isinstance(first_val, list):
        # Legacy flat format — broadcast to every system
        for key, val in raw.items():
            if not isinstance(val, list):
                raise ValueError(
                    f"chain_intervals: mixed format detected. Key '{key}' has "
                    f"value type {type(val).__name__}, expected list (legacy flat format)."
                )
        print(f"  NOTE: Broadcasting legacy flat chain_intervals to all {len(systems)} systems.")
        return {system: json.loads(json.dumps(raw)) for system in systems}

    raise ValueError(
        f"chain_intervals: unrecognized format. First value is {type(first_val).__name__}, "
        f"expected dict (per-system) or list (legacy flat)."
    )


# ─── Config loading ──────────────────────────────────────────────────────────

def load_config(ini_path):
    """
    Load and validate an INI config file, returning a structured dict.

    Required sections: [systems], [analysis]
    Optional sections: [paths], [selections], [rmsd], [rmsf], [hbonds], [plotting]
    """
    if not os.path.isfile(ini_path):
        print(f"ERROR: Config file not found: {ini_path}")
        sys.exit(1)

    cp = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    cp.read(ini_path, encoding='utf-8')

    cfg = {}

    # ── [systems] — required ──────────────────────────────────────────────
    if not cp.has_section('systems'):
        print("ERROR: [systems] section is required in the INI file.")
        sys.exit(1)

    cfg['systems'] = _parse_json_or_path(cp.get('systems', 'systems'))
    cfg['variations'] = _parse_json_or_path(cp.get('systems', 'variations'))
    cfg['num_replicates'] = cp.getint('systems', 'num_replicates', fallback=3)
    cfg['traj_format'] = cp.get('systems', 'traj_format', fallback='dcd').strip()
    cfg['top_format'] = cp.get('systems', 'top_format', fallback='top').strip()
    cfg['start_frame'] = cp.getint('systems', 'start_frame', fallback=500)

    # ── [paths] — optional, defaults to cwd ───────────────────────────────
    cfg['input_dir'] = cp.get('paths', 'input_dir', fallback='.').strip()
    cfg['output_dir'] = cp.get('paths', 'output_dir', fallback='./results').strip()
    cfg['plot_dir'] = cp.get('paths', 'plot_dir', fallback='').strip()
    if not cfg['plot_dir']:
        cfg['plot_dir'] = os.path.join(cfg['output_dir'], 'plots')

    # ── [analysis] — required (which analyses to run) ─────────────────────
    if not cp.has_section('analysis'):
        print("ERROR: [analysis] section is required in the INI file.")
        sys.exit(1)

    cfg['run_rmsd'] = _parse_bool(cp.get('analysis', 'run_rmsd', fallback='false'))
    cfg['run_rmsf'] = _parse_bool(cp.get('analysis', 'run_rmsf', fallback='false'))
    cfg['run_2d_rmsd'] = _parse_bool(cp.get('analysis', 'run_2d_rmsd', fallback='false'))
    cfg['run_rog'] = _parse_bool(cp.get('analysis', 'run_rog', fallback='false'))
    cfg['run_hbonds'] = _parse_bool(cp.get('analysis', 'run_hbonds', fallback='false'))
    cfg['run_plots'] = _parse_bool(cp.get('analysis', 'run_plots', fallback='true'))

    # Per-analysis plot toggles (default to true for backward compatibility)
    cfg['plot_rmsd'] = _parse_bool(cp.get('analysis', 'plot_rmsd', fallback='true'))
    cfg['plot_rmsf'] = _parse_bool(cp.get('analysis', 'plot_rmsf', fallback='true'))
    cfg['plot_2d_rmsd'] = _parse_bool(cp.get('analysis', 'plot_2d_rmsd', fallback='true'))
    cfg['plot_rog'] = _parse_bool(cp.get('analysis', 'plot_rog', fallback='true'))
    cfg['plot_hbonds'] = _parse_bool(cp.get('analysis', 'plot_hbonds', fallback='true'))

    # Global time unit — try [analysis] first, fall back to [rmsd] for backward compat
    _time_unit = cp.get('analysis', 'time_unit', fallback='').strip()
    if not _time_unit:
        _time_unit = cp.get('rmsd', 'time_unit', fallback='ns').strip()
    cfg['time_unit'] = _time_unit
    validate_time_unit(cfg['time_unit'])

    # ── [selections] — shared selections ──────────────────────────────────
    cfg['target_selection'] = cp.get('selections', 'target_selection',
                                     fallback='protein and backbone').strip()

    # ref_selection can be a JSON list or a plain string (backward compat).
    _raw_ref_sel = cp.get('selections', 'ref_selection',
                          fallback='protein and backbone').strip()
    try:
        _parsed_ref_sel = json.loads(_raw_ref_sel)
        if isinstance(_parsed_ref_sel, list):
            cfg['ref_selection'] = [str(s).strip() for s in _parsed_ref_sel]
        else:
            cfg['ref_selection'] = [str(_parsed_ref_sel).strip()]
    except (json.JSONDecodeError, TypeError):
        cfg['ref_selection'] = [_raw_ref_sel]

    cfg['rog_selection'] = cp.get('selections', 'rog_selection',
                                  fallback='protein and backbone').strip()

    # wrap_selection — controls PBC wrapping behaviour.
    # 'auto' (default): unwrap/centre protein, wrap everything else.
    # 'none' / None:     disable all PBC wrapping.
    # any other string:  MDAnalysis selection for atoms to wrap.
    _raw_wrap_sel = cp.get('selections', 'wrap_selection', fallback='auto').strip()
    if not _raw_wrap_sel or _raw_wrap_sel.lower() == 'none':
        cfg['wrap_selection'] = None
    else:
        cfg['wrap_selection'] = _raw_wrap_sel

    # ── [rmsd] — RMSD-specific ────────────────────────────────────────────
    cfg['time_interval_between_frames'] = _parse_optional_float(
        cp.get('rmsd', 'time_interval_between_frames', fallback='none'))

    # time_interval_between_frames is mandatory when RMSD is enabled — it
    # cannot be inferred from trajectory files (especially DCD).
    if cfg['run_rmsd'] and cfg['time_interval_between_frames'] is None:
        print(
            "ERROR: 'time_interval_between_frames' is not set or is set to 'none' "
            "in the [rmsd] section of your INI file.\n"
            "  This value CANNOT be inferred from trajectory files and must be "
            "provided explicitly (in picoseconds).\n"
            "  Without it the time axis of RMSD (and any DCD-based correction) "
            "will be unreliable.\n"
            "  Please set  time_interval_between_frames = <value>  in [rmsd] "
            "and re-run."
        )
        sys.exit(1)

    _raw_group_selections = cp.get('rmsd', 'group_selections', fallback='none')
    try:
        cfg['group_selections'] = _parse_optional_json(_raw_group_selections)
    except json.JSONDecodeError as e:
        raw = (_raw_group_selections or '').strip()
        # Support legacy syntax like: [chainid A, chainid B, chainid C]
        if raw.startswith('[') and raw.endswith(']'):
            inner = raw[1:-1].strip()
            if inner:
                parsed = []
                for token in inner.split(','):
                    item = token.strip().strip('"\'')
                    if item:
                        parsed.append(item)
                cfg['group_selections'] = parsed if parsed else None
            else:
                cfg['group_selections'] = None
        else:
            print(
                "ERROR: Invalid JSON for [rmsd] group_selections.\n"
                "  Use a JSON list, e.g.: [\"segid A\", \"segid B\", \"segid C\"]\n"
                f"  Parser error: {e}"
            )
            sys.exit(1)

    if cfg['group_selections'] is not None:
        if not isinstance(cfg['group_selections'], list):
            print(
                "ERROR: [rmsd] group_selections must be a JSON list of selection strings.\n"
                f"  Got: {type(cfg['group_selections']).__name__}"
            )
            sys.exit(1)
        cfg['group_selections'] = [str(s).strip() for s in cfg['group_selections'] if str(s).strip()]
        if not cfg['group_selections']:
            cfg['group_selections'] = None

    # ── [rmsf] — RMSF-specific ────────────────────────────────────────────
    raw_chain_intervals = _parse_optional_json(
        cp.get('rmsf', 'chain_intervals', fallback='none'))
    cfg['chain_intervals'] = _normalize_chain_intervals(
        raw_chain_intervals, cfg['systems'])

    cfg['rmsf_group_selections'] = _parse_optional_json(
        cp.get('rmsf', 'group_selections', fallback='none'))

    if cfg['rmsf_group_selections'] is not None:
        if not isinstance(cfg['rmsf_group_selections'], list):
            print(
                "ERROR: [rmsf] group_selections must be a JSON list of selection strings.\n"
                f"  Got: {type(cfg['rmsf_group_selections']).__name__}"
            )
            sys.exit(1)
        cfg['rmsf_group_selections'] = [str(s).strip() for s in cfg['rmsf_group_selections'] if str(s).strip()]
        if not cfg['rmsf_group_selections']:
            cfg['rmsf_group_selections'] = None

    # ── [hbonds] — H-bond specific ────────────────────────────────────────
    cfg['d_a_cutoff'] = float(cp.get('hbonds', 'd_a_cutoff', fallback='3.5'))
    cfg['d_h_a_angle_cutoff'] = float(cp.get('hbonds', 'd_h_a_angle_cutoff', fallback='150.0'))
    cfg['update_selections'] = _parse_bool(cp.get('hbonds', 'update_selections', fallback='true'))
    cfg['acceptors_sel'] = _parse_optional_str(cp.get('hbonds', 'acceptors_sel', fallback='none'))
    cfg['hydrogens_sel'] = _parse_optional_str(cp.get('hbonds', 'hydrogens_sel', fallback='none'))
    cfg['between_pairs'] = _parse_optional_json(cp.get('hbonds', 'between_pairs', fallback='none'))
    if cfg['between_pairs'] is not None:
        if not isinstance(cfg['between_pairs'], list):
            print(
                "ERROR: [hbonds] between_pairs must be a JSON list of [selection1, selection2] pairs."
            )
            sys.exit(1)
        cleaned_pairs = []
        for pair in cfg['between_pairs']:
            if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                print(f"ERROR: Invalid between_pairs entry: {pair!r}. Expected [selection1, selection2].")
                sys.exit(1)
            a, b = str(pair[0]).strip(), str(pair[1]).strip()
            if not a or not b:
                print(f"ERROR: between_pairs entries cannot be empty: {pair!r}")
                sys.exit(1)
            cleaned_pairs.append([a, b])
        cfg['between_pairs'] = cleaned_pairs

    # ── [parallelization] — parallel backend settings ───────────────────
    _parallel_backend = cp.get('parallelization', 'backend', fallback='serial').strip().lower()
    if _parallel_backend not in VALID_BACKENDS:
        print(f"ERROR: Unknown parallel backend '{_parallel_backend}'. "
              f"Accepted values: {', '.join(VALID_BACKENDS)}")
        sys.exit(1)
    cfg['parallel_backend'] = _parallel_backend
    cfg['n_workers'] = _parse_optional_int(
        cp.get('parallelization', 'n_workers', fallback='none'))

    # ── [plotting] — plot settings ────────────────────────────────────────
    cfg['dpi'] = int(cp.get('plotting', 'dpi', fallback='400'))
    cfg['rmsd_show_kde'] = _parse_bool(cp.get('plotting', 'rmsd_show_kde', fallback='true'))
    cfg['rog_show_kde'] = _parse_bool(cp.get('plotting', 'rog_show_kde', fallback='true'))
    cfg['rmsf_plot_type'] = cp.get('plotting', 'rmsf_plot_type', fallback='line').strip()
    cfg['twod_rmsd_cmap'] = cp.get('plotting', 'twod_rmsd_cmap', fallback='viridis').strip()
    cfg['hbonds_top_n'] = int(cp.get('plotting', 'hbonds_top_n', fallback='20'))

    # ── [plot_groups] — named comparison groups ─────────────────────────
    cfg['plot_groups'] = {}
    cfg['replicate_mode'] = 'separate'  # default
    if cp.has_section('plot_groups'):
        cfg['replicate_mode'] = cp.get('plot_groups', 'replicate_mode', fallback='separate').strip().lower()
        if cfg['replicate_mode'] not in ('separate', 'average'):
            print(f"WARNING: Unknown replicate_mode '{cfg['replicate_mode']}'. Falling back to 'separate'.")
            cfg['replicate_mode'] = 'separate'
        for key in cp.options('plot_groups'):
            if key == 'replicate_mode':
                continue
            raw_val = cp.get('plot_groups', key).strip()
            if not raw_val or raw_val.lower() == 'none':
                continue
            try:
                members = json.loads(raw_val)
                if not isinstance(members, list):
                    print(f"WARNING: plot_groups.{key} must be a JSON list of [system, variation] pairs. Skipping.")
                    continue
                # Validate each member
                validated = []
                for m in members:
                    if not isinstance(m, list) or len(m) != 2:
                        print(f"WARNING: plot_groups.{key}: each entry must be [system, variation]. Got {m!r}. Skipping entry.")
                        continue
                    validated.append(tuple(m))
                if validated:
                    cfg['plot_groups'][key] = validated
            except json.JSONDecodeError as e:
                print(f"WARNING: plot_groups.{key}: invalid JSON: {e}. Skipping.")

    return cfg


# ─── Data directory setup ────────────────────────────────────────────────────

# Per-analysis subdirectory names — must match Nextflow publishDir layout.
_ANALYSIS_SUBDIRS = {
    'rmsd':    'rmsd',
    'rmsf':    'rmsf',
    '2d_rmsd': '2d_rmsd',
    'rog':     'rog',
    'hbonds':  'hbonds',
}


def setup_workdir(cfg):
    """
    Create the working directory structure with symlinks to the input data.
    Returns the work_dir path where analyses will run.

    Pickle files are written into per-analysis subdirectories under work_dir
    (e.g. ``work/rmsd/``, ``work/rmsf/``) so that the layout matches the
    Nextflow ``publishDir`` output.
    """
    work_dir = os.path.abspath(os.path.join(cfg['output_dir'], 'work'))
    os.makedirs(work_dir, exist_ok=True)

    input_dir = os.path.abspath(cfg['input_dir'])

    for system in cfg['systems']:
        for variation in cfg['variations'][system]:
            data_dir = os.path.join(work_dir, system, variation)
            os.makedirs(data_dir, exist_ok=True)

            src_dir = os.path.join(input_dir, system, variation)
            if not os.path.isdir(src_dir):
                print(f"  WARNING: Source directory not found: {src_dir}")
                continue

            for fname in os.listdir(src_dir):
                src = os.path.join(src_dir, fname)
                dst = os.path.join(data_dir, fname)
                if os.path.isfile(src) and not os.path.exists(dst):
                    os.symlink(os.path.abspath(src), dst)

    # Pre-create per-analysis subdirectories
    for subdir in _ANALYSIS_SUBDIRS.values():
        os.makedirs(os.path.join(work_dir, subdir), exist_ok=True)

    return work_dir


def _clean_ephemeral_files(work_dir):
    """Remove aligned-trajectory side-effect files (rmsfit_*) from work_dir.

    Nextflow processes are isolated and never publish these files, so
    executor.py cleans them up after all analyses complete to keep the
    output tree consistent with Nextflow runs.
    """
    removed = []
    for pattern in ('rmsfit_*.dcd', 'rmsfit_*.xtc'):
        for fpath in glob.glob(os.path.join(work_dir, pattern)):
            os.remove(fpath)
            removed.append(os.path.basename(fpath))
    if removed:
        print(f"  Cleaned {len(removed)} ephemeral aligned-trajectory file(s).")
    return removed


# ─── Analysis runners ────────────────────────────────────────────────────────

def run_rmsd(cfg, dry_run=False):
    """Run RMSD analysis for all systems, iterating over ref_selections."""
    from run_rms_analysis import run_rms_analysis

    print(f"\n{'='*60}")
    print("  RMSD ANALYSIS")
    print(f"{'='*60}")

    if dry_run:
        print("  [DRY RUN] Would run RMSD for:", cfg['systems'])
        print("  [DRY RUN] ref_selections:", cfg['ref_selection'])
        return

    ref_selections = cfg['ref_selection']  # already a list
    for ref_idx, ref_sel in enumerate(ref_selections):
        ref_suffix = f'_ref{ref_idx}' if len(ref_selections) > 1 else ''
        if ref_suffix:
            print(f"\n  --- ref_selection {ref_idx}: '{ref_sel}' ---")

        run_rms_analysis(
            systems=cfg['systems'],
            variations=cfg['variations'],
            num_replicates=cfg['num_replicates'],
            analysis='RMSD',
            target_selection=cfg['target_selection'],
            ref_selection=ref_sel,
            start_frame=cfg['start_frame'],
            traj_format=cfg['traj_format'],
            top_format=cfg['top_format'],
            group_selections=cfg['group_selections'],
            time_interval_between_frames=cfg['time_interval_between_frames'],
            time_unit=cfg['time_unit'],
            ref_suffix=ref_suffix,
            wrap_selection=cfg['wrap_selection'],
            output_dir=os.path.join(cfg['_work_dir'], _ANALYSIS_SUBDIRS['rmsd']),
            parallel_backend=cfg['parallel_backend'],
            n_workers=cfg['n_workers'],
            strict_wrapping=cfg['strict_wrapping'],
        )


def run_rmsf(cfg, dry_run=False):
    """Run RMSF analysis for all systems."""
    from run_rms_analysis import run_rms_analysis

    print(f"\n{'='*60}")
    print("  RMSF ANALYSIS")
    print(f"{'='*60}")

    if dry_run:
        print("  [DRY RUN] Would run RMSF for:", cfg['systems'])
        if cfg['chain_intervals']:
            print("  Chain intervals:", cfg['chain_intervals'])
        return

    # RMSF uses only the first ref_selection (no iteration)
    ref_sel = cfg['ref_selection'][0] if cfg['ref_selection'] else 'protein and backbone'

    run_rms_analysis(
        systems=cfg['systems'],
        variations=cfg['variations'],
        num_replicates=cfg['num_replicates'],
        analysis='RMSF',
        target_selection=cfg['target_selection'],
        ref_selection=ref_sel,
        start_frame=cfg['start_frame'],
        traj_format=cfg['traj_format'],
        top_format=cfg['top_format'],
        group_selections=cfg['rmsf_group_selections'],
        chain_intervals=cfg['chain_intervals'],
        wrap_selection=cfg['wrap_selection'],
        output_dir=os.path.join(cfg['_work_dir'], _ANALYSIS_SUBDIRS['rmsf']),
        parallel_backend=cfg['parallel_backend'],
        n_workers=cfg['n_workers'],
        strict_wrapping=cfg['strict_wrapping'],
    )


def run_2d_rmsd(cfg, dry_run=False):
    """Run 2D-RMSD analysis for all systems."""
    from run_rms_analysis import run_rms_analysis

    print(f"\n{'='*60}")
    print("  2D-RMSD ANALYSIS")
    print(f"{'='*60}")
    print("  NOTE: 2D-RMSD is compute-intensive (O(N²) in frames).")
    print("  For production workloads with multiple systems/replicates,")
    print("  prefer the Nextflow pipeline: nextflow run main.nf")
    print()

    if dry_run:
        print("  [DRY RUN] Would run 2D-RMSD for:", cfg['systems'])
        return

    # 2D-RMSD uses only the first ref_selection (no iteration)
    ref_sel = cfg['ref_selection'][0] if cfg['ref_selection'] else 'protein and backbone'

    run_rms_analysis(
        systems=cfg['systems'],
        variations=cfg['variations'],
        num_replicates=cfg['num_replicates'],
        analysis='2D-RMSD',
        target_selection=cfg['target_selection'],
        ref_selection=ref_sel,
        start_frame=cfg['start_frame'],
        traj_format=cfg['traj_format'],
        top_format=cfg['top_format'],
        wrap_selection=cfg['wrap_selection'],
        output_dir=os.path.join(cfg['_work_dir'], _ANALYSIS_SUBDIRS['2d_rmsd']),
        strict_wrapping=cfg['strict_wrapping'],
    )


def run_rog(cfg, dry_run=False):
    """Run Radius of Gyration analysis for all systems."""
    from run_rog_analysis import run_rog_analysis

    print(f"\n{'='*60}")
    print("  RADIUS OF GYRATION ANALYSIS")
    print(f"{'='*60}")

    if dry_run:
        print("  [DRY RUN] Would run RoG for:", cfg['systems'])
        return

    run_rog_analysis(
        systems=cfg['systems'],
        variations=cfg['variations'],
        num_replicates=cfg['num_replicates'],
        start_frame=cfg['start_frame'],
        traj_format=cfg['traj_format'],
        top_format=cfg['top_format'],
        selection=cfg['rog_selection'],
        time_unit=cfg['time_unit'],
        wrap_selection=cfg['wrap_selection'],
        output_dir=os.path.join(cfg['_work_dir'], _ANALYSIS_SUBDIRS['rog']),
        strict_wrapping=cfg['strict_wrapping'],
    )


def run_hbonds(cfg, dry_run=False):
    """Run Hydrogen Bonds analysis for all systems."""
    from run_hbonds_analysis import run_hbonds_analysis

    print(f"\n{'='*60}")
    print("  HYDROGEN BONDS ANALYSIS")
    print(f"{'='*60}")
    print("  NOTE: H-bond analysis can be very slow for large trajectories.")
    print("  For production workloads with multiple systems/replicates,")
    print("  prefer the Nextflow pipeline: nextflow run main.nf")
    print()

    if dry_run:
        print("  [DRY RUN] Would run H-bonds for:", cfg['systems'])
        return

    run_hbonds_analysis(
        systems=cfg['systems'],
        variations=cfg['variations'],
        num_replicates=cfg['num_replicates'],
        d_a_cutoff=cfg['d_a_cutoff'],
        d_h_a_angle_cutoff=cfg['d_h_a_angle_cutoff'],
        start_frame=cfg['start_frame'],
        traj_format=cfg['traj_format'],
        top_format=cfg['top_format'],
        acceptors_sel=cfg['acceptors_sel'],
        hydrogens_sel=cfg['hydrogens_sel'],
        between_pairs=cfg['between_pairs'],
        update_selections=cfg['update_selections'],
        wrap_selection=cfg['wrap_selection'],
        output_dir=os.path.join(cfg['_work_dir'], _ANALYSIS_SUBDIRS['hbonds']),
        parallel_backend=cfg['parallel_backend'],
        n_workers=cfg['n_workers'],
        strict_wrapping=cfg['strict_wrapping'],
    )


# ─── Plotting runners ────────────────────────────────────────────────────────

def _collect_pickles(work_dir, prefix, cfg, analysis_type=None):
    """
    Find all pickle files matching a prefix for the configured systems.
    Returns a sorted list of absolute paths.

    Parameters
    ----------
    analysis_type : str or None
        Key into ``_ANALYSIS_SUBDIRS`` (e.g. ``'rmsd'``).  When provided,
        pickles are searched inside ``work_dir/<subdir>/``.
    """
    search_dir = os.path.join(work_dir, _ANALYSIS_SUBDIRS[analysis_type]) if analysis_type else work_dir
    patterns = []
    for system in cfg['systems']:
        for variation in cfg['variations'][system]:
            for rep in range(1, cfg['num_replicates'] + 1):
                patterns.append(os.path.join(search_dir, f'{prefix}_{system}_{variation}_rep{rep}*.pkl'))

    pickles = []
    for pat in patterns:
        pickles.extend(glob.glob(pat))
    return sorted(set(pickles))


def plot_rmsd_results(cfg, work_dir, dry_run=False):
    """Plot all RMSD results — individual + comparison group overlays."""
    from plotting.plot_rmsd import plot_rmsd, plot_rmsd_comparison, plot_rmsd_comparison_average

    pickles = _collect_pickles(work_dir, 'rmsd_plot', cfg, analysis_type='rmsd')
    if not pickles:
        print("  No RMSD pickle files found to plot.")
        return

    plot_dir = os.path.join(cfg['plot_dir'], 'rmsd')
    os.makedirs(plot_dir, exist_ok=True)

    # ── Individual plots ──────────────────────────────────────────────────
    for pkl in pickles:
        if dry_run:
            print(f"  [DRY RUN] Would plot: {os.path.basename(pkl)}")
            continue
        print(f"  Plotting: {os.path.basename(pkl)}")
        plot_rmsd(pkl, output_dir=plot_dir, dpi=cfg['dpi'], show_kde=cfg['rmsd_show_kde'])

    # ── Comparison group plots ────────────────────────────────────────────
    if cfg.get('plot_groups'):
        rmsd_dir = os.path.join(work_dir, _ANALYSIS_SUBDIRS['rmsd'])
        _plot_rmsd_comparison_groups(cfg, rmsd_dir, plot_dir, pickles, dry_run)


def plot_rmsf_results(cfg, work_dir, dry_run=False):
    """Plot all RMSF results — individual + comparison group overlays."""
    from plotting.plot_rmsf import plot_rmsf, plot_rmsf_comparison, plot_rmsf_comparison_average

    pickles = _collect_pickles(work_dir, 'rmsf_plot', cfg, analysis_type='rmsf')
    if not pickles:
        print("  No RMSF pickle files found to plot.")
        return

    plot_dir = os.path.join(cfg['plot_dir'], 'rmsf')
    os.makedirs(plot_dir, exist_ok=True)

    # ── Individual plots ──────────────────────────────────────────────────
    for pkl in pickles:
        if dry_run:
            print(f"  [DRY RUN] Would plot: {os.path.basename(pkl)}")
            continue
        print(f"  Plotting: {os.path.basename(pkl)}")
        plot_rmsf(pkl, output_dir=plot_dir, dpi=cfg['dpi'], plot_type=cfg['rmsf_plot_type'])

    # ── Comparison group plots ────────────────────────────────────────────
    if cfg.get('plot_groups'):
        rmsf_dir = os.path.join(work_dir, _ANALYSIS_SUBDIRS['rmsf'])
        _plot_rmsf_comparison_groups(cfg, rmsf_dir, plot_dir, pickles, dry_run)


# ─── Comparison-group helper functions ────────────────────────────────────────

def _find_pickle_for_member(work_dir, prefix, system, variation, rep, sel_idx=None, ref_idx=None):
    """Find a pickle file matching the given system/variation/rep/selection/ref."""
    ref_part = f'_ref{ref_idx}' if ref_idx is not None else ''
    if sel_idx is not None:
        pattern = os.path.join(work_dir, f'{prefix}_{system}_{variation}_rep{rep}{ref_part}_sel{sel_idx}.pkl')
    else:
        pattern = os.path.join(work_dir, f'{prefix}_{system}_{variation}_rep{rep}{ref_part}.pkl')
    matches = glob.glob(pattern)
    if matches:
        return matches[0]
    # Fallback: try wildcard (handles both sel and non-sel cases)
    pattern = os.path.join(work_dir, f'{prefix}_{system}_{variation}_rep{rep}{ref_part}*.pkl')
    matches = glob.glob(pattern)
    return matches[0] if matches else None


def _detect_selection_indices(work_dir, prefix, cfg):
    """Detect all selection indices present in RMSD pickle names.
    Returns a sorted list of sel_idx integers, or [None] if no _selN suffix found."""
    sel_indices = set()
    for system in cfg['systems']:
        for variation in cfg['variations'][system]:
            for rep in range(1, cfg['num_replicates'] + 1):
                pattern = os.path.join(work_dir, f'{prefix}_{system}_{variation}_rep{rep}*_sel*.pkl')
                for match in glob.glob(pattern):
                    basename = os.path.basename(match)
                    # Extract selN from basename
                    import re
                    m = re.search(r'_sel(\d+)\.pkl$', basename)
                    if m:
                        sel_indices.add(int(m.group(1)))
    if sel_indices:
        return sorted(sel_indices)
    return [None]


def _detect_ref_indices(work_dir, prefix, cfg):
    """Detect all ref_selection indices present in RMSD pickle names.
    Returns a sorted list of ref_idx integers, or [None] if no _refN suffix found."""
    ref_indices = set()
    for system in cfg['systems']:
        for variation in cfg['variations'][system]:
            for rep in range(1, cfg['num_replicates'] + 1):
                pattern = os.path.join(work_dir, f'{prefix}_{system}_{variation}_rep{rep}_ref*.pkl')
                for match in glob.glob(pattern):
                    basename = os.path.basename(match)
                    import re
                    m = re.search(r'_ref(\d+)', basename)
                    if m:
                        ref_indices.add(int(m.group(1)))
    if ref_indices:
        return sorted(ref_indices)
    return [None]


def _get_selection_label_from_pickle(pkl_path):
    """Read the 'selection' and 'ref_selection' keys from an RMSD pickle dict."""
    import pickle
    try:
        with open(pkl_path, 'rb') as f:
            data = pickle.load(f)
        if isinstance(data, dict):
            return data.get('selection', ''), data.get('ref_selection', '')
    except Exception:
        pass
    return '', ''


def _plot_rmsd_comparison_groups(cfg, work_dir, plot_dir, all_pickles, dry_run):
    """Generate RMSD comparison plots for each named plot_group."""
    from plotting.plot_rmsd import plot_rmsd_comparison, plot_rmsd_comparison_average

    replicate_mode = cfg.get('replicate_mode', 'separate')
    sel_indices = _detect_selection_indices(work_dir, 'rmsd_plot', cfg)
    ref_indices = _detect_ref_indices(work_dir, 'rmsd_plot', cfg)

    for group_name, members in cfg['plot_groups'].items():
        for ref_idx in ref_indices:
            ref_suffix = f' (ref{ref_idx})' if ref_idx is not None else ''

            for sel_idx in sel_indices:
                sel_suffix = f' (sel{sel_idx})' if sel_idx is not None else ''

                if replicate_mode == 'average':
                    # Average across replicates
                    pickle_groups = []
                    labels = []
                    selection_label = ''
                    ref_selection_label = ''

                    for system, variation in members:
                        rep_pickles = []
                        for rep in range(1, cfg['num_replicates'] + 1):
                            pkl = _find_pickle_for_member(work_dir, 'rmsd_plot', system, variation, rep, sel_idx, ref_idx)
                            if pkl and os.path.exists(pkl):
                                rep_pickles.append(pkl)
                                if not selection_label:
                                    selection_label, ref_selection_label = _get_selection_label_from_pickle(pkl)
                        if rep_pickles:
                            pickle_groups.append(rep_pickles)
                            labels.append(f'{system} / {variation}')

                    if not pickle_groups:
                        print(f"  WARNING: No pickles found for group '{group_name}'{ref_suffix}{sel_suffix}. Skipping.")
                        continue

                    gname = f'{group_name}{ref_suffix}{sel_suffix}'.strip()
                    if dry_run:
                        print(f"  [DRY RUN] Would plot averaged RMSD comparison: {gname}")
                        continue

                    # Build subtitle showing selection and ref_selection
                    subtitle = selection_label
                    if ref_selection_label and ref_selection_label != selection_label:
                        subtitle = f'{selection_label} (aligned on: {ref_selection_label})'

                    print(f"  Plotting averaged RMSD comparison: {gname}")
                    plot_rmsd_comparison_average(
                        pickle_groups, labels, gname,
                        output_dir=plot_dir, dpi=cfg['dpi'],
                        show_kde=cfg['rmsd_show_kde'],
                        selection_label=subtitle,
                    )
                else:
                    # Separate: one comparison plot per replicate
                    for rep in range(1, cfg['num_replicates'] + 1):
                        pkl_files = []
                        labels = []
                        selection_label = ''
                        ref_selection_label = ''

                        for system, variation in members:
                            pkl = _find_pickle_for_member(work_dir, 'rmsd_plot', system, variation, rep, sel_idx, ref_idx)
                            if pkl and os.path.exists(pkl):
                                pkl_files.append(pkl)
                                labels.append(f'{system} / {variation}')
                                if not selection_label:
                                    selection_label, ref_selection_label = _get_selection_label_from_pickle(pkl)

                        if not pkl_files:
                            continue

                        rep_suffix = f'_rep{rep}' if cfg['num_replicates'] > 1 else ''
                        gname = f'{group_name}{ref_suffix}{sel_suffix}{rep_suffix}'.strip()

                        if dry_run:
                            print(f"  [DRY RUN] Would plot RMSD comparison: {gname}")
                            continue

                        # Build subtitle showing selection and ref_selection
                        subtitle = selection_label
                        if ref_selection_label and ref_selection_label != selection_label:
                            subtitle = f'{selection_label} (aligned on: {ref_selection_label})'

                        print(f"  Plotting RMSD comparison: {gname}")
                        plot_rmsd_comparison(
                            pkl_files, labels, gname,
                            output_dir=plot_dir, dpi=cfg['dpi'],
                            show_kde=cfg['rmsd_show_kde'],
                            selection_label=subtitle,
                        )


def _plot_rmsf_comparison_groups(cfg, work_dir, plot_dir, all_pickles, dry_run):
    """Generate RMSF comparison plots for each named plot_group."""
    from plotting.plot_rmsf import plot_rmsf_comparison, plot_rmsf_comparison_average

    replicate_mode = cfg.get('replicate_mode', 'separate')

    # Detect if chain-split pickles exist
    chain_ids = set()
    for pkl in all_pickles:
        basename = os.path.basename(pkl)
        import re
        m = re.search(r'_chain([A-Za-z0-9]+)\.pkl$', basename)
        if m:
            chain_ids.add(m.group(1))

    # If no chains, use [None] to mean "whole protein"
    chain_id_list = sorted(chain_ids) if chain_ids else [None]

    for group_name, members in cfg['plot_groups'].items():
        for chain_id in chain_id_list:
            chain_suffix = f'_chain{chain_id}' if chain_id else ''

            if replicate_mode == 'average':
                pickle_groups = []
                labels = []

                for system, variation in members:
                    rep_pickles = []
                    for rep in range(1, cfg['num_replicates'] + 1):
                        if chain_id:
                            pattern = os.path.join(work_dir, f'rmsf_plot_{system}_{variation}_rep{rep}_chain{chain_id}.pkl')
                        else:
                            pattern = os.path.join(work_dir, f'rmsf_plot_{system}_{variation}_rep{rep}.pkl')
                        matches = glob.glob(pattern)
                        if matches:
                            rep_pickles.append(matches[0])
                    if rep_pickles:
                        pickle_groups.append(rep_pickles)
                        labels.append(f'{system} / {variation}')

                if not pickle_groups:
                    continue

                gname = f'{group_name}{chain_suffix}'
                if dry_run:
                    print(f"  [DRY RUN] Would plot averaged RMSF comparison: {gname}")
                    continue

                print(f"  Plotting averaged RMSF comparison: {gname}")
                plot_rmsf_comparison_average(
                    pickle_groups, labels, gname,
                    output_dir=plot_dir, dpi=cfg['dpi'],
                    selection_label=cfg.get('target_selection'),
                )
            else:
                for rep in range(1, cfg['num_replicates'] + 1):
                    pkl_files = []
                    labels = []

                    for system, variation in members:
                        if chain_id:
                            pattern = os.path.join(work_dir, f'rmsf_plot_{system}_{variation}_rep{rep}_chain{chain_id}.pkl')
                        else:
                            pattern = os.path.join(work_dir, f'rmsf_plot_{system}_{variation}_rep{rep}.pkl')
                        matches = glob.glob(pattern)
                        if matches:
                            pkl_files.append(matches[0])
                            labels.append(f'{system} / {variation}')

                    if not pkl_files:
                        continue

                    rep_suffix = f'_rep{rep}' if cfg['num_replicates'] > 1 else ''
                    gname = f'{group_name}{chain_suffix}{rep_suffix}'

                    if dry_run:
                        print(f"  [DRY RUN] Would plot RMSF comparison: {gname}")
                        continue

                    print(f"  Plotting RMSF comparison: {gname}")
                    plot_rmsf_comparison(
                        pkl_files, labels, gname,
                        output_dir=plot_dir, dpi=cfg['dpi'],
                        selection_label=cfg.get('target_selection'),
                    )


def plot_2d_rmsd_results(cfg, work_dir, dry_run=False):
    """Plot all 2D-RMSD results."""
    from plotting.plot_2d_rmsd import plot_2d_rmsd

    pickles = _collect_pickles(work_dir, '2d_rmsd_plot', cfg, analysis_type='2d_rmsd')
    if not pickles:
        print("  No 2D-RMSD pickle files found to plot.")
        return

    plot_dir = os.path.join(cfg['plot_dir'], '2d_rmsd')
    os.makedirs(plot_dir, exist_ok=True)

    for pkl in pickles:
        if dry_run:
            print(f"  [DRY RUN] Would plot: {os.path.basename(pkl)}")
            continue
        print(f"  Plotting: {os.path.basename(pkl)}")
        plot_2d_rmsd(pkl, output_dir=plot_dir, dpi=cfg['dpi'], cmap=cfg['twod_rmsd_cmap'])


def plot_rog_results(cfg, work_dir, dry_run=False):
    """Plot all RoG results."""
    from plotting.plot_rog import plot_rog

    pickles = _collect_pickles(work_dir, 'rog_plot', cfg, analysis_type='rog')
    if not pickles:
        print("  No RoG pickle files found to plot.")
        return

    plot_dir = os.path.join(cfg['plot_dir'], 'rog')
    os.makedirs(plot_dir, exist_ok=True)

    for pkl in pickles:
        if dry_run:
            print(f"  [DRY RUN] Would plot: {os.path.basename(pkl)}")
            continue
        print(f"  Plotting: {os.path.basename(pkl)}")
        plot_rog(pkl, output_dir=plot_dir, dpi=cfg['dpi'], show_kde=cfg['rog_show_kde'])


def plot_hbonds_results(cfg, work_dir, dry_run=False):
    """Plot all H-bond results."""
    from plotting.plot_hbonds import plot_hbonds

    pickles = _collect_pickles(work_dir, 'hbonds_plot', cfg, analysis_type='hbonds')
    if not pickles:
        print("  No H-bonds pickle files found to plot.")
        return

    plot_dir = os.path.join(cfg['plot_dir'], 'hbonds')
    os.makedirs(plot_dir, exist_ok=True)

    for pkl in pickles:
        if dry_run:
            print(f"  [DRY RUN] Would plot: {os.path.basename(pkl)}")
            continue
        print(f"  Plotting: {os.path.basename(pkl)}")
        plot_hbonds(pkl, output_dir=plot_dir, dpi=cfg['dpi'],
                    time_unit=cfg['time_unit'], top_n=cfg['hbonds_top_n'])


# ─── Pickle validation for plot-only scenarios ──────────────────────────────

# Mapping from analysis name to its pickle file prefix
_ANALYSIS_PICKLE_PREFIX = {
    'rmsd':    'rmsd_plot',
    'rmsf':    'rmsf_plot',
    '2d_rmsd': '2d_rmsd_plot',
    'rog':     'rog_plot',
    'hbonds':  'hbonds_plot',
}

_ANALYSIS_DISPLAY_NAME = {
    'rmsd':    'RMSD',
    'rmsf':    'RMSF',
    '2d_rmsd': '2D-RMSD',
    'rog':     'Radius of Gyration',
    'hbonds':  'Hydrogen Bonds',
}


def validate_plot_pickles(cfg, work_dir):
    """
    For each analysis where plotting is requested but analysis is disabled,
    verify that the required pickle files already exist in work_dir.

    Returns a list of error messages. An empty list means all checks passed.
    """
    errors = []

    for analysis_key in ('rmsd', 'rmsf', '2d_rmsd', 'rog', 'hbonds'):
        plot_flag = cfg.get(f'plot_{analysis_key}', False)
        run_flag = cfg.get(f'run_{analysis_key}', False)

        if not plot_flag or run_flag:
            # Either plotting is off, or analysis will run first — nothing to check
            continue

        # Plotting is ON, analysis is OFF → pickles must already exist
        prefix = _ANALYSIS_PICKLE_PREFIX[analysis_key]
        display_name = _ANALYSIS_DISPLAY_NAME[analysis_key]
        found_pickles = _collect_pickles(work_dir, prefix, cfg, analysis_type=analysis_key)

        if not found_pickles:
            # Build a detailed error message listing every expected file
            expected = []
            for system in cfg['systems']:
                for variation in cfg['variations'][system]:
                    for rep in range(1, cfg['num_replicates'] + 1):
                        expected_pattern = f"{prefix}_{system}_{variation}_rep{rep}*.pkl"
                        expected.append(os.path.join(work_dir, expected_pattern))

            msg_lines = [
                f"  {display_name} plotting is enabled (plot_{analysis_key} = true) "
                f"but {display_name} analysis is disabled (run_{analysis_key} = false).",
                f"  No pre-computed pickle files were found in the work directory.",
                f"  Work directory searched: {work_dir}",
                f"  Expected file pattern(s):",
            ]
            for exp in expected:
                msg_lines.append(f"    - {exp}")
            msg_lines.append(
                f"  Resolution: Either enable run_{analysis_key} = true to compute "
                f"the analysis first, or ensure that the pickle files from a "
                f"previous run are present in the work directory above."
            )
            errors.append("\n".join(msg_lines))

    return errors


# ─── Validation ──────────────────────────────────────────────────────────────

def validate_inputs(cfg):
    """
    Pre-flight checks: verify that expected input files exist before
    starting analyses. Returns (ok_count, missing_list).
    """
    input_dir = os.path.abspath(cfg['input_dir'])
    missing = []
    found = 0

    for system in cfg['systems']:
        if system not in cfg['variations']:
            missing.append(f"System '{system}' listed in systems but has no entry in variations.")
            continue

        for variation in cfg['variations'][system]:
            for rep in range(1, cfg['num_replicates'] + 1):
                top = os.path.join(input_dir, system, variation,
                                   f'{system}_system_{variation}.{cfg["top_format"]}')
                traj, traj_candidates = resolve_trajectory_file(
                    system,
                    variation,
                    rep,
                    cfg['traj_format'],
                    base_dir=input_dir,
                )

                if not os.path.isfile(top):
                    missing.append(f"Topology not found: {top}")
                else:
                    found += 1

                if not os.path.isfile(traj):
                    attempted = ', '.join(
                        os.path.join(input_dir, rel) for rel in traj_candidates
                    )
                    missing.append(f"Trajectory not found (tried: {attempted})")
                else:
                    found += 1

    return found, missing


def print_config_summary(cfg):
    """Print a human-readable summary of the loaded configuration."""
    print(f"\n{'─'*60}")
    print("  CONFIGURATION SUMMARY")
    print(f"{'─'*60}")
    print(f"  Systems:     {cfg['systems']}")
    for s in cfg['systems']:
        print(f"    {s}: {cfg['variations'].get(s, [])}")
    print(f"  Replicates:  {cfg['num_replicates']}")
    print(f"  Format:      {cfg['traj_format']}")
    print(f"  Start frame: {cfg['start_frame']}")
    print(f"  Input dir:   {os.path.abspath(cfg['input_dir'])}")
    print(f"  Output dir:  {os.path.abspath(cfg['output_dir'])}")
    print(f"  Plot dir:    {os.path.abspath(cfg['plot_dir'])}")
    print()
    analyses = []
    if cfg['run_rmsd']:    analyses.append('RMSD')
    if cfg['run_rmsf']:    analyses.append('RMSF')
    if cfg['run_2d_rmsd']: analyses.append('2D-RMSD')
    if cfg['run_rog']:     analyses.append('RoG')
    if cfg['run_hbonds']:  analyses.append('H-bonds')
    print(f"  Analyses:    {', '.join(analyses) if analyses else '(none)'}")
    print(f"  Plotting:    {'yes' if cfg['run_plots'] else 'no (master switch off)'}")
    if cfg['run_plots']:
        plots = []
        if cfg['plot_rmsd']:    plots.append('RMSD')
        if cfg['plot_rmsf']:    plots.append('RMSF')
        if cfg['plot_2d_rmsd']: plots.append('2D-RMSD')
        if cfg['plot_rog']:     plots.append('RoG')
        if cfg['plot_hbonds']:  plots.append('H-bonds')
        print(f"  Plot types:  {', '.join(plots) if plots else '(none)'}")

    if cfg['run_rmsd'] and cfg['time_interval_between_frames']:
        print(f"  RMSD dt:     {cfg['time_interval_between_frames']} ps → {cfg['time_unit']}")
    if cfg['run_rmsd'] and cfg.get('group_selections'):
        print(f"  RMSD sels:   {cfg['group_selections']}")
    if cfg['run_rmsd'] and len(cfg.get('ref_selection', [])) > 1:
        print(f"  RMSD ref:    {cfg['ref_selection']}")
    if cfg['run_rmsf'] and cfg['chain_intervals']:
        print(f"  RMSF chains: {cfg['chain_intervals']}")
    if cfg['run_hbonds']:
        if cfg['between_pairs']:
            print(f"  H-bonds:     pair-focused ({len(cfg['between_pairs'])} pairs)")
        elif cfg['acceptors_sel']:
            print(f"  H-bonds:     atom-focused (A={cfg['acceptors_sel']}, H={cfg['hydrogens_sel']})")
    if cfg.get('plot_groups'):
        print(f"  Plot groups: {list(cfg['plot_groups'].keys())} (mode={cfg['replicate_mode']})")

    # Parallelization summary
    _pb = cfg.get('parallel_backend', 'serial')
    _nw = cfg.get('n_workers')
    if _pb != 'serial':
        print(f"  Parallel:    backend='{_pb}', n_workers={_nw if _nw else 'auto'}")
        print(f"               Applies to: RMSD, H-bonds")
        print(f"               Serial only: RMSF, 2D-RMSD, RoG")
    else:
        print(f"  Parallel:    serial (no parallelization)")

    print(f"{'─'*60}\n")


# ─── Main execution ──────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Run the full MD analysis pipeline from an INI configuration file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Example:
    python executor.py pipeline.ini
    python executor.py pipeline.ini --dry-run
    python executor.py pipeline.ini --force

The INI file should contain sections: [systems], [analysis], and optionally
[paths], [selections], [rmsd], [rmsf], [hbonds], [plotting].
See example_pipeline.ini for a complete template."""
    )

    parser.add_argument('config', type=str, help='Path to the .INI configuration file')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print what would be done without executing analyses')
    parser.add_argument('--force', action='store_true',
                        help='Delete existing pickle files before running (re-compute everything)')
    parser.add_argument('--analysis-only', action='store_true',
                        help='Run analyses but skip plotting')
    parser.add_argument('--plot-only', action='store_true',
                        help='Skip analyses and only generate plots from existing pickles')
    parser.add_argument('--strict-wrapping', action='store_true',
                        help='Fail if trajectory lacks fragment data for PBC wrapping (default: warn and use residue-based fallback)')

    args = parser.parse_args()

    # Load config
    cfg = load_config(args.config)
    cfg['strict_wrapping'] = args.strict_wrapping
    if args.analysis_only:
        cfg['run_plots'] = False
    if args.plot_only:
        cfg['run_plots'] = True
        # In plot-only mode, enable all per-analysis plot flags so that
        # every analysis type that has pickles will be plotted.
        for key in ('plot_rmsd', 'plot_rmsf', 'plot_2d_rmsd', 'plot_rog', 'plot_hbonds'):
            cfg[key] = True

    print_config_summary(cfg)

    # Validate inputs (unless plot-only)
    if not args.plot_only:
        found, missing = validate_inputs(cfg)
        if missing:
            print("WARNING: Some expected input files are missing:")
            for m in missing:
                print(f"  - {m}")
            print()
            if found == 0:
                print("ERROR: No input files found at all. Check your [paths] input_dir and file naming.")
                sys.exit(1)
        else:
            print(f"Input validation passed: {found} files found.\n")

    # Set up work directory with symlinks
    work_dir = setup_workdir(cfg)
    cfg['_work_dir'] = work_dir          # used by analysis runners for output_dir
    original_dir = os.getcwd()
    os.chdir(work_dir)
    print(f"Working directory: {work_dir}\n")

    # Validate that pickle files exist for plot-without-analysis scenarios
    if cfg['run_plots'] and not args.dry_run:
        pickle_errors = validate_plot_pickles(cfg, work_dir)
        if pickle_errors:
            print(f"\n{'!'*60}")
            print("  MISSING PICKLE FILES — Cannot generate requested plots")
            print(f"{'!'*60}\n")
            for err in pickle_errors:
                print(err)
                print()
            print(
                "The pipeline cannot proceed because plotting was requested "
                "for analyses that were not run, and no pre-existing results "
                "were found.\n"
                "Exiting."
            )
            os.chdir(original_dir)
            sys.exit(1)

    start_time = _time.time()

    try:
        # ── Force mode: clean up existing pickles ─────────────────────────
        if args.force and not args.dry_run and not args.plot_only:
            print("Force mode: removing existing pickle files in work directory...")
            for subdir in _ANALYSIS_SUBDIRS.values():
                subdir_path = os.path.join(work_dir, subdir)
                if not os.path.isdir(subdir_path):
                    continue
                for pkl in glob.glob(os.path.join(subdir_path, '*.pkl')):
                    os.remove(pkl)
                    print(f"  Removed: {subdir}/{os.path.basename(pkl)}")
            # Also remove aligned trajectory files from the work root
            for fit in glob.glob(os.path.join(work_dir, 'rmsfit_*.dcd')) + \
                        glob.glob(os.path.join(work_dir, 'rmsfit_*.xtc')):
                os.remove(fit)
                print(f"  Removed: {os.path.basename(fit)}")
            print()

        # ── Run analyses ──────────────────────────────────────────────────
        if not args.plot_only:
            if cfg['run_rmsd']:
                run_rmsd(cfg, dry_run=args.dry_run)

            if cfg['run_rmsf']:
                run_rmsf(cfg, dry_run=args.dry_run)

            if cfg['run_2d_rmsd']:
                run_2d_rmsd(cfg, dry_run=args.dry_run)

            if cfg['run_rog']:
                run_rog(cfg, dry_run=args.dry_run)

            if cfg['run_hbonds']:
                run_hbonds(cfg, dry_run=args.dry_run)

            # Clean up ephemeral aligned-trajectory files (rmsfit_*)
            # so the output tree matches Nextflow's isolated processes.
            if not args.dry_run:
                _clean_ephemeral_files(work_dir)

        # ── Run plots ─────────────────────────────────────────────────────
        if cfg['run_plots']:
            print(f"\n{'='*60}")
            print("  GENERATING PLOTS")
            print(f"{'='*60}")

            if cfg['plot_rmsd']:
                plot_rmsd_results(cfg, work_dir, dry_run=args.dry_run)

            if cfg['plot_rmsf']:
                plot_rmsf_results(cfg, work_dir, dry_run=args.dry_run)

            if cfg['plot_2d_rmsd']:
                plot_2d_rmsd_results(cfg, work_dir, dry_run=args.dry_run)

            if cfg['plot_rog']:
                plot_rog_results(cfg, work_dir, dry_run=args.dry_run)

            if cfg['plot_hbonds']:
                plot_hbonds_results(cfg, work_dir, dry_run=args.dry_run)

            # Clean up intermediate files, preserving analysis outputs
            if not args.dry_run:
                cleanup_work_directory(work_dir, verbose=True)

    finally:
        os.chdir(original_dir)

    elapsed = _time.time() - start_time
    minutes, seconds = divmod(elapsed, 60)

    print(f"\n{'='*60}")
    print(f"  PIPELINE COMPLETE — {int(minutes)}m {seconds:.1f}s elapsed")
    print(f"{'='*60}")
    print(f"  Pickle files: {work_dir}/")
    print(f"  Plots:        {os.path.abspath(cfg['plot_dir'])}/")
    print()

    return 0


if __name__ == '__main__':
    sys.exit(main())

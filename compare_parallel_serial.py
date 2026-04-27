#!/usr/bin/env python3
"""Compare pickle outputs from serial and parallel MDAnalysis runs.

This script validates that parallelized analysis results (using the
``multiprocessing`` or ``dask`` backends in MDAnalysis 2.8+) are
numerically equivalent to serial results.

**This is NOT the same as ``compare_outputs.py``**, which compares
Python (executor.py) vs Nextflow-orchestrated runs.  This script
compares two runs of the *same* orchestration method that differ only
in their ``[parallelization]`` settings.

Usage — automatic mode
----------------------
    # Run both serial and parallel pipelines, then compare:
    python compare_parallel_serial.py \\
        --config test_comparison.ini \\
        --backend multiprocessing \\
        --n-workers 4

    # Speed optimisation — reuse existing serial results:
    python compare_parallel_serial.py \\
        --config test_comparison.ini \\
        --backend multiprocessing \\
        --n-workers 4 \\
        --dir-serial results_comparison/work

Usage — manual mode
-------------------
    # Pre-existing result directories (e.g. from two executor.py runs):
    python compare_parallel_serial.py \\
        --dir-serial   results_serial/work \\
        --dir-parallel results_parallel/work

Automatic mode creates temporary output directories
(``_cps_serial/`` and ``_cps_parallel/``) which are kept after the run
so you can inspect the results.  Pass ``--clean`` to remove them.

Exit codes
----------
  0  All matched pickles are numerically equivalent.
  1  At least one pickle pair has unacceptable differences.
  2  Structural mismatch (files missing from one tree).
"""
from __future__ import annotations

import argparse
import configparser
import importlib
import os
import pickle
import shutil
import subprocess
import sys
import tempfile
import time as _time
from dataclasses import dataclass, field
from typing import Any

import numpy as np


# ─── Robust unpickler ────────────────────────────────────────────────────────
# See compare_outputs.py for a full explanation.  Pickle files created by
# ``python script.py`` store custom classes under ``__main__``, which breaks
# when a different script tries to load them.

_ANALYSIS_MODULES = (
    'run_rog_analysis',
    'run_rms_analysis',
    'run_hbonds_analysis',
)


class _AnalysisUnpickler(pickle.Unpickler):
    """Unpickler that redirects ``__main__`` class lookups to analysis modules."""

    def find_class(self, module: str, name: str) -> type:
        if module == '__main__':
            for mod_name in _ANALYSIS_MODULES:
                try:
                    mod = importlib.import_module(mod_name)
                    return getattr(mod, name)
                except (ImportError, AttributeError):
                    continue
        return super().find_class(module, name)


def _safe_pickle_load(filepath: str) -> Any:
    """Load a pickle file, handling ``__main__`` module references gracefully."""
    with open(filepath, 'rb') as f:
        try:
            return pickle.load(f)
        except (AttributeError, ModuleNotFoundError):
            f.seek(0)
            return _AnalysisUnpickler(f).load()


# ─── Helpers ──────────────────────────────────────────────────────────────────

_SKIP_PATTERNS = {
    'pipeline_info',
    'plots',
    'rmsfit_',
    '.nextflow',
    '__pycache__',
}


def _should_skip(path: str) -> bool:
    """Return True if *path* should be excluded from comparison."""
    base = os.path.basename(path)
    if not base.endswith('.pkl'):
        return True
    return any(pat in path for pat in _SKIP_PATTERNS)


def _collect_pickle_tree(root: str) -> dict[str, str]:
    """Return ``{relative_path: absolute_path}`` for all ``.pkl`` files under *root*."""
    tree: dict[str, str] = {}
    root = os.path.abspath(root)
    for dirpath, _dirs, filenames in os.walk(root):
        for fname in filenames:
            full = os.path.join(dirpath, fname)
            rel = os.path.relpath(full, root)
            if _should_skip(rel):
                continue
            tree[rel] = full
    return tree


# ─── Comparison engine ────────────────────────────────────────────────────────

@dataclass
class ComparisonResult:
    """Result of comparing a single pickle pair."""
    filename: str
    status: str = 'MATCH'       # MATCH | MISMATCH | ONLY_SERIAL | ONLY_PARALLEL
    details: str = ''
    max_abs_diff: float | None = None
    max_rel_diff: float | None = None

    @property
    def ok(self) -> bool:
        return self.status == 'MATCH'


# Maximum recursion depth for _compare_values to prevent stack overflow
# on deeply nested or circular MDAnalysis objects.
_MAX_COMPARE_DEPTH = 30


def _compare_values(a: Any, b: Any, path: str = '', rtol: float = 1e-5,
                    atol: float = 1e-8, _depth: int = 0,
                    _seen: set | None = None) -> list[str]:
    """Recursively compare two Python objects, returning a list of differences.

    Uses *_depth* and *_seen* (set of ``id()`` pairs) to guard against
    infinite recursion caused by circular references in MDAnalysis
    ``Topology`` objects and similar structures.
    """
    diffs: list[str] = []

    # --- depth guard ---
    if _depth > _MAX_COMPARE_DEPTH:
        return diffs  # treat as equal at excessive depth

    # --- cycle guard (only for mutable, identity-bearing objects) ---
    if _seen is None:
        _seen = set()

    if type(a) is not type(b):
        diffs.append(f"{path}: type mismatch ({type(a).__name__} vs {type(b).__name__})")
        return diffs

    if isinstance(a, np.ndarray):
        if a.shape != b.shape:
            diffs.append(f"{path}: shape mismatch {a.shape} vs {b.shape}")
        elif a.dtype.kind in ('U', 'S', 'O'):
            # String or object arrays — element-wise equality
            if not np.array_equal(a, b):
                mismatches = np.sum(a != b)
                diffs.append(f"{path}: {mismatches} element(s) differ in string/object array")
        elif not np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True):
            abs_diff = np.abs(a - b)
            max_abs = float(np.nanmax(abs_diff))
            with np.errstate(divide='ignore', invalid='ignore'):
                rel_diff = np.where(np.abs(b) > 0, abs_diff / np.abs(b), 0.0)
            max_rel = float(np.nanmax(rel_diff))
            diffs.append(
                f"{path}: arrays differ (max |Δ| = {max_abs:.2e}, "
                f"max relative = {max_rel:.2e})"
            )
        return diffs

    if isinstance(a, dict):
        pair_id = (id(a), id(b))
        if pair_id in _seen:
            return diffs  # circular reference — skip
        _seen.add(pair_id)
        all_keys = set(a) | set(b)
        for k in sorted(all_keys, key=str):
            if k not in a:
                diffs.append(f"{path}['{k}']: only in parallel")
            elif k not in b:
                diffs.append(f"{path}['{k}']: only in serial")
            else:
                diffs.extend(_compare_values(a[k], b[k], f"{path}['{k}']", rtol, atol,
                                             _depth + 1, _seen))
        return diffs

    if isinstance(a, (list, tuple)):
        if len(a) != len(b):
            diffs.append(f"{path}: length mismatch ({len(a)} vs {len(b)})")
        else:
            for i, (ai, bi) in enumerate(zip(a, b)):
                diffs.extend(_compare_values(ai, bi, f"{path}[{i}]", rtol, atol,
                                             _depth + 1, _seen))
        return diffs

    # Handle MDAnalysis result objects that store data in .results
    if hasattr(a, 'results') and hasattr(b, 'results'):
        pair_id = (id(a), id(b))
        if pair_id in _seen:
            return diffs
        _seen.add(pair_id)
        diffs.extend(_compare_values(vars(a.results), vars(b.results),
                                     f"{path}.results", rtol, atol,
                                     _depth + 1, _seen))
        return diffs

    # For objects with __dict__, compare attribute-by-attribute
    if hasattr(a, '__dict__') and hasattr(b, '__dict__'):
        pair_id = (id(a), id(b))
        if pair_id in _seen:
            return diffs  # circular reference — skip
        _seen.add(pair_id)
        return _compare_values(a.__dict__, b.__dict__, path or type(a).__name__,
                               rtol, atol, _depth + 1, _seen)

    # Scalar comparison
    if isinstance(a, (int, float, np.integer, np.floating)):
        if not np.isclose(a, b, rtol=rtol, atol=atol, equal_nan=True):
            diffs.append(f"{path}: {a} vs {b} (Δ = {abs(a - b):.2e})")
        return diffs

    if a != b:
        diffs.append(f"{path}: {a!r} != {b!r}")
    return diffs


def compare_pickle_pair(path_serial: str, path_parallel: str,
                        filename: str, rtol: float, atol: float) -> ComparisonResult:
    """Load and compare a matched pair of pickle files."""
    try:
        data_serial = _safe_pickle_load(path_serial)
        data_parallel = _safe_pickle_load(path_parallel)
    except Exception as e:
        return ComparisonResult(filename, 'MISMATCH', f'Error loading pickle: {e}')

    diffs = _compare_values(data_serial, data_parallel, filename, rtol, atol)

    if diffs:
        # Extract numeric diff stats
        max_abs = None
        max_rel = None
        for d in diffs:
            if 'max |Δ|' in d:
                import re
                m = re.search(r'max \|Δ\| = ([\d.eE+-]+)', d)
                if m:
                    val = float(m.group(1))
                    max_abs = max(max_abs or 0, val)
                m = re.search(r'max relative = ([\d.eE+-]+)', d)
                if m:
                    val = float(m.group(1))
                    max_rel = max(max_rel or 0, val)

        return ComparisonResult(
            filename, 'MISMATCH',
            '\n'.join(f'  - {d}' for d in diffs),
            max_abs_diff=max_abs,
            max_rel_diff=max_rel,
        )

    return ComparisonResult(filename, 'MATCH')


def compare_trees(dir_serial: str, dir_parallel: str,
                  rtol: float = 1e-5, atol: float = 1e-8) -> list[ComparisonResult]:
    """Compare all pickle files between serial and parallel result trees."""
    serial_tree = _collect_pickle_tree(dir_serial)
    parallel_tree = _collect_pickle_tree(dir_parallel)

    all_keys = sorted(set(serial_tree) | set(parallel_tree))
    results: list[ComparisonResult] = []

    for key in all_keys:
        if key in serial_tree and key not in parallel_tree:
            results.append(ComparisonResult(key, 'ONLY_SERIAL'))
        elif key not in serial_tree and key in parallel_tree:
            results.append(ComparisonResult(key, 'ONLY_PARALLEL'))
        else:
            result = compare_pickle_pair(
                serial_tree[key], parallel_tree[key], key, rtol, atol
            )
            results.append(result)

    return results


# ─── Report ──────────────────────────────────────────────────────────────────

def format_report(results: list[ComparisonResult],
                  dir_serial: str, dir_parallel: str) -> str:
    """Format a Markdown comparison report."""
    lines: list[str] = []
    lines.append('# Parallel vs Serial Comparison Report')
    lines.append('')
    lines.append(f'**Serial directory:**  `{os.path.abspath(dir_serial)}`')
    lines.append(f'**Parallel directory:** `{os.path.abspath(dir_parallel)}`')
    lines.append('')

    matched = [r for r in results if r.status == 'MATCH']
    mismatched = [r for r in results if r.status == 'MISMATCH']
    only_serial = [r for r in results if r.status == 'ONLY_SERIAL']
    only_parallel = [r for r in results if r.status == 'ONLY_PARALLEL']

    lines.append(f'## Summary')
    lines.append('')
    lines.append(f'| Status | Count |')
    lines.append(f'|--------|-------|')
    lines.append(f'| ✅ Match | {len(matched)} |')
    lines.append(f'| ❌ Mismatch | {len(mismatched)} |')
    lines.append(f'| ⚠️  Only in serial | {len(only_serial)} |')
    lines.append(f'| ⚠️  Only in parallel | {len(only_parallel)} |')
    lines.append(f'| **Total** | **{len(results)}** |')
    lines.append('')

    if mismatched:
        lines.append('## Mismatches')
        lines.append('')
        for r in mismatched:
            lines.append(f'### `{r.filename}`')
            if r.max_abs_diff is not None:
                lines.append(f'- Max absolute difference: {r.max_abs_diff:.2e}')
            if r.max_rel_diff is not None:
                lines.append(f'- Max relative difference: {r.max_rel_diff:.2e}')
            lines.append('')
            lines.append(f'```')
            lines.append(r.details)
            lines.append(f'```')
            lines.append('')

    if only_serial:
        lines.append('## Files only in serial')
        lines.append('')
        for r in only_serial:
            lines.append(f'- `{r.filename}`')
        lines.append('')

    if only_parallel:
        lines.append('## Files only in parallel')
        lines.append('')
        for r in only_parallel:
            lines.append(f'- `{r.filename}`')
        lines.append('')

    if not mismatched and not only_serial and not only_parallel:
        lines.append('## Result')
        lines.append('')
        lines.append('**All pickle files are numerically equivalent.** ✅')
        lines.append('')
        lines.append('Parallel execution with MDAnalysis 2.8+ native backends ')
        lines.append('produces results identical to serial execution within the ')
        lines.append('specified tolerances.')

    return '\n'.join(lines)


# ─── Auto-mode helpers ────────────────────────────────────────────────────────

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def _create_temp_ini(
    original_path: str,
    output_dir: str,
    backend: str = 'serial',
    n_workers: int | None = None,
) -> str:
    """Return path to a temporary INI derived from *original_path*.

    The copy overrides ``[paths] output_dir``, ``[parallelization]`` settings,
    and forces ``run_plots = false`` so the comparison runs as fast as possible.
    """
    cp = configparser.RawConfigParser(inline_comment_prefixes=())
    cp.read(original_path)

    # Override output directory
    if not cp.has_section('paths'):
        cp.add_section('paths')
    cp.set('paths', 'output_dir', os.path.abspath(output_dir))

    # Override parallelization
    if not cp.has_section('parallelization'):
        cp.add_section('parallelization')
    cp.set('parallelization', 'backend', backend)
    cp.set('parallelization', 'n_workers',
           str(n_workers) if n_workers else 'none')

    # Disable plotting for speed
    if cp.has_section('analysis'):
        cp.set('analysis', 'run_plots', 'false')

    fd, tmp_path = tempfile.mkstemp(suffix='.ini', prefix=f'cps_{backend}_')
    with os.fdopen(fd, 'w') as f:
        cp.write(f)
    return tmp_path


def _run_executor(ini_path: str) -> int:
    """Run ``executor.py --force`` with the given INI and return its exit code.

    stdout / stderr are streamed straight to the terminal so the user can
    follow progress.
    """
    executor = os.path.join(_SCRIPT_DIR, 'executor.py')
    cmd = [sys.executable, executor, ini_path, '--force']
    return subprocess.call(cmd)


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Compare serial and parallel MDAnalysis analysis outputs.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Automatic mode — run both pipelines and compare:
    python compare_parallel_serial.py \\
        --config test_comparison.ini \\
        --backend multiprocessing --n-workers 4

Reuse existing serial results (fastest):
    python compare_parallel_serial.py \\
        --config test_comparison.ini \\
        --backend multiprocessing --n-workers 4 \\
        --dir-serial results_comparison/work

Manual mode — compare pre-existing directories:
    python compare_parallel_serial.py \\
        --dir-serial  results_serial/work \\
        --dir-parallel results_parallel/work

This script is NOT for comparing Python vs Nextflow runs.
Use compare_outputs.py for that purpose."""
    )

    # ── Auto-mode arguments ───────────────────────────────────────────────
    auto = parser.add_argument_group('automatic mode')
    auto.add_argument(
        '--config', metavar='INI',
        help='Path to an executor.py INI config file.  When given, the '
             'script runs the pipeline twice (serial + parallel) and '
             'compares the outputs automatically.')
    auto.add_argument(
        '--backend', default='multiprocessing',
        choices=('multiprocessing', 'dask'),
        help='Parallel backend for the non-serial run '
             '(default: multiprocessing)')
    auto.add_argument(
        '--n-workers', type=int, default=None,
        help='Number of workers for the parallel run '
             '(default: CPU count / 2)')
    auto.add_argument(
        '--clean', action='store_true',
        help='Remove temporary result directories after comparison')

    # ── Manual-mode arguments ─────────────────────────────────────────────
    manual = parser.add_argument_group('manual mode')
    manual.add_argument(
        '--dir-serial', metavar='DIR',
        help='Directory containing serial run pickle outputs')
    manual.add_argument(
        '--dir-parallel', metavar='DIR',
        help='Directory containing parallel run pickle outputs')

    # ── Common arguments ──────────────────────────────────────────────────
    parser.add_argument(
        '--rtol', type=float, default=1e-5,
        help='Relative tolerance for numeric comparison (default: 1e-5)')
    parser.add_argument(
        '--atol', type=float, default=1e-8,
        help='Absolute tolerance for numeric comparison (default: 1e-8)')
    parser.add_argument(
        '-o', '--output', default=None,
        help='Output report file (default: stdout)')

    args = parser.parse_args()

    # ── Resolve mode and validate ─────────────────────────────────────────
    auto_mode = args.config is not None
    manual_mode = (args.dir_serial is not None and
                   args.dir_parallel is not None)

    if not auto_mode and not manual_mode:
        parser.error(
            'Either --config (automatic mode) or both --dir-serial and '
            '--dir-parallel (manual mode) are required.')

    if auto_mode and args.dir_parallel:
        parser.error(
            '--dir-parallel cannot be used together with --config.  '
            'In automatic mode, the parallel results are produced by '
            'the script itself.')

    cleanup_dirs: list[str] = []     # dirs to remove if --clean

    if auto_mode:
        if not os.path.isfile(args.config):
            print(f'ERROR: Config file not found: {args.config}',
                  file=sys.stderr)
            return 2

        t0 = _time.time()

        # ── Serial results ────────────────────────────────────────────
        if args.dir_serial:
            # Reuse existing serial results
            dir_serial = args.dir_serial
            if not os.path.isdir(dir_serial):
                print(f'ERROR: Serial directory not found: {dir_serial}',
                      file=sys.stderr)
                return 2
            print(f'\nReusing existing serial results: {dir_serial}\n')
        else:
            serial_outdir = os.path.abspath('_cps_serial')
            cleanup_dirs.append(serial_outdir)
            print(f'\n{"="*60}')
            print('  SERIAL RUN')
            print(f'{"="*60}\n')
            tmp_ini = _create_temp_ini(args.config, serial_outdir,
                                       'serial', None)
            try:
                rc = _run_executor(tmp_ini)
            finally:
                os.unlink(tmp_ini)
            if rc != 0:
                print('ERROR: Serial pipeline failed.', file=sys.stderr)
                return 2
            dir_serial = os.path.join(serial_outdir, 'work')

        # ── Parallel results ──────────────────────────────────────────
        parallel_outdir = os.path.abspath('_cps_parallel')
        cleanup_dirs.append(parallel_outdir)
        print(f'\n{"="*60}')
        print(f'  PARALLEL RUN  (backend={args.backend}, '
              f'n_workers={args.n_workers or "auto"})')
        print(f'{"="*60}\n')
        tmp_ini = _create_temp_ini(args.config, parallel_outdir,
                                   args.backend, args.n_workers)
        try:
            rc = _run_executor(tmp_ini)
        finally:
            os.unlink(tmp_ini)
        if rc != 0:
            print('ERROR: Parallel pipeline failed.', file=sys.stderr)
            return 2

        dir_parallel = os.path.join(parallel_outdir, 'work')

        elapsed = _time.time() - t0
        m, s = divmod(elapsed, 60)
        print(f'\nBoth runs completed in {int(m)}m {s:.1f}s.\n')
    else:
        # Manual mode
        dir_serial = args.dir_serial
        dir_parallel = args.dir_parallel

    # ── Validate directories ──────────────────────────────────────────────
    if not os.path.isdir(dir_serial):
        print(f'ERROR: Serial directory not found: {dir_serial}',
              file=sys.stderr)
        return 2
    if not os.path.isdir(dir_parallel):
        print(f'ERROR: Parallel directory not found: {dir_parallel}',
              file=sys.stderr)
        return 2

    # ── Run comparison ────────────────────────────────────────────────────
    results = compare_trees(dir_serial, dir_parallel,
                            rtol=args.rtol, atol=args.atol)

    if not results:
        print('No pickle files found in either directory.', file=sys.stderr)
        return 2

    # ── Generate report ───────────────────────────────────────────────────
    report = format_report(results, dir_serial, dir_parallel)

    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(report)
        print(f'Report written to {args.output}', file=sys.stderr)
    else:
        print(report)

    # ── Cleanup ───────────────────────────────────────────────────────────
    if args.clean:
        for d in cleanup_dirs:
            if os.path.isdir(d):
                shutil.rmtree(d)
                print(f'Cleaned: {d}', file=sys.stderr)

    # ── Exit code ─────────────────────────────────────────────────────────
    has_mismatches = any(r.status == 'MISMATCH' for r in results)
    has_structural = any(r.status in ('ONLY_SERIAL', 'ONLY_PARALLEL')
                         for r in results)

    if has_structural:
        return 2
    if has_mismatches:
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3
"""Compare pickle outputs from executor.py and Nextflow runs.

Usage
-----
    python compare_outputs.py \\
        --dir-python  results/work \\
        --dir-nextflow results_nf

Walks both directory trees, pairs pickles by filename, and does a semantic
comparison using ``numpy.allclose`` for numeric arrays and ``==`` for other
Python objects.  Produces a Markdown report (printed to stdout or written
to a file).

Exit codes
----------
  0  All matched pickles are equivalent.
  1  At least one pickle pair differs.
  2  Structural mismatch (files present in one tree but not the other).
"""
from __future__ import annotations

import argparse
import importlib
import os
import pickle
import sys
from dataclasses import dataclass
from typing import Any

import numpy as np


# ─── Robust unpickler ────────────────────────────────────────────────────────
# Pickle files created by scripts run as ``python script.py`` (e.g. via
# Nextflow) store custom classes under the ``__main__`` module.  When a
# *different* script loads those pickles, ``__main__`` points to the wrong
# module and ``pickle.load`` fails.  This unpickler intercepts such lookups
# and resolves them from the known analysis modules instead.

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

# Files / directories to skip during comparison.
_SKIP_PATTERNS = {
    'pipeline_info',
    'plots',
    'rmsfit_',
    '.nextflow',
    '__pycache__',
}


def _should_skip(path: str) -> bool:
    """Return True if *path* (basename) should be excluded from comparison."""
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
    status: str = 'MATCH'         # MATCH | MISMATCH | ONLY_PYTHON | ONLY_NEXTFLOW
    details: str = ''

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

    # --- cycle guard ---
    if _seen is None:
        _seen = set()

    if type(a) is not type(b):
        diffs.append(f"{path}: type mismatch ({type(a).__name__} vs {type(b).__name__})")
        return diffs

    if isinstance(a, np.ndarray):
        if a.shape != b.shape:
            diffs.append(f"{path}: shape mismatch {a.shape} vs {b.shape}")
        elif not np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True):
            max_diff = float(np.nanmax(np.abs(a - b)))
            diffs.append(f"{path}: arrays differ (max abs diff = {max_diff:.2e})")
        return diffs

    if isinstance(a, dict):
        pair_id = (id(a), id(b))
        if pair_id in _seen:
            return diffs  # circular reference — skip
        _seen.add(pair_id)
        all_keys = set(a) | set(b)
        for k in sorted(all_keys, key=str):
            subpath = f"{path}.{k}" if path else str(k)
            if k not in a:
                diffs.append(f"{subpath}: missing in left (Python)")
            elif k not in b:
                diffs.append(f"{subpath}: missing in right (Nextflow)")
            else:
                diffs.extend(_compare_values(a[k], b[k], subpath, rtol, atol,
                                             _depth + 1, _seen))
        return diffs

    if isinstance(a, (list, tuple)):
        if len(a) != len(b):
            diffs.append(f"{path}: length mismatch ({len(a)} vs {len(b)})")
            return diffs
        for i, (va, vb) in enumerate(zip(a, b)):
            diffs.extend(_compare_values(va, vb, f"{path}[{i}]", rtol, atol,
                                         _depth + 1, _seen))
        return diffs

    if isinstance(a, float):
        if not np.isclose(a, b, rtol=rtol, atol=atol, equal_nan=True):
            diffs.append(f"{path}: {a} != {b}")
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

    # Generic attribute-bearing objects — compare __dict__
    if hasattr(a, '__dict__') and hasattr(b, '__dict__'):
        pair_id = (id(a), id(b))
        if pair_id in _seen:
            return diffs  # circular reference — skip
        _seen.add(pair_id)
        diffs.extend(_compare_values(a.__dict__, b.__dict__, path, rtol, atol,
                                     _depth + 1, _seen))
        return diffs

    # Scalar fallback
    if a != b:
        diffs.append(f"{path}: {a!r} != {b!r}")
    return diffs


def compare_pickle_pair(path_a: str, path_b: str, filename: str,
                        rtol: float = 1e-5, atol: float = 1e-8) -> ComparisonResult:
    """Load two pickle files and compare their contents semantically."""
    try:
        obj_a = _safe_pickle_load(path_a)
    except Exception as exc:
        return ComparisonResult(filename, 'MISMATCH', f"Failed to load Python pickle: {exc}")

    try:
        obj_b = _safe_pickle_load(path_b)
    except Exception as exc:
        return ComparisonResult(filename, 'MISMATCH', f"Failed to load Nextflow pickle: {exc}")

    diffs = _compare_values(obj_a, obj_b, rtol=rtol, atol=atol)
    if diffs:
        detail_str = '; '.join(diffs[:5])
        if len(diffs) > 5:
            detail_str += f' ... and {len(diffs) - 5} more'
        return ComparisonResult(filename, 'MISMATCH', detail_str)
    return ComparisonResult(filename, 'MATCH', 'All values within tolerance')


# ─── Report ──────────────────────────────────────────────────────────────────

def _render_markdown(results: list[ComparisonResult], dir_py: str, dir_nf: str) -> str:
    """Render a Markdown comparison report."""
    lines: list[str] = [
        '# Output Comparison Report',
        '',
        f'**Python dir**: `{dir_py}`  ',
        f'**Nextflow dir**: `{dir_nf}`  ',
        '',
    ]

    n_match = sum(1 for r in results if r.status == 'MATCH')
    n_mismatch = sum(1 for r in results if r.status == 'MISMATCH')
    n_only_py = sum(1 for r in results if r.status == 'ONLY_PYTHON')
    n_only_nf = sum(1 for r in results if r.status == 'ONLY_NEXTFLOW')

    lines.append(f'**Summary**: {n_match} match, {n_mismatch} mismatch, '
                 f'{n_only_py} only-python, {n_only_nf} only-nextflow')
    lines.append('')

    if not results:
        lines.append('_No pickle files found to compare._')
        return '\n'.join(lines)

    lines.append('| File | Status | Details |')
    lines.append('|------|--------|---------|')
    for r in sorted(results, key=lambda x: x.filename):
        icon = {'MATCH': ':white_check_mark:', 'MISMATCH': ':x:',
                'ONLY_PYTHON': ':warning:', 'ONLY_NEXTFLOW': ':warning:'}.get(r.status, '')
        # Escape pipe characters in details
        safe_details = r.details.replace('|', '\\|')
        lines.append(f'| {r.filename} | {icon} {r.status} | {safe_details} |')

    lines.append('')
    return '\n'.join(lines)


# ─── Main ────────────────────────────────────────────────────────────────────

def run_comparison(dir_python: str, dir_nextflow: str, *,
                   rtol: float = 1e-5, atol: float = 1e-8) -> list[ComparisonResult]:
    """Compare all pickles between two output directories."""
    tree_py = _collect_pickle_tree(dir_python)
    tree_nf = _collect_pickle_tree(dir_nextflow)

    all_keys = sorted(set(tree_py) | set(tree_nf))
    results: list[ComparisonResult] = []

    for key in all_keys:
        if key in tree_py and key in tree_nf:
            results.append(compare_pickle_pair(tree_py[key], tree_nf[key], key,
                                                rtol=rtol, atol=atol))
        elif key in tree_py:
            results.append(ComparisonResult(key, 'ONLY_PYTHON', 'Not found in Nextflow output'))
        else:
            results.append(ComparisonResult(key, 'ONLY_NEXTFLOW', 'Not found in Python output'))

    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Compare pickle outputs from executor.py and Nextflow runs.')
    parser.add_argument('--dir-python', required=True,
                        help='Root output directory from executor.py (e.g. results/work)')
    parser.add_argument('--dir-nextflow', required=True,
                        help='Root output directory from Nextflow (e.g. results_nf)')
    parser.add_argument('--rtol', type=float, default=1e-5,
                        help='Relative tolerance for numeric comparison (default: 1e-5)')
    parser.add_argument('--atol', type=float, default=1e-8,
                        help='Absolute tolerance for numeric comparison (default: 1e-8)')
    parser.add_argument('-o', '--output', default=None,
                        help='Write Markdown report to file (default: stdout)')
    args = parser.parse_args()

    results = run_comparison(args.dir_python, args.dir_nextflow,
                             rtol=args.rtol, atol=args.atol)
    report = _render_markdown(results, args.dir_python, args.dir_nextflow)

    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(report)
        print(f"Report written to {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(report)

    # Exit code
    if any(not r.ok for r in results):
        sys.exit(1 if any(r.status == 'MISMATCH' for r in results) else 2)


if __name__ == '__main__':
    main()

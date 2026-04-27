"""
parallelization.py — Helpers for MDAnalysis 2.8+ native parallel backends.

Provides a small configuration layer that resolves user-facing options
(backend name, worker count) into keyword arguments suitable for
``AnalysisBase.run()``.

Supported backends
------------------
* ``'serial'``           — default, single-process execution.
* ``'multiprocessing'``  — stdlib ``multiprocessing`` (no extra deps).
* ``'dask'``             — Dask distributed (requires ``dask[distributed]``).

Only analyses whose MDAnalysis class explicitly supports non-serial
backends will use them.  RMSF, 2D-RMSD, and RoG always run serially.

Trajectory transformations
--------------------------
MDAnalysis checks ``TransformationBase.parallelizable`` at run time.
If any attached transformation is incompatible with the chosen backend,
``_configure_backend()`` raises ``ValueError``.  The helper
:func:`safe_run` catches this and retries with ``backend='serial'``.
"""

from __future__ import annotations

import multiprocessing
import os
import warnings
from dataclasses import dataclass, field
from typing import Any

# ─── Constants ────────────────────────────────────────────────────────────────

VALID_BACKENDS = ('serial', 'multiprocessing', 'dask')
"""Backends recognised by MDAnalysis 2.8+."""

# Analyses that support split-apply-combine parallelisation as of MDA 2.8.
# Others (RMSF, 2D-RMSD / DistanceMatrix) must always run serially.
PARALLELIZABLE_ANALYSES = frozenset({'RMSD', 'hbonds'})
"""Analysis keys whose MDAnalysis classes support non-serial backends."""

SERIAL_ONLY_ANALYSES = frozenset({'RMSF', '2D-RMSD', 'RoG'})
"""Analysis keys that must always run serially."""


# ─── Configuration dataclass ─────────────────────────────────────────────────

@dataclass
class ParallelConfig:
    """User-facing parallelisation settings.

    Attributes
    ----------
    backend : str
        One of ``'serial'``, ``'multiprocessing'``, or ``'dask'``.
    n_workers : int or None
        Number of worker processes.  ``None`` → auto-detect
        (``min(cpu_count, 4)`` for safety).
    """

    backend: str = 'serial'
    n_workers: int | None = None

    def __post_init__(self) -> None:
        self.backend = self.backend.strip().lower()
        if self.backend not in VALID_BACKENDS:
            raise ValueError(
                f"Unknown parallel backend '{self.backend}'. "
                f"Accepted values: {', '.join(VALID_BACKENDS)}"
            )
        if self.backend == 'dask':
            _check_dask_available()
        if self.n_workers is not None and self.n_workers < 1:
            raise ValueError(f"n_workers must be >= 1, got {self.n_workers}")

    @property
    def is_parallel(self) -> bool:
        """Return ``True`` when a non-serial backend is requested."""
        return self.backend != 'serial'


# ─── Public helpers ───────────────────────────────────────────────────────────

def resolve_n_workers(n_workers: int | None, backend: str) -> int:
    """Determine the effective number of workers.

    * For ``'serial'`` always returns ``1``.
    * When *n_workers* is ``None``, defaults to
      ``min(cpu_count, 4)`` to avoid saturating shared machines.
    * Respects the ``MD_ANALYSIS_N_WORKERS`` environment variable as
      an override.
    """
    if backend == 'serial':
        return 1

    # Environment variable override
    env_val = os.environ.get('MD_ANALYSIS_N_WORKERS')
    if env_val is not None:
        try:
            env_workers = int(env_val)
            if env_workers >= 1:
                return env_workers
        except ValueError:
            pass

    if n_workers is not None:
        return n_workers

    cpu_count = multiprocessing.cpu_count()
    return min(cpu_count, 4)


def get_run_kwargs(
    config: ParallelConfig,
    analysis_key: str,
) -> dict[str, Any]:
    """Build ``**kwargs`` for ``AnalysisBase.run()``.

    Parameters
    ----------
    config : ParallelConfig
        User-facing parallel settings.
    analysis_key : str
        One of the keys in :data:`PARALLELIZABLE_ANALYSES` or
        :data:`SERIAL_ONLY_ANALYSES`.

    Returns
    -------
    dict
        Contains ``backend`` and, for non-serial backends, ``n_workers``.
        Safe to unpack into ``.run(**kwargs)``.
    """
    # Serial-only analyses are never parallelised, regardless of config.
    if analysis_key not in PARALLELIZABLE_ANALYSES:
        if config.is_parallel:
            print(
                f"  INFO: {analysis_key} does not support parallel backends — "
                f"running serially."
            )
        return {'backend': 'serial'}

    effective_workers = resolve_n_workers(config.n_workers, config.backend)

    if config.backend == 'serial':
        return {'backend': 'serial'}

    kwargs: dict[str, Any] = {
        'backend': config.backend,
        'n_workers': effective_workers,
    }

    print(
        f"  Parallel: backend='{config.backend}', n_workers={effective_workers}"
    )
    # No verbose/progressbar in parallel mode — MDAnalysis disables it.
    return kwargs


def safe_run(analysis_obj, run_kwargs: dict[str, Any],
             *, start: int | None = None, stop: int | None = None,
             step: int | None = None) -> Any:
    """Run an MDAnalysis analysis with automatic serial fallback.

    Attempts ``analysis_obj.run(**run_kwargs)`` first.  If the backend
    raises ``ValueError`` (typically because trajectory transformations
    are not parallelisable), retries with ``backend='serial'``.

    Parameters
    ----------
    analysis_obj
        An ``AnalysisBase`` subclass instance (e.g. ``rms.RMSD``).
    run_kwargs : dict
        Keyword arguments for ``.run()`` (from :func:`get_run_kwargs`).
    start, stop, step : int or None
        Frame-slicing arguments forwarded to ``.run()``.

    Returns
    -------
    The analysis object (same as ``analysis_obj.run()`` returns).
    """
    # Merge frame-slicing into kwargs
    call_kw = dict(run_kwargs)
    if start is not None:
        call_kw['start'] = start
    if stop is not None:
        call_kw['stop'] = stop
    if step is not None:
        call_kw['step'] = step

    backend = call_kw.get('backend', 'serial')

    if backend == 'serial':
        return analysis_obj.run(**call_kw)

    try:
        return analysis_obj.run(**call_kw)
    except ValueError as exc:
        # MDAnalysis raises ValueError when transforms are not parallelisable
        if 'parallelizable' in str(exc).lower() or 'transform' in str(exc).lower():
            warnings.warn(
                f"Parallel backend '{backend}' failed due to non-parallelisable "
                f"trajectory transformations. Falling back to serial execution.\n"
                f"  Original error: {exc}",
                RuntimeWarning,
                stacklevel=2,
            )
            # Retry serially
            serial_kw = {k: v for k, v in call_kw.items()
                         if k not in ('backend', 'n_workers')}
            serial_kw['backend'] = 'serial'
            return analysis_obj.run(**serial_kw)
        raise  # Re-raise if it's a different ValueError


def _check_dask_available() -> None:
    """Raise ``ImportError`` with a helpful message if dask is not installed."""
    try:
        import dask  # noqa: F401
    except ImportError:
        raise ImportError(
            "The 'dask' parallel backend requires dask[distributed]. "
            "Install it with:  pip install 'dask[distributed]'  "
            "or:  uv sync --extra parallel"
        ) from None


def parallel_config_from_dict(d: dict) -> ParallelConfig:
    """Create a :class:`ParallelConfig` from a dict (e.g. parsed INI section).

    Recognised keys (all optional):
      * ``backend`` — str, default ``'serial'``
      * ``n_workers`` — int or ``None``, default ``None``
    """
    backend = d.get('backend', 'serial')
    n_workers_raw = d.get('n_workers')
    n_workers: int | None = None
    if n_workers_raw is not None:
        raw = str(n_workers_raw).strip().lower()
        if raw in ('none', 'auto', ''):
            n_workers = None
        else:
            n_workers = int(raw)
    return ParallelConfig(backend=backend, n_workers=n_workers)

"""
Tests for parallelization.py — ParallelConfig, backend resolution,
run_kwargs generation, and safe_run fallback logic.

Covers:
  - ParallelConfig validation (valid/invalid backends, n_workers)
  - resolve_n_workers defaults and environment variable override
  - get_run_kwargs for parallelizable vs serial-only analyses
  - safe_run with serial backend
  - safe_run fallback when parallel backend raises ValueError
  - parallel_config_from_dict parsing
  - dask availability check
"""
import os
import sys
import warnings
import pytest
from unittest.mock import MagicMock, patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from parallelization import (
    ParallelConfig,
    resolve_n_workers,
    get_run_kwargs,
    safe_run,
    parallel_config_from_dict,
    VALID_BACKENDS,
    PARALLELIZABLE_ANALYSES,
    SERIAL_ONLY_ANALYSES,
)


# ─── ParallelConfig ──────────────────────────────────────────────────────────

class TestParallelConfig:

    def test_default_serial(self):
        cfg = ParallelConfig()
        assert cfg.backend == 'serial'
        assert cfg.n_workers is None
        assert not cfg.is_parallel

    def test_multiprocessing(self):
        cfg = ParallelConfig(backend='multiprocessing', n_workers=4)
        assert cfg.backend == 'multiprocessing'
        assert cfg.n_workers == 4
        assert cfg.is_parallel

    def test_invalid_backend(self):
        with pytest.raises(ValueError, match="Unknown parallel backend"):
            ParallelConfig(backend='mpi')

    def test_invalid_n_workers(self):
        with pytest.raises(ValueError, match="n_workers must be >= 1"):
            ParallelConfig(backend='serial', n_workers=0)

    def test_negative_n_workers(self):
        with pytest.raises(ValueError, match="n_workers must be >= 1"):
            ParallelConfig(backend='serial', n_workers=-1)

    def test_backend_case_insensitive(self):
        cfg = ParallelConfig(backend='SERIAL')
        assert cfg.backend == 'serial'

    def test_backend_whitespace(self):
        cfg = ParallelConfig(backend='  multiprocessing  ')
        assert cfg.backend == 'multiprocessing'

    def test_dask_without_package(self):
        """Dask backend should raise ImportError if dask is not installed."""
        with patch.dict('sys.modules', {'dask': None}):
            with pytest.raises(ImportError, match="dask"):
                ParallelConfig(backend='dask')


# ─── resolve_n_workers ───────────────────────────────────────────────────────

class TestResolveNWorkers:

    def test_serial_always_1(self):
        assert resolve_n_workers(8, 'serial') == 1
        assert resolve_n_workers(None, 'serial') == 1

    def test_explicit_value(self):
        assert resolve_n_workers(4, 'multiprocessing') == 4

    def test_auto_capped_at_4(self):
        with patch('parallelization.multiprocessing.cpu_count', return_value=16):
            result = resolve_n_workers(None, 'multiprocessing')
            assert result == 4

    def test_auto_less_than_4(self):
        with patch('parallelization.multiprocessing.cpu_count', return_value=2):
            result = resolve_n_workers(None, 'multiprocessing')
            assert result == 2

    def test_env_override(self):
        with patch.dict(os.environ, {'MD_ANALYSIS_N_WORKERS': '7'}):
            result = resolve_n_workers(None, 'multiprocessing')
            assert result == 7

    def test_env_invalid_ignored(self):
        with patch.dict(os.environ, {'MD_ANALYSIS_N_WORKERS': 'invalid'}):
            with patch('parallelization.multiprocessing.cpu_count', return_value=8):
                result = resolve_n_workers(None, 'multiprocessing')
                assert result == 4  # falls through to auto


# ─── get_run_kwargs ──────────────────────────────────────────────────────────

class TestGetRunKwargs:

    def test_serial_always_serial(self):
        cfg = ParallelConfig(backend='serial')
        kwargs = get_run_kwargs(cfg, 'RMSD')
        assert kwargs == {'backend': 'serial'}

    def test_parallel_rmsd(self):
        cfg = ParallelConfig(backend='multiprocessing', n_workers=2)
        kwargs = get_run_kwargs(cfg, 'RMSD')
        assert kwargs['backend'] == 'multiprocessing'
        assert kwargs['n_workers'] == 2

    def test_parallel_hbonds(self):
        cfg = ParallelConfig(backend='multiprocessing', n_workers=3)
        kwargs = get_run_kwargs(cfg, 'hbonds')
        assert kwargs['backend'] == 'multiprocessing'
        assert kwargs['n_workers'] == 3

    def test_serial_only_rmsf(self):
        cfg = ParallelConfig(backend='multiprocessing', n_workers=4)
        kwargs = get_run_kwargs(cfg, 'RMSF')
        assert kwargs == {'backend': 'serial'}

    def test_serial_only_2d_rmsd(self):
        cfg = ParallelConfig(backend='multiprocessing', n_workers=4)
        kwargs = get_run_kwargs(cfg, '2D-RMSD')
        assert kwargs == {'backend': 'serial'}

    def test_serial_only_rog(self):
        cfg = ParallelConfig(backend='multiprocessing', n_workers=4)
        kwargs = get_run_kwargs(cfg, 'RoG')
        assert kwargs == {'backend': 'serial'}


# ─── safe_run ────────────────────────────────────────────────────────────────

class TestSafeRun:

    def test_serial_run(self):
        """Serial backend calls .run() directly."""
        mock_analysis = MagicMock()
        mock_analysis.run.return_value = mock_analysis

        result = safe_run(mock_analysis, {'backend': 'serial'}, start=10, step=1)

        mock_analysis.run.assert_called_once_with(
            backend='serial', start=10, step=1
        )
        assert result is mock_analysis

    def test_parallel_run_success(self):
        """Parallel backend succeeds without fallback."""
        mock_analysis = MagicMock()
        mock_analysis.run.return_value = mock_analysis

        result = safe_run(
            mock_analysis,
            {'backend': 'multiprocessing', 'n_workers': 2},
            start=5, stop=100
        )

        mock_analysis.run.assert_called_once_with(
            backend='multiprocessing', n_workers=2, start=5, stop=100
        )

    def test_parallel_fallback_on_transform_error(self):
        """Falls back to serial when transforms are not parallelizable."""
        mock_analysis = MagicMock()

        # First call raises ValueError about parallelizable transforms
        # Second call (serial) succeeds
        mock_analysis.run.side_effect = [
            ValueError("Transformation is not parallelizable"),
            mock_analysis,
        ]

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = safe_run(
                mock_analysis,
                {'backend': 'multiprocessing', 'n_workers': 2},
                start=5
            )

        # Should have been called twice — once parallel, once serial
        assert mock_analysis.run.call_count == 2
        second_call = mock_analysis.run.call_args
        assert second_call.kwargs.get('backend') == 'serial'
        assert 'n_workers' not in second_call.kwargs

        # Should have issued a warning
        assert len(w) == 1
        assert 'Falling back to serial' in str(w[0].message)

    def test_non_transform_valueerror_reraises(self):
        """Non-transform ValueError is re-raised, not caught."""
        mock_analysis = MagicMock()
        mock_analysis.run.side_effect = ValueError("Bad selection string")

        with pytest.raises(ValueError, match="Bad selection"):
            safe_run(
                mock_analysis,
                {'backend': 'multiprocessing', 'n_workers': 2}
            )


# ─── parallel_config_from_dict ───────────────────────────────────────────────

class TestParallelConfigFromDict:

    def test_empty_dict(self):
        cfg = parallel_config_from_dict({})
        assert cfg.backend == 'serial'
        assert cfg.n_workers is None

    def test_with_values(self):
        cfg = parallel_config_from_dict({'backend': 'multiprocessing', 'n_workers': '4'})
        assert cfg.backend == 'multiprocessing'
        assert cfg.n_workers == 4

    def test_none_string(self):
        cfg = parallel_config_from_dict({'backend': 'serial', 'n_workers': 'none'})
        assert cfg.n_workers is None

    def test_auto_string(self):
        cfg = parallel_config_from_dict({'backend': 'serial', 'n_workers': 'auto'})
        assert cfg.n_workers is None


# ─── Constants ───────────────────────────────────────────────────────────────

class TestConstants:

    def test_valid_backends(self):
        assert 'serial' in VALID_BACKENDS
        assert 'multiprocessing' in VALID_BACKENDS
        assert 'dask' in VALID_BACKENDS

    def test_parallelizable_analyses(self):
        assert 'RMSD' in PARALLELIZABLE_ANALYSES
        assert 'hbonds' in PARALLELIZABLE_ANALYSES

    def test_serial_only_analyses(self):
        assert 'RMSF' in SERIAL_ONLY_ANALYSES
        assert '2D-RMSD' in SERIAL_ONLY_ANALYSES
        assert 'RoG' in SERIAL_ONLY_ANALYSES

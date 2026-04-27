"""Tests for Python ↔ Nextflow output parity pipeline.

Covers:
- ini_to_nf_params: run_plots master switch, complete parameter mapping,
  NF param ↔ main.nf param name alignment
- compare_outputs: directory structure matching, _should_skip edge cases,
  multi-diff truncation
- compare_parallel_serial: tree comparison utilities
- CLI scripts: top_format passthrough to analysis functions
- Output structure: pickle naming, per-analysis subdirs, _ANALYSIS_SUBDIRS
  consistency with NF publishDir layout
"""
import json
import os
import pickle
import re
import sys
import textwrap

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ini_to_nf_params import _cfg_to_nf_params, _dump_yaml
from compare_outputs import (
    _collect_pickle_tree,
    _compare_values,
    _should_skip,
    _safe_pickle_load,
    _AnalysisUnpickler,
    compare_pickle_pair,
    run_comparison,
    _render_markdown,
    ComparisonResult,
)
from compare_parallel_serial import (
    compare_trees,
    _compare_values as cps_compare_values,
    ComparisonResult as CPSResult,
)
from executor import _ANALYSIS_SUBDIRS


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _make_cfg(**overrides):
    """Create a complete executor-style cfg dict for testing."""
    cfg = {
        'systems': ['SysA'],
        'variations': {'SysA': ['v1']},
        'num_replicates': 1,
        'traj_format': 'dcd',
        'top_format': 'top',
        'start_frame': 0,
        'input_dir': '/data',
        'output_dir': '/results',
        'run_rmsd': True, 'run_rmsf': True, 'run_2d_rmsd': False,
        'run_rog': True, 'run_hbonds': False,
        'run_plots': True,
        'plot_rmsd': True, 'plot_rmsf': True, 'plot_2d_rmsd': True,
        'plot_rog': True, 'plot_hbonds': True,
        'target_selection': 'protein and backbone',
        'ref_selection': ['protein and backbone'],
        'rog_selection': 'protein and backbone',
        'wrap_selection': 'auto',
        'time_interval_between_frames': 2.0,
        'time_unit': 'ns',
        'group_selections': None,
        'chain_intervals': None,
        'd_a_cutoff': 3.5,
        'd_h_a_angle_cutoff': 150.0,
        'update_selections': True,
        'acceptors_sel': None,
        'hydrogens_sel': None,
        'between_pairs': None,
        'dpi': 400,
        'hbonds_top_n': 20,
        'parallel_backend': 'serial',
        'n_workers': None,
    }
    cfg.update(overrides)
    return cfg


def _build_parity_tree(tmp_path, side, analysis_data):
    """Build a directory tree with pickle files for comparison testing.

    Parameters
    ----------
    side : str
        'py' or 'nf' (just used for subdir naming).
    analysis_data : dict[str, dict]
        Mapping of ``{analysis_subdir}/{filename}.pkl`` → pickle data.

    Returns the root directory path.
    """
    root = tmp_path / side
    for rel_path, data in analysis_data.items():
        full = root / rel_path
        full.parent.mkdir(parents=True, exist_ok=True)
        with open(full, 'wb') as f:
            pickle.dump(data, f)
    return str(root)


# ═══════════════════════════════════════════════════════════════════════════════
#  ini_to_nf_params — run_plots master switch
# ═══════════════════════════════════════════════════════════════════════════════

class TestRunPlotsMasterSwitch:
    """Verify that run_plots=false disables all plot flags in NF params."""

    def test_run_plots_false_disables_all_plot_flags(self):
        cfg = _make_cfg(run_plots=False, plot_rmsd=True, plot_rmsf=True,
                        plot_2d_rmsd=True, plot_rog=True, plot_hbonds=True)
        params = _cfg_to_nf_params(cfg)
        for key in ('plot_rmsd', 'plot_rmsf', 'plot_2d_rmsd', 'plot_rog', 'plot_hbonds'):
            assert params[key] is False, \
                f"{key} should be False when run_plots=False"

    def test_run_plots_true_preserves_individual_flags(self):
        cfg = _make_cfg(run_plots=True, plot_rmsd=True, plot_rmsf=False)
        params = _cfg_to_nf_params(cfg)
        assert params['plot_rmsd'] is True
        assert params['plot_rmsf'] is False

    def test_run_plots_default_true_preserves_flags(self):
        """When run_plots is absent (defaults to True), flags should pass through."""
        cfg = _make_cfg(plot_rmsd=True, plot_rmsf=False)
        # Simulate missing run_plots key (cfg.get('run_plots', True) → True)
        cfg.pop('run_plots', None)
        params = _cfg_to_nf_params(cfg)
        assert params['plot_rmsd'] is True
        assert params['plot_rmsf'] is False


# ═══════════════════════════════════════════════════════════════════════════════
#  ini_to_nf_params — complete parameter mapping
# ═══════════════════════════════════════════════════════════════════════════════

class TestCompleteParameterMapping:
    """Verify that _cfg_to_nf_params maps all parameters the NF pipeline needs."""

    def test_all_analysis_toggles_present(self):
        cfg = _make_cfg(run_rmsd=True, run_rmsf=False, run_2d_rmsd=True,
                        run_rog=False, run_hbonds=True)
        params = _cfg_to_nf_params(cfg)
        assert params['run_rmsd'] is True
        assert params['run_rmsf'] is False
        assert params['run_2d_rmsd'] is True
        assert params['run_rog'] is False
        assert params['run_hbonds'] is True

    def test_hbonds_params_mapped(self):
        cfg = _make_cfg(d_a_cutoff=4.0, d_h_a_angle_cutoff=120.0,
                        update_selections=False)
        params = _cfg_to_nf_params(cfg)
        assert params['d_a_cutoff'] == 4.0
        assert params['d_h_a_angle_cutoff'] == 120.0
        assert params['update_selections'] is False

    def test_acceptors_sel_mapped(self):
        cfg = _make_cfg(acceptors_sel='protein and name O*',
                        hydrogens_sel='nucleic and name H*')
        params = _cfg_to_nf_params(cfg)
        assert params['acceptors_sel'] == 'protein and name O*'
        assert params['hydrogens_sel'] == 'nucleic and name H*'

    def test_between_pairs_mapped_as_json(self):
        pairs = [['protein', 'resname URA']]
        cfg = _make_cfg(between_pairs=pairs)
        params = _cfg_to_nf_params(cfg)
        assert json.loads(params['between_pairs']) == pairs

    def test_chain_intervals_mapped_as_json(self):
        ci = {'SysA': {'A': [1, 100]}}
        cfg = _make_cfg(chain_intervals=ci)
        params = _cfg_to_nf_params(cfg)
        assert json.loads(params['chain_intervals']) == ci

    def test_hbonds_top_n_mapped(self):
        cfg = _make_cfg(hbonds_top_n=30)
        params = _cfg_to_nf_params(cfg)
        assert params['hbonds_top_n'] == 30

    def test_rog_selection_mapped(self):
        cfg = _make_cfg(rog_selection='protein')
        params = _cfg_to_nf_params(cfg)
        assert params['rog_selection'] == 'protein'

    def test_outdir_is_absolute(self):
        cfg = _make_cfg(output_dir='./relative_path')
        params = _cfg_to_nf_params(cfg)
        assert os.path.isabs(params['outdir'])

    def test_input_dir_is_absolute(self):
        cfg = _make_cfg(input_dir='./data')
        params = _cfg_to_nf_params(cfg)
        assert os.path.isabs(params['input_dir'])


# ═══════════════════════════════════════════════════════════════════════════════
#  ini_to_nf_params — param names match main.nf declarations
# ═══════════════════════════════════════════════════════════════════════════════

class TestNfParamNameAlignment:
    """Verify that every key from _cfg_to_nf_params matches a params.X in main.nf."""

    @pytest.fixture
    def nf_param_names(self):
        """Extract all params.X names declared in main.nf."""
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        nf_path = os.path.join(project_root, 'main.nf')
        with open(nf_path, 'r') as f:
            content = f.read()
        # Match params.foo declarations (top-level param defaults)
        names = set(re.findall(r'params\.(\w+)\s*=', content))
        return names

    def test_all_mapped_params_exist_in_nf(self, nf_param_names):
        """Every key produced by _cfg_to_nf_params should be a valid NF param."""
        # Build a maximally populated cfg
        cfg = _make_cfg(
            group_selections=['protein', 'nucleic'],
            chain_intervals={'SysA': {'A': [1, 100]}},
            time_interval_between_frames=2.0,
            acceptors_sel='protein and name O*',
            hydrogens_sel='nucleic and name H*',
            between_pairs=[['protein', 'nucleic']],
            parallel_backend='multiprocessing',
            n_workers=4,
        )
        params = _cfg_to_nf_params(cfg)

        unrecognized = set(params.keys()) - nf_param_names
        assert unrecognized == set(), \
            f"These params are produced by _cfg_to_nf_params but not declared in main.nf: {unrecognized}"


# ═══════════════════════════════════════════════════════════════════════════════
#  Output structure: _ANALYSIS_SUBDIRS ↔ NF publishDir consistency
# ═══════════════════════════════════════════════════════════════════════════════

class TestAnalysisSubdirConsistency:
    """Verify executor subdirs match Nextflow publishDir layout."""

    @pytest.fixture
    def nf_publish_dirs(self):
        """Extract publishDir analysis subdirectory names from main.nf."""
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        nf_path = os.path.join(project_root, 'main.nf')
        with open(nf_path, 'r') as f:
            content = f.read()
        # publishDir "${params.outdir}/rmsd" → extract 'rmsd'
        dirs = set(re.findall(
            r'publishDir\s+"\$\{params\.outdir\}/(\w+)"', content))
        return dirs

    def test_executor_subdirs_match_nf(self, nf_publish_dirs):
        """Every NF publishDir for analyses should exist in _ANALYSIS_SUBDIRS."""
        executor_dirs = set(_ANALYSIS_SUBDIRS.values())
        # NF also has plots/ subdirs — filter those out
        nf_analysis_dirs = {d for d in nf_publish_dirs if d != 'plots'}
        missing_in_executor = nf_analysis_dirs - executor_dirs
        assert missing_in_executor == set(), \
            f"NF has publishDir for {missing_in_executor} but executor has no matching subdir"

    def test_executor_subdirs_covered_by_nf(self, nf_publish_dirs):
        """Every executor analysis subdir should have a matching NF publishDir."""
        executor_dirs = set(_ANALYSIS_SUBDIRS.values())
        nf_analysis_dirs = {d for d in nf_publish_dirs if d != 'plots'}
        missing_in_nf = executor_dirs - nf_analysis_dirs
        assert missing_in_nf == set(), \
            f"Executor has {missing_in_nf} but NF has no matching publishDir"


# ═══════════════════════════════════════════════════════════════════════════════
#  compare_outputs — _should_skip edge cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestShouldSkip:
    """Cover all _should_skip patterns."""

    def test_skips_non_pkl(self):
        assert _should_skip('rmsd/data.csv') is True

    def test_skips_pipeline_info(self):
        assert _should_skip('pipeline_info/trace.pkl') is True

    def test_skips_plots(self):
        assert _should_skip('plots/rmsd/plot.pkl') is True

    def test_skips_rmsfit(self):
        assert _should_skip('rmsfit_aligned.pkl') is True

    def test_skips_nextflow_dir(self):
        assert _should_skip('.nextflow/cache.pkl') is True

    def test_skips_pycache(self):
        assert _should_skip('__pycache__/module.pkl') is True

    def test_accepts_analysis_pkl(self):
        assert _should_skip('rmsd/rmsd_plot_SysA_v1_rep1.pkl') is False

    def test_accepts_nested_analysis_pkl(self):
        assert _should_skip('rog/rog_plot_SysA_v1_rep1.pkl') is False


# ═══════════════════════════════════════════════════════════════════════════════
#  compare_outputs — multi-diff truncation
# ═══════════════════════════════════════════════════════════════════════════════

class TestMultiDiffTruncation:
    """compare_pickle_pair should truncate details when >5 diffs."""

    def test_more_than_5_diffs_truncated(self, tmp_path):
        a_data = {f'key_{i}': float(i) for i in range(10)}
        b_data = {f'key_{i}': float(i * 100) for i in range(10)}
        a = tmp_path / 'a.pkl'
        b = tmp_path / 'b.pkl'
        a.write_bytes(pickle.dumps(a_data))
        b.write_bytes(pickle.dumps(b_data))
        result = compare_pickle_pair(str(a), str(b), 'test.pkl')
        assert result.status == 'MISMATCH'
        assert 'and' in result.details and 'more' in result.details


# ═══════════════════════════════════════════════════════════════════════════════
#  compare_outputs — full parity tree matching
# ═══════════════════════════════════════════════════════════════════════════════

class TestParityTreeMatching:
    """End-to-end tests for compare_outputs.run_comparison with realistic trees."""

    def test_matching_multi_analysis_tree(self, tmp_path):
        """Identical trees with multiple analysis types should all MATCH."""
        data_rmsd = {'rmsd_data': np.array([1.0, 2.0, 3.0])}
        data_rmsf = {'rmsf_obj': np.array([0.5, 0.6])}
        data_rog = {'rog_obj': np.array([15.0, 15.1])}
        files = {
            'rmsd/rmsd_plot_SysA_v1_rep1.pkl': data_rmsd,
            'rmsf/rmsf_plot_SysA_v1_rep1.pkl': data_rmsf,
            'rog/rog_plot_SysA_v1_rep1.pkl': data_rog,
        }
        py_dir = _build_parity_tree(tmp_path, 'py', files)
        nf_dir = _build_parity_tree(tmp_path, 'nf', files)
        results = run_comparison(py_dir, nf_dir)
        assert len(results) == 3
        assert all(r.ok for r in results)

    def test_missing_analysis_in_one_side(self, tmp_path):
        """When one side has an extra analysis, it shows as structural mismatch."""
        py_dir = _build_parity_tree(tmp_path, 'py', {
            'rmsd/rmsd_plot_SysA_v1_rep1.pkl': {'v': 1},
            'rog/rog_plot_SysA_v1_rep1.pkl': {'v': 2},
        })
        nf_dir = _build_parity_tree(tmp_path, 'nf', {
            'rmsd/rmsd_plot_SysA_v1_rep1.pkl': {'v': 1},
        })
        results = run_comparison(py_dir, nf_dir)
        only_py = [r for r in results if r.status == 'ONLY_PYTHON']
        assert len(only_py) == 1
        assert 'rog' in only_py[0].filename

    def test_numeric_tolerance_within_epsilon(self, tmp_path):
        """Values within rtol/atol should MATCH."""
        py_dir = _build_parity_tree(tmp_path, 'py', {
            'rmsd/test.pkl': {'v': np.array([1.0, 2.0, 3.0])},
        })
        nf_dir = _build_parity_tree(tmp_path, 'nf', {
            'rmsd/test.pkl': {'v': np.array([1.0, 2.0, 3.0 + 1e-7])},
        })
        results = run_comparison(py_dir, nf_dir)
        assert all(r.ok for r in results)

    def test_pipeline_info_skipped(self, tmp_path):
        """pipeline_info pickles should be ignored during comparison."""
        py_dir = _build_parity_tree(tmp_path, 'py', {
            'rmsd/test.pkl': {'v': 1},
        })
        nf_root = tmp_path / 'nf'
        # Create both a real pickle and a pipeline_info one
        rmsd_dir = nf_root / 'rmsd'
        rmsd_dir.mkdir(parents=True)
        with open(rmsd_dir / 'test.pkl', 'wb') as f:
            pickle.dump({'v': 1}, f)
        pi_dir = nf_root / 'pipeline_info'
        pi_dir.mkdir()
        with open(pi_dir / 'spurious.pkl', 'wb') as f:
            pickle.dump({'should': 'be ignored'}, f)

        results = run_comparison(py_dir, str(nf_root))
        assert len(results) == 1  # Only rmsd/test.pkl
        assert results[0].ok


# ═══════════════════════════════════════════════════════════════════════════════
#  compare_parallel_serial — string/object array handling
# ═══════════════════════════════════════════════════════════════════════════════

class TestCPSCompareValues:
    """Cover compare_parallel_serial._compare_values edge cases."""

    def test_string_arrays_equal(self):
        a = np.array(['A', 'B', 'C'])
        b = np.array(['A', 'B', 'C'])
        diffs = cps_compare_values(a, b)
        assert diffs == []

    def test_string_arrays_differ(self):
        a = np.array(['A', 'B', 'C'])
        b = np.array(['A', 'X', 'C'])
        diffs = cps_compare_values(a, b)
        assert len(diffs) == 1
        assert 'string/object array' in diffs[0]

    def test_objects_with_dict_compared(self):
        class Obj:
            def __init__(self, x):
                self.x = x
        diffs = cps_compare_values(Obj(1), Obj(1))
        assert diffs == []

    def test_integer_numpy_scalar_mismatch(self):
        diffs = cps_compare_values(np.int64(5), np.int64(10))
        assert len(diffs) == 1


# ═══════════════════════════════════════════════════════════════════════════════
#  CLI scripts: top_format passthrough
# ═══════════════════════════════════════════════════════════════════════════════

class TestTopFormatPassthrough:
    """Verify that --top-format is forwarded to the analysis function."""

    def test_rms_main_passes_top_format(self, monkeypatch):
        """run_rms_analysis.main() should pass top_format to run_rms_analysis()."""
        from run_rms_analysis import main as rms_main
        captured = {}

        def fake_run_rms_analysis(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr('run_rms_analysis.run_rms_analysis', fake_run_rms_analysis)
        monkeypatch.setattr('sys.argv', [
            'run_rms_analysis',
            '--systems', '["SysA"]',
            '--variations', '{"SysA": ["v1"]}',
            '--analysis', 'RMSD',
            '--target-selection', 'protein',
            '--ref-selection', 'protein',
            '--top-format', 'psf',
        ])
        rms_main()
        assert captured.get('top_format') == 'psf'

    def test_rog_main_passes_top_format(self, monkeypatch):
        """run_rog_analysis.main() should pass top_format to run_rog_analysis()."""
        from run_rog_analysis import main as rog_main
        captured = {}

        def fake_run_rog_analysis(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr('run_rog_analysis.run_rog_analysis', fake_run_rog_analysis)
        monkeypatch.setattr('sys.argv', [
            'run_rog_analysis',
            '--systems', '["SysA"]',
            '--variations', '{"SysA": ["v1"]}',
            '--top-format', 'psf',
        ])
        rog_main()
        assert captured.get('top_format') == 'psf'

    def test_hbonds_main_passes_top_format(self, monkeypatch):
        """run_hbonds_analysis.main() should pass top_format to run_hbonds_analysis()."""
        from run_hbonds_analysis import main as hbonds_main
        captured = {}

        def fake_run_hbonds_analysis(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr('run_hbonds_analysis.run_hbonds_analysis', fake_run_hbonds_analysis)
        monkeypatch.setattr('sys.argv', [
            'run_hbonds_analysis',
            '--systems', '["SysA"]',
            '--variations', '{"SysA": ["v1"]}',
            '--between-pairs', '[["protein", "nucleic"]]',
            '--top-format', 'psf',
        ])
        hbonds_main()
        assert captured.get('top_format') == 'psf'


# ═══════════════════════════════════════════════════════════════════════════════
#  Pickle naming: executor vs Nextflow produce same filenames
# ═══════════════════════════════════════════════════════════════════════════════

class TestPickleNamingConvention:
    """Verify that both modes would produce identically named pickle files."""

    def test_rmsd_pickle_name(self):
        # Executor: run_rms_analysis writes rmsd_plot_{sys}_{var}_rep{N}.pkl
        # NF: writes the same then renames from rep1 → rep{N}
        expected = 'rmsd_plot_SysA_wild_rep1.pkl'
        assert re.match(r'rmsd_plot_\w+_\w+_rep\d+\.pkl', expected)

    def test_rmsf_pickle_name(self):
        expected = 'rmsf_plot_SysA_wild_rep1.pkl'
        assert re.match(r'rmsf_plot_\w+_\w+_rep\d+\.pkl', expected)

    def test_rmsf_chain_pickle_name(self):
        expected = 'rmsf_plot_SysA_wild_rep1_chainA.pkl'
        assert re.match(r'rmsf_plot_\w+_\w+_rep\d+_chain\w+\.pkl', expected)

    def test_2d_rmsd_pickle_name(self):
        expected = '2d_rmsd_plot_SysA_wild_rep1.pkl'
        assert re.match(r'2d_rmsd_plot_\w+_\w+_rep\d+\.pkl', expected)

    def test_rog_pickle_name(self):
        expected = 'rog_plot_SysA_wild_rep1.pkl'
        assert re.match(r'rog_plot_\w+_\w+_rep\d+\.pkl', expected)

    def test_hbonds_pickle_name(self):
        expected = 'hbonds_plot_SysA_wild_rep1.pkl'
        assert re.match(r'hbonds_plot_\w+_\w+_rep\d+\.pkl', expected)


# ═══════════════════════════════════════════════════════════════════════════════
#  _render_markdown — empty and edge cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestRenderMarkdownEdgeCases:
    def test_all_statuses_counted(self):
        results = [
            ComparisonResult('a.pkl', 'MATCH', 'ok'),
            ComparisonResult('b.pkl', 'MISMATCH', 'bad'),
            ComparisonResult('c.pkl', 'ONLY_PYTHON', 'left'),
            ComparisonResult('d.pkl', 'ONLY_NEXTFLOW', 'right'),
        ]
        md = _render_markdown(results, '/py', '/nf')
        assert '1 match' in md
        assert '1 mismatch' in md
        assert '1 only-python' in md
        assert '1 only-nextflow' in md

    def test_sorted_by_filename(self):
        results = [
            ComparisonResult('z.pkl', 'MATCH', ''),
            ComparisonResult('a.pkl', 'MATCH', ''),
        ]
        md = _render_markdown(results, '/py', '/nf')
        a_pos = md.index('a.pkl')
        z_pos = md.index('z.pkl')
        assert a_pos < z_pos


# ═══════════════════════════════════════════════════════════════════════════════
#  compare_outputs — corrupt first pickle file
# ═══════════════════════════════════════════════════════════════════════════════

class TestCorruptFirstPickle:
    def test_corrupt_first_file(self, tmp_path):
        """compare_pickle_pair should report MISMATCH if first file is corrupt."""
        a = tmp_path / 'a.pkl'
        b = tmp_path / 'b.pkl'
        a.write_bytes(b'not-a-pickle')
        b.write_bytes(pickle.dumps({'ok': True}))
        result = compare_pickle_pair(str(a), str(b), 'test.pkl')
        assert result.status == 'MISMATCH'
        assert 'Failed to load Python pickle' in result.details


# ═══════════════════════════════════════════════════════════════════════════════
#  compare_parallel_serial — tree comparison with realistic structures
# ═══════════════════════════════════════════════════════════════════════════════

class TestCPSTreeComparison:
    """compare_trees with per-analysis subdirectory structure."""

    def test_matching_trees_with_subdirs(self, tmp_path):
        data = {'v': np.array([1.0, 2.0])}
        serial = _build_parity_tree(tmp_path, 'serial', {
            'rmsd/rmsd_plot_SysA_v1_rep1.pkl': data,
            'rog/rog_plot_SysA_v1_rep1.pkl': data,
        })
        parallel = _build_parity_tree(tmp_path, 'parallel', {
            'rmsd/rmsd_plot_SysA_v1_rep1.pkl': data,
            'rog/rog_plot_SysA_v1_rep1.pkl': data,
        })
        results = compare_trees(serial, parallel)
        assert len(results) == 2
        assert all(r.ok for r in results)

    def test_only_serial_detected(self, tmp_path):
        serial = _build_parity_tree(tmp_path, 'serial', {
            'rmsd/test.pkl': {'v': 1},
        })
        parallel = _build_parity_tree(tmp_path, 'parallel', {})
        (tmp_path / 'parallel').mkdir(exist_ok=True)
        results = compare_trees(serial, str(tmp_path / 'parallel'))
        assert any(r.status == 'ONLY_SERIAL' for r in results)

    def test_only_parallel_detected(self, tmp_path):
        serial = _build_parity_tree(tmp_path, 'serial', {})
        (tmp_path / 'serial').mkdir(exist_ok=True)
        parallel = _build_parity_tree(tmp_path, 'parallel', {
            'rmsd/test.pkl': {'v': 1},
        })
        results = compare_trees(str(tmp_path / 'serial'), parallel)
        assert any(r.status == 'ONLY_PARALLEL' for r in results)


# ═══════════════════════════════════════════════════════════════════════════════
#  RoGResults.__module__ fix — pickle portability
# ═══════════════════════════════════════════════════════════════════════════════

class TestRoGResultsModuleAttribute:
    """Verify RoGResults.__module__ is set to 'run_rog_analysis' so pickles
    serialized under __main__ (Nextflow) can be loaded by any script."""

    def test_module_attribute_is_canonical(self):
        from run_rog_analysis import RoGResults
        assert RoGResults.__module__ == 'run_rog_analysis', \
            "RoGResults.__module__ must be 'run_rog_analysis' for cross-script pickle compat"

    def test_pickle_roundtrip_stores_canonical_module(self, tmp_path):
        """A pickled RoGResults should reference 'run_rog_analysis', not '__main__'."""
        from run_rog_analysis import RoGResults
        obj = RoGResults(
            frames=np.array([0, 1, 2]),
            times=np.array([0.0, 0.001, 0.002]),
            rog_values=np.array([15.0, 15.1, 15.2]),
        )
        pkl_path = tmp_path / 'rog_test.pkl'
        with open(pkl_path, 'wb') as f:
            pickle.dump(obj, f)

        # Read raw bytes and check the module reference
        raw = pkl_path.read_bytes()
        assert b'run_rog_analysis' in raw
        assert b'__main__' not in raw


# ═══════════════════════════════════════════════════════════════════════════════
#  _safe_pickle_load / _AnalysisUnpickler — __main__ redirection
# ═══════════════════════════════════════════════════════════════════════════════

def _dump_as_main(obj, path, cls):
    """Pickle *obj* so its class is stored under ``__main__`` — simulating
    the scenario when Nextflow runs an analysis script directly.

    ``pickle.dump`` verifies the class is accessible in the module referenced
    by ``cls.__module__``, so we must temporarily inject the class into
    ``sys.modules['__main__']`` as well.
    """
    import sys
    main_mod = sys.modules['__main__']
    original_module = cls.__module__
    had_attr = hasattr(main_mod, cls.__name__)
    try:
        cls.__module__ = '__main__'
        setattr(main_mod, cls.__name__, cls)
        with open(path, 'wb') as f:
            pickle.dump(obj, f)
    finally:
        cls.__module__ = original_module
        if not had_attr:
            delattr(main_mod, cls.__name__)


class TestSafePickleLoad:
    """Verify _safe_pickle_load handles pickles with __main__ class references."""

    def test_loads_normal_pickle(self, tmp_path):
        """Standard pickles load without issues."""
        data = {'key': np.array([1.0, 2.0])}
        pkl = tmp_path / 'normal.pkl'
        pkl.write_bytes(pickle.dumps(data))
        result = _safe_pickle_load(str(pkl))
        assert 'key' in result
        assert np.allclose(result['key'], [1.0, 2.0])

    def test_loads_pickle_with_main_module_rog(self, tmp_path):
        """Simulate a Nextflow-produced RoG pickle that references __main__.RoGResults."""
        from run_rog_analysis import RoGResults
        obj = RoGResults(
            frames=np.array([0, 1]),
            times=np.array([0.0, 0.001]),
            rog_values=np.array([15.0, 15.1]),
        )
        pkl = tmp_path / 'rog_main.pkl'
        _dump_as_main(obj, pkl, RoGResults)

        # Verify the raw pickle contains __main__
        raw = pkl.read_bytes()
        assert b'__main__' in raw

        # _safe_pickle_load should still succeed
        loaded = _safe_pickle_load(str(pkl))
        assert hasattr(loaded, 'rog_data')
        assert loaded.rog_data.shape == (2, 3)
        assert np.allclose(loaded.rog_data[:, 2], [15.0, 15.1])

    def test_compare_pickle_pair_handles_main_module(self, tmp_path):
        """compare_pickle_pair should succeed even when one pickle uses __main__."""
        from run_rog_analysis import RoGResults
        obj = RoGResults(
            frames=np.array([0, 1, 2]),
            times=np.array([0.0, 0.001, 0.002]),
            rog_values=np.array([15.0, 15.1, 15.2]),
        )

        # "Python" pickle — uses canonical module name
        py_pkl = tmp_path / 'py.pkl'
        with open(py_pkl, 'wb') as f:
            pickle.dump(obj, f)

        # "NF" pickle — simulate __main__ reference
        nf_pkl = tmp_path / 'nf.pkl'
        _dump_as_main(obj, nf_pkl, RoGResults)

        result = compare_pickle_pair(str(py_pkl), str(nf_pkl), 'rog_test.pkl')
        assert result.status == 'MATCH', f"Expected MATCH but got {result.status}: {result.details}"

    def test_corrupt_pickle_still_fails(self, tmp_path):
        """Truly corrupt pickles should still produce MISMATCH."""
        good = tmp_path / 'good.pkl'
        good.write_bytes(pickle.dumps({'ok': True}))
        bad = tmp_path / 'bad.pkl'
        bad.write_bytes(b'definitely-not-a-pickle')
        result = compare_pickle_pair(str(bad), str(good), 'test.pkl')
        assert result.status == 'MISMATCH'


class TestAnalysisUnpickler:
    """Unit tests for the _AnalysisUnpickler class."""

    def test_redirects_main_to_run_rog_analysis(self, tmp_path):
        """__main__.RoGResults should resolve to run_rog_analysis.RoGResults."""
        from run_rog_analysis import RoGResults
        obj = RoGResults(
            frames=np.array([0]),
            times=np.array([0.0]),
            rog_values=np.array([15.0]),
        )
        pkl = tmp_path / 'test.pkl'
        _dump_as_main(obj, pkl, RoGResults)

        with open(pkl, 'rb') as f:
            loaded = _AnalysisUnpickler(f).load()
        assert type(loaded).__name__ == 'RoGResults'
        assert type(loaded).__module__ == 'run_rog_analysis'

    def test_non_main_module_unaffected(self, tmp_path):
        """Pickles with normal module refs should load as usual."""
        data = {'normal': 42}
        pkl = tmp_path / 'normal.pkl'
        pkl.write_bytes(pickle.dumps(data))
        with open(pkl, 'rb') as f:
            loaded = _AnalysisUnpickler(f).load()
        assert loaded == data

"""Tests for ini_to_nf_params.py — INI to Nextflow YAML converter."""
import os
import sys
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ini_to_nf_params import _cfg_to_nf_params, _dump_yaml


def _make_cfg(**overrides):
    """Create a minimal executor-style cfg dict."""
    cfg = {
        'systems': ['SysA'],
        'variations': {'SysA': ['v1']},
        'num_replicates': 1,
        'traj_format': 'dcd',
        'top_format': 'top',
        'start_frame': 0,
        'input_dir': '/data',
        'output_dir': '/results',
        'run_rmsd': True, 'run_rmsf': False, 'run_2d_rmsd': False,
        'run_rog': True, 'run_hbonds': False,
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
    }
    cfg.update(overrides)
    return cfg


class TestCfgToNfParams:
    def test_basic_mapping(self):
        cfg = _make_cfg()
        params = _cfg_to_nf_params(cfg)
        assert params['num_replicates'] == 1
        assert params['run_rmsd'] is True
        assert params['run_rmsf'] is False
        assert params['traj_format'] == 'dcd'

    def test_systems_serialized_as_json(self):
        cfg = _make_cfg()
        params = _cfg_to_nf_params(cfg)
        assert '"SysA"' in params['systems']

    def test_wrap_selection_none_becomes_none_str(self):
        cfg = _make_cfg(wrap_selection=None)
        params = _cfg_to_nf_params(cfg)
        assert params['wrap_selection'] == 'none'

    def test_ref_selection_takes_first(self):
        cfg = _make_cfg(ref_selection=['sel_a', 'sel_b'])
        params = _cfg_to_nf_params(cfg)
        assert params['ref_selection'] == 'sel_a'

    def test_omits_none_time_interval(self):
        cfg = _make_cfg(time_interval_between_frames=None)
        params = _cfg_to_nf_params(cfg)
        assert 'time_interval_between_frames' not in params

    def test_includes_group_selections_when_set(self):
        cfg = _make_cfg(group_selections=['protein', 'resname URA'])
        params = _cfg_to_nf_params(cfg)
        assert 'group_selections' in params

    def test_parallel_backend_default(self):
        """parallel_backend should default to 'serial' when absent."""
        cfg = _make_cfg()
        params = _cfg_to_nf_params(cfg)
        assert params['parallel_backend'] == 'serial'

    def test_parallel_backend_explicit(self):
        """Explicit parallel_backend should be passed through."""
        cfg = _make_cfg(parallel_backend='multiprocessing')
        params = _cfg_to_nf_params(cfg)
        assert params['parallel_backend'] == 'multiprocessing'

    def test_n_workers_included_when_set(self):
        """n_workers should be included when not None."""
        cfg = _make_cfg(n_workers=4)
        params = _cfg_to_nf_params(cfg)
        assert params['n_workers'] == 4

    def test_n_workers_omitted_when_none(self):
        """n_workers should be omitted when None."""
        cfg = _make_cfg(n_workers=None)
        params = _cfg_to_nf_params(cfg)
        assert 'n_workers' not in params


class TestDumpYaml:
    def test_bool_lowercase(self):
        yaml = _dump_yaml({'run_rmsd': True, 'run_rmsf': False})
        assert 'run_rmsd: true' in yaml
        assert 'run_rmsf: false' in yaml

    def test_string_quoting(self):
        yaml = _dump_yaml({'sel': 'protein and backbone'})
        assert "'protein and backbone'" in yaml

    def test_numeric_unquoted(self):
        yaml = _dump_yaml({'dpi': 400, 'cutoff': 3.5})
        assert 'dpi: 400' in yaml
        assert 'cutoff: 3.5' in yaml

    def test_null_value(self):
        yaml = _dump_yaml({'x': None})
        assert 'x: null' in yaml

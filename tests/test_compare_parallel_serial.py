"""
Tests for compare_parallel_serial.py — pickle comparison engine.

Covers:
  - _collect_pickle_tree file discovery
  - _compare_values for arrays, dicts, scalars, objects
  - compare_pickle_pair for matched pairs
  - compare_trees for directory comparison
  - format_report output structure
  - ComparisonResult properties
  - _create_temp_ini INI override generation
  - main() CLI mode validation
"""
import os
import sys
import pickle
import pytest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from compare_parallel_serial import (
    _collect_pickle_tree,
    _compare_values,
    _create_temp_ini,
    compare_pickle_pair,
    compare_trees,
    format_report,
    main,
    ComparisonResult,
)


# ─── _collect_pickle_tree ────────────────────────────────────────────────────

class TestCollectPickleTree:

    def test_finds_pkl_files(self, tmp_path):
        (tmp_path / 'rmsd_plot.pkl').write_bytes(b'data')
        (tmp_path / 'sub').mkdir()
        (tmp_path / 'sub' / 'rmsf_plot.pkl').write_bytes(b'data')
        tree = _collect_pickle_tree(str(tmp_path))
        assert 'rmsd_plot.pkl' in tree
        assert os.path.join('sub', 'rmsf_plot.pkl') in tree

    def test_ignores_non_pkl(self, tmp_path):
        (tmp_path / 'data.csv').write_text('a,b')
        (tmp_path / 'rmsd.pkl').write_bytes(b'data')
        tree = _collect_pickle_tree(str(tmp_path))
        assert len(tree) == 1
        assert 'rmsd.pkl' in tree

    def test_skips_pycache(self, tmp_path):
        cache = tmp_path / '__pycache__'
        cache.mkdir()
        (cache / 'module.pkl').write_bytes(b'data')
        tree = _collect_pickle_tree(str(tmp_path))
        assert len(tree) == 0

    def test_skips_pipeline_info(self, tmp_path):
        info = tmp_path / 'pipeline_info'
        info.mkdir()
        (info / 'trace.pkl').write_bytes(b'data')
        tree = _collect_pickle_tree(str(tmp_path))
        assert len(tree) == 0

    def test_empty_dir(self, tmp_path):
        tree = _collect_pickle_tree(str(tmp_path))
        assert tree == {}


# ─── _compare_values ─────────────────────────────────────────────────────────

class TestCompareValues:

    def test_identical_arrays(self):
        a = np.array([1.0, 2.0, 3.0])
        diffs = _compare_values(a, a.copy())
        assert diffs == []

    def test_close_arrays(self):
        a = np.array([1.0, 2.0, 3.0])
        b = a + 1e-10
        diffs = _compare_values(a, b, rtol=1e-5, atol=1e-8)
        assert diffs == []

    def test_different_arrays(self):
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([1.0, 2.0, 4.0])
        diffs = _compare_values(a, b)
        assert len(diffs) == 1
        assert 'differ' in diffs[0]

    def test_shape_mismatch(self):
        a = np.array([1.0, 2.0])
        b = np.array([1.0, 2.0, 3.0])
        diffs = _compare_values(a, b)
        assert len(diffs) == 1
        assert 'shape' in diffs[0]

    def test_dict_identical(self):
        d = {'a': 1, 'b': 'hello'}
        diffs = _compare_values(d, d.copy())
        assert diffs == []

    def test_dict_missing_key(self):
        a = {'x': 1}
        b = {'x': 1, 'y': 2}
        diffs = _compare_values(a, b)
        assert len(diffs) == 1
        assert 'only in' in diffs[0]

    def test_list_identical(self):
        diffs = _compare_values([1, 2, 3], [1, 2, 3])
        assert diffs == []

    def test_list_length_mismatch(self):
        diffs = _compare_values([1, 2], [1, 2, 3])
        assert len(diffs) == 1
        assert 'length' in diffs[0]

    def test_scalar_mismatch(self):
        diffs = _compare_values(1.0, 2.0)
        assert len(diffs) == 1

    def test_string_mismatch(self):
        diffs = _compare_values('hello', 'world')
        assert len(diffs) == 1

    def test_type_mismatch(self):
        diffs = _compare_values(1, 'one')
        assert len(diffs) == 1
        assert 'type mismatch' in diffs[0]

    def test_nan_handling(self):
        a = np.array([1.0, np.nan, 3.0])
        b = np.array([1.0, np.nan, 3.0])
        diffs = _compare_values(a, b)
        assert diffs == []

    def test_object_comparison(self):
        class Obj:
            def __init__(self, val):
                self.val = val
        diffs = _compare_values(Obj(1), Obj(1))
        assert diffs == []
        diffs = _compare_values(Obj(1), Obj(2))
        assert len(diffs) == 1

    def test_circular_reference_does_not_recurse_infinitely(self):
        """Objects with circular __dict__ references must not cause RecursionError."""
        class CircularObj:
            pass
        a = CircularObj()
        b = CircularObj()
        # Create circular reference: a.child.parent -> a
        a.child = CircularObj()
        a.child.parent = a
        a.child.value = 42
        b.child = CircularObj()
        b.child.parent = b
        b.child.value = 42
        # Should complete without RecursionError
        diffs = _compare_values(a, b)
        assert isinstance(diffs, list)

    def test_deeply_nested_objects_respect_depth_limit(self):
        """Deeply nested structures beyond _MAX_COMPARE_DEPTH are treated as equal."""
        from compare_parallel_serial import _MAX_COMPARE_DEPTH
        # Build a chain deeper than _MAX_COMPARE_DEPTH
        class Node:
            pass
        a_root = Node()
        b_root = Node()
        a_cur, b_cur = a_root, b_root
        for i in range(_MAX_COMPARE_DEPTH + 10):
            a_next, b_next = Node(), Node()
            a_cur.child = a_next
            b_cur.child = b_next
            a_cur, b_cur = a_next, b_next
        # Leaf values differ, but they're beyond the depth limit
        a_cur.value = 1
        b_cur.value = 999
        diffs = _compare_values(a_root, b_root)
        # The difference at the leaf should NOT be reported because it's beyond depth
        assert all('value' not in d for d in diffs)

    def test_dict_circular_reference(self):
        """Dicts with circular references must not cause RecursionError."""
        a = {'key': 'val'}
        b = {'key': 'val'}
        a['self'] = a  # circular
        b['self'] = b
        diffs = _compare_values(a, b)
        assert isinstance(diffs, list)


# ─── ComparisonResult ────────────────────────────────────────────────────────

class TestComparisonResult:

    def test_match_is_ok(self):
        r = ComparisonResult('test.pkl', 'MATCH')
        assert r.ok

    def test_mismatch_not_ok(self):
        r = ComparisonResult('test.pkl', 'MISMATCH', 'some diff')
        assert not r.ok

    def test_only_serial_not_ok(self):
        r = ComparisonResult('test.pkl', 'ONLY_SERIAL')
        assert not r.ok


# ─── compare_pickle_pair ─────────────────────────────────────────────────────

class TestComparePicklePair:

    def test_matching_pair(self, tmp_path):
        data = {'rmsd_data': np.array([1.0, 2.0, 3.0]), 'time_unit': 'ns'}
        p1 = tmp_path / 'a.pkl'
        p2 = tmp_path / 'b.pkl'
        with open(p1, 'wb') as f:
            pickle.dump(data, f)
        with open(p2, 'wb') as f:
            pickle.dump(data, f)
        result = compare_pickle_pair(str(p1), str(p2), 'test.pkl', 1e-5, 1e-8)
        assert result.status == 'MATCH'

    def test_mismatching_pair(self, tmp_path):
        p1 = tmp_path / 'a.pkl'
        p2 = tmp_path / 'b.pkl'
        with open(p1, 'wb') as f:
            pickle.dump({'val': np.array([1.0])}, f)
        with open(p2, 'wb') as f:
            pickle.dump({'val': np.array([99.0])}, f)
        result = compare_pickle_pair(str(p1), str(p2), 'test.pkl', 1e-5, 1e-8)
        assert result.status == 'MISMATCH'

    def test_corrupt_pickle(self, tmp_path):
        p1 = tmp_path / 'a.pkl'
        p2 = tmp_path / 'b.pkl'
        p1.write_bytes(b'not a pickle')
        with open(p2, 'wb') as f:
            pickle.dump({'val': 1}, f)
        result = compare_pickle_pair(str(p1), str(p2), 'test.pkl', 1e-5, 1e-8)
        assert result.status == 'MISMATCH'
        assert 'Error' in result.details


# ─── compare_trees ───────────────────────────────────────────────────────────

class TestCompareTrees:

    def test_matching_trees(self, tmp_path):
        serial = tmp_path / 'serial'
        parallel = tmp_path / 'parallel'
        serial.mkdir()
        parallel.mkdir()

        data = {'rmsd': np.array([1.0, 2.0, 3.0])}
        for d in [serial, parallel]:
            with open(d / 'rmsd_plot.pkl', 'wb') as f:
                pickle.dump(data, f)

        results = compare_trees(str(serial), str(parallel))
        assert len(results) == 1
        assert results[0].status == 'MATCH'

    def test_file_only_in_serial(self, tmp_path):
        serial = tmp_path / 'serial'
        parallel = tmp_path / 'parallel'
        serial.mkdir()
        parallel.mkdir()

        with open(serial / 'rmsd_plot.pkl', 'wb') as f:
            pickle.dump({'data': 1}, f)

        results = compare_trees(str(serial), str(parallel))
        assert len(results) == 1
        assert results[0].status == 'ONLY_SERIAL'

    def test_file_only_in_parallel(self, tmp_path):
        serial = tmp_path / 'serial'
        parallel = tmp_path / 'parallel'
        serial.mkdir()
        parallel.mkdir()

        with open(parallel / 'extra.pkl', 'wb') as f:
            pickle.dump({'data': 1}, f)

        results = compare_trees(str(serial), str(parallel))
        assert len(results) == 1
        assert results[0].status == 'ONLY_PARALLEL'


# ─── format_report ───────────────────────────────────────────────────────────

class TestFormatReport:

    def test_all_match_report(self, tmp_path):
        results = [ComparisonResult('test.pkl', 'MATCH')]
        report = format_report(results, '/serial', '/parallel')
        assert 'numerically equivalent' in report
        assert 'Match' in report

    def test_mismatch_report(self):
        results = [ComparisonResult('test.pkl', 'MISMATCH', 'values differ')]
        report = format_report(results, '/s', '/p')
        assert 'Mismatch' in report
        assert 'values differ' in report

    def test_structural_mismatch_report(self):
        results = [ComparisonResult('a.pkl', 'ONLY_SERIAL')]
        report = format_report(results, '/s', '/p')
        assert 'only in serial' in report.lower()


# ─── _create_temp_ini ────────────────────────────────────────────────────────

class TestCreateTempIni:
    """Validate the temp-INI generator used by auto mode."""

    @pytest.fixture()
    def sample_ini(self, tmp_path):
        ini = tmp_path / 'test.ini'
        ini.write_text(
            '[systems]\n'
            'systems = ["A"]\n'
            'variations = {"A": ["wild"]}\n'
            '\n'
            '[paths]\n'
            'input_dir = /data\n'
            'output_dir = results\n'
            '\n'
            '[analysis]\n'
            'run_rmsd = true\n'
            'run_plots = true\n'
            '\n'
            '[parallelization]\n'
            'backend = serial\n'
            'n_workers = none\n'
        )
        return str(ini)

    def test_overrides_backend(self, sample_ini, tmp_path):
        import configparser
        result = _create_temp_ini(sample_ini, str(tmp_path / 'out'),
                                  backend='multiprocessing', n_workers=4)
        try:
            cp = configparser.RawConfigParser()
            cp.read(result)
            assert cp.get('parallelization', 'backend') == 'multiprocessing'
            assert cp.get('parallelization', 'n_workers') == '4'
        finally:
            os.unlink(result)

    def test_overrides_output_dir(self, sample_ini, tmp_path):
        import configparser
        out = str(tmp_path / 'custom_out')
        result = _create_temp_ini(sample_ini, out)
        try:
            cp = configparser.RawConfigParser()
            cp.read(result)
            assert cp.get('paths', 'output_dir') == os.path.abspath(out)
        finally:
            os.unlink(result)

    def test_forces_plots_off(self, sample_ini, tmp_path):
        import configparser
        result = _create_temp_ini(sample_ini, str(tmp_path / 'out'))
        try:
            cp = configparser.RawConfigParser()
            cp.read(result)
            assert cp.get('analysis', 'run_plots') == 'false'
        finally:
            os.unlink(result)

    def test_preserves_other_settings(self, sample_ini, tmp_path):
        import configparser
        result = _create_temp_ini(sample_ini, str(tmp_path / 'out'))
        try:
            cp = configparser.RawConfigParser()
            cp.read(result)
            assert cp.get('systems', 'systems') == '["A"]'
            assert cp.get('analysis', 'run_rmsd') == 'true'
        finally:
            os.unlink(result)

    def test_n_workers_none_writes_none(self, sample_ini, tmp_path):
        import configparser
        result = _create_temp_ini(sample_ini, str(tmp_path / 'out'),
                                  n_workers=None)
        try:
            cp = configparser.RawConfigParser()
            cp.read(result)
            assert cp.get('parallelization', 'n_workers') == 'none'
        finally:
            os.unlink(result)

    def test_creates_missing_parallelization_section(self, tmp_path):
        """INI without [parallelization] should still work."""
        import configparser
        ini = tmp_path / 'minimal.ini'
        ini.write_text('[systems]\nsystems = ["X"]\n')
        result = _create_temp_ini(str(ini), str(tmp_path / 'out'),
                                  backend='dask', n_workers=8)
        try:
            cp = configparser.RawConfigParser()
            cp.read(result)
            assert cp.get('parallelization', 'backend') == 'dask'
            assert cp.get('parallelization', 'n_workers') == '8'
        finally:
            os.unlink(result)


# ─── main() CLI validation ───────────────────────────────────────────────────

class TestMainCLI:
    """Test that main() validates argument combinations correctly."""

    def test_no_args_exits_with_error(self, monkeypatch):
        monkeypatch.setattr('sys.argv', ['compare_parallel_serial.py'])
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 2

    def test_config_with_dir_parallel_rejected(self, monkeypatch, tmp_path):
        ini = tmp_path / 'test.ini'
        ini.write_text('[systems]\nsystems = ["X"]\n')
        monkeypatch.setattr('sys.argv', [
            'compare_parallel_serial.py',
            '--config', str(ini),
            '--backend', 'multiprocessing',
            '--dir-parallel', '/tmp/x',
        ])
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 2

    def test_only_dir_serial_rejected(self, monkeypatch, tmp_path):
        monkeypatch.setattr('sys.argv', [
            'compare_parallel_serial.py',
            '--dir-serial', str(tmp_path),
        ])
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 2

    def test_config_missing_file_returns_2(self, monkeypatch):
        monkeypatch.setattr('sys.argv', [
            'compare_parallel_serial.py',
            '--config', '/nonexistent/path.ini',
            '--backend', 'multiprocessing',
        ])
        assert main() == 2

    def test_manual_mode_missing_dirs_returns_2(self, monkeypatch):
        monkeypatch.setattr('sys.argv', [
            'compare_parallel_serial.py',
            '--dir-serial', '/nonexistent/serial',
            '--dir-parallel', '/nonexistent/parallel',
        ])
        assert main() == 2

    def test_manual_mode_matching_pickles_returns_0(self, monkeypatch, tmp_path):
        serial = tmp_path / 'serial'
        parallel = tmp_path / 'parallel'
        serial.mkdir()
        parallel.mkdir()
        data = {'v': np.array([1.0, 2.0])}
        for d in (serial, parallel):
            with open(d / 'test.pkl', 'wb') as f:
                pickle.dump(data, f)
        monkeypatch.setattr('sys.argv', [
            'compare_parallel_serial.py',
            '--dir-serial', str(serial),
            '--dir-parallel', str(parallel),
        ])
        assert main() == 0

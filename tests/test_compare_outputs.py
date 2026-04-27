"""Tests for compare_outputs.py — semantic pickle comparison tool."""
import os
import pickle
import pytest
import numpy as np
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from compare_outputs import (
    _collect_pickle_tree,
    _compare_values,
    compare_pickle_pair,
    run_comparison,
    ComparisonResult,
    _render_markdown,
)


# ─── _collect_pickle_tree ────────────────────────────────────────────────────

class TestCollectPickleTree:
    def test_finds_pkl_files(self, tmp_path):
        (tmp_path / 'rmsd').mkdir()
        (tmp_path / 'rmsd' / 'a.pkl').write_bytes(b'x')
        (tmp_path / 'rmsd' / 'b.txt').write_bytes(b'x')   # ignored
        tree = _collect_pickle_tree(str(tmp_path))
        assert 'rmsd/a.pkl' in tree
        assert 'rmsd/b.txt' not in tree

    def test_skips_pipeline_info(self, tmp_path):
        (tmp_path / 'pipeline_info').mkdir()
        (tmp_path / 'pipeline_info' / 'c.pkl').write_bytes(b'x')
        tree = _collect_pickle_tree(str(tmp_path))
        assert len(tree) == 0

    def test_skips_rmsfit(self, tmp_path):
        (tmp_path / 'rmsfit_test.pkl').write_bytes(b'x')
        tree = _collect_pickle_tree(str(tmp_path))
        assert len(tree) == 0

    def test_empty_dir(self, tmp_path):
        tree = _collect_pickle_tree(str(tmp_path))
        assert tree == {}


# ─── _compare_values ─────────────────────────────────────────────────────────

class TestCompareValues:
    def test_identical_dicts(self):
        assert _compare_values({'a': 1}, {'a': 1}) == []

    def test_different_scalars(self):
        diffs = _compare_values(42, 99)
        assert len(diffs) == 1

    def test_numpy_allclose(self):
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([1.0, 2.0, 3.0 + 1e-9])
        assert _compare_values(a, b) == []

    def test_numpy_too_different(self):
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([1.0, 2.0, 4.0])
        diffs = _compare_values(a, b)
        assert len(diffs) == 1
        assert 'arrays differ' in diffs[0]

    def test_shape_mismatch(self):
        a = np.array([1, 2])
        b = np.array([1, 2, 3])
        diffs = _compare_values(a, b)
        assert 'shape mismatch' in diffs[0]

    def test_nested_dict(self):
        a = {'x': {'y': np.array([1.0])}}
        b = {'x': {'y': np.array([2.0])}}
        diffs = _compare_values(a, b)
        assert len(diffs) == 1
        assert 'x.y' in diffs[0]

    def test_type_mismatch(self):
        diffs = _compare_values(42, '42')
        assert 'type mismatch' in diffs[0]

    def test_list_comparison(self):
        assert _compare_values([1, 2], [1, 2]) == []
        diffs = _compare_values([1, 2], [1, 3])
        assert len(diffs) == 1

    def test_float_tolerance(self):
        assert _compare_values(1.0, 1.0 + 1e-10) == []
        diffs = _compare_values(1.0, 2.0)
        assert len(diffs) == 1

    def test_circular_reference_does_not_recurse_infinitely(self):
        """Objects with circular __dict__ references must not cause RecursionError."""
        class CircularObj:
            pass
        a = CircularObj()
        b = CircularObj()
        a.child = CircularObj()
        a.child.parent = a
        a.child.value = 42
        b.child = CircularObj()
        b.child.parent = b
        b.child.value = 42
        diffs = _compare_values(a, b)
        assert isinstance(diffs, list)

    def test_deeply_nested_objects_respect_depth_limit(self):
        """Deeply nested structures beyond _MAX_COMPARE_DEPTH are treated as equal."""
        from compare_outputs import _MAX_COMPARE_DEPTH
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
        a_cur.value = 1
        b_cur.value = 999
        diffs = _compare_values(a_root, b_root)
        assert all('value' not in d for d in diffs)

    def test_dict_circular_reference(self):
        """Dicts with circular references must not cause RecursionError."""
        a = {'key': 'val'}
        b = {'key': 'val'}
        a['self'] = a
        b['self'] = b
        diffs = _compare_values(a, b)
        assert isinstance(diffs, list)


# ─── compare_pickle_pair ─────────────────────────────────────────────────────

class TestComparePicklePair:
    def test_matching_pickles(self, tmp_path):
        data = {'values': np.array([1.0, 2.0, 3.0])}
        a = tmp_path / 'a.pkl'
        b = tmp_path / 'b.pkl'
        a.write_bytes(pickle.dumps(data))
        b.write_bytes(pickle.dumps(data))
        result = compare_pickle_pair(str(a), str(b), 'test.pkl')
        assert result.ok
        assert result.status == 'MATCH'

    def test_mismatching_pickles(self, tmp_path):
        a = tmp_path / 'a.pkl'
        b = tmp_path / 'b.pkl'
        a.write_bytes(pickle.dumps({'v': np.array([1.0])}))
        b.write_bytes(pickle.dumps({'v': np.array([99.0])}))
        result = compare_pickle_pair(str(a), str(b), 'test.pkl')
        assert not result.ok
        assert result.status == 'MISMATCH'

    def test_corrupt_pickle(self, tmp_path):
        a = tmp_path / 'a.pkl'
        b = tmp_path / 'b.pkl'
        a.write_bytes(pickle.dumps({'ok': True}))
        b.write_bytes(b'not a pickle')
        result = compare_pickle_pair(str(a), str(b), 'test.pkl')
        assert result.status == 'MISMATCH'
        assert 'Failed to load' in result.details


# ─── run_comparison ──────────────────────────────────────────────────────────

class TestRunComparison:
    def _setup_dirs(self, tmp_path, py_files, nf_files):
        py_dir = tmp_path / 'python'
        nf_dir = tmp_path / 'nextflow'
        for name, data in py_files.items():
            p = py_dir / name
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(pickle.dumps(data))
        for name, data in nf_files.items():
            p = nf_dir / name
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(pickle.dumps(data))
        return str(py_dir), str(nf_dir)

    def test_all_match(self, tmp_path):
        data = {'v': np.array([1.0, 2.0])}
        py, nf = self._setup_dirs(tmp_path,
            {'rmsd/test.pkl': data}, {'rmsd/test.pkl': data})
        results = run_comparison(py, nf)
        assert all(r.ok for r in results)

    def test_only_python(self, tmp_path):
        py, nf = self._setup_dirs(tmp_path,
            {'rmsd/test.pkl': {'x': 1}}, {})
        results = run_comparison(py, nf)
        assert results[0].status == 'ONLY_PYTHON'

    def test_only_nextflow(self, tmp_path):
        py, nf = self._setup_dirs(tmp_path,
            {}, {'rmsd/test.pkl': {'x': 1}})
        results = run_comparison(py, nf)
        assert results[0].status == 'ONLY_NEXTFLOW'


# ─── Markdown report ────────────────────────────────────────────────────────

class TestRenderMarkdown:
    def test_contains_header(self):
        md = _render_markdown([], '/a', '/b')
        assert '# Output Comparison Report' in md

    def test_contains_summary(self):
        results = [ComparisonResult('f.pkl', 'MATCH', 'ok')]
        md = _render_markdown(results, '/a', '/b')
        assert '1 match' in md
        assert '0 mismatch' in md

    def test_empty_results(self):
        md = _render_markdown([], '/a', '/b')
        assert 'No pickle files' in md

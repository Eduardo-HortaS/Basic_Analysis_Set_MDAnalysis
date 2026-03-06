"""Additional tests to improve code coverage across the pipeline.

Covers:
- ini_to_nf_params.main() CLI
- compare_parallel_serial.main() CLI
- compare_outputs.main() CLI
- executor: validate_inputs, print_config_summary, dry-run paths,
  _parse_json_or_path (file branch), _parse_optional_int
- run_rog_analysis.main() / run_hbonds_analysis.main() argument parsing
"""
import os
import pickle
import sys
import textwrap

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from executor import (
    _parse_json_or_path,
    _parse_optional_int,
    validate_inputs,
    print_config_summary,
    run_rmsd,
    run_rmsf,
    run_2d_rmsd,
    run_rog,
    run_hbonds,
)
from ini_to_nf_params import _dump_yaml, main as nf_main
from compare_outputs import main as co_main, _render_markdown, ComparisonResult
from compare_parallel_serial import main as cps_main, format_report


# ─── Helper: write a minimal valid INI ────────────────────────────────────────

def _write_ini(tmp_path, extra_sections='', *, filename='test.ini'):
    """Write a minimal valid INI config and return its path."""
    ini = tmp_path / filename
    ini.write_text(textwrap.dedent(f"""\
        [systems]
        systems = ["SysA"]
        variations = {{"SysA": ["wild"]}}
        num_replicates = 1
        traj_format = dcd
        top_format = top
        start_frame = 0

        [paths]
        input_dir = {tmp_path}
        output_dir = {tmp_path / 'results'}

        [analysis]
        run_rmsd = true
        run_rmsf = false
        run_rog = false
        run_hbonds = false
        run_plots = false

        [rmsd]
        time_interval_between_frames = 2.0

        {extra_sections}
    """))
    return str(ini)


def _make_cfg(tmp_path, **overrides):
    """Return a cfg dict suitable for executor helpers."""
    cfg = {
        'systems': ['SysA'],
        'variations': {'SysA': ['wild']},
        'num_replicates': 1,
        'traj_format': 'dcd',
        'top_format': 'top',
        'start_frame': 0,
        'input_dir': str(tmp_path),
        'output_dir': str(tmp_path / 'results'),
        'plot_dir': str(tmp_path / 'results' / 'plots'),
        'run_rmsd': True, 'run_rmsf': False, 'run_2d_rmsd': False,
        'run_rog': False, 'run_hbonds': False, 'run_plots': False,
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
        'rmsd_show_kde': True,
        'rog_show_kde': True,
        'rmsf_plot_type': 'line',
        'twod_rmsd_cmap': 'viridis',
        'parallel_backend': 'serial',
        'n_workers': None,
        'plot_groups': {},
        'replicate_mode': 'separate',
    }
    cfg.update(overrides)
    return cfg


# ═══════════════════════════════════════════════════════════════════════════════
#  executor.py — parsing helpers
# ═══════════════════════════════════════════════════════════════════════════════

class TestParseJsonOrPathFromFile:
    """Cover the file-reading branch of _parse_json_or_path."""

    def test_reads_json_from_file(self, tmp_path):
        f = tmp_path / 'data.json'
        f.write_text('["SysA", "SysB"]')
        result = _parse_json_or_path(str(f))
        assert result == ['SysA', 'SysB']

    def test_reads_dict_from_file(self, tmp_path):
        f = tmp_path / 'data.json'
        f.write_text('{"SysA": ["wild"]}')
        result = _parse_json_or_path(str(f))
        assert result == {'SysA': ['wild']}


class TestParseOptionalInt:
    """Cover _parse_optional_int."""

    def test_returns_int(self):
        assert _parse_optional_int('  4  ') == 4

    def test_returns_none_for_none_string(self):
        assert _parse_optional_int('none') is None

    def test_returns_none_for_empty(self):
        assert _parse_optional_int('') is None


# ═══════════════════════════════════════════════════════════════════════════════
#  executor.py — validate_inputs
# ═══════════════════════════════════════════════════════════════════════════════

class TestValidateInputs:
    def test_all_files_missing(self, tmp_path):
        cfg = _make_cfg(tmp_path)
        found, missing = validate_inputs(cfg)
        assert found == 0
        assert len(missing) > 0

    def test_all_files_present(self, tmp_path):
        cfg = _make_cfg(tmp_path)
        # Create expected files
        data_dir = tmp_path / 'SysA' / 'wild'
        data_dir.mkdir(parents=True)
        (data_dir / 'SysA_system_wild.top').write_text('x')
        (data_dir / 'SysA_production_wild_rep_1.dcd').write_text('x')
        found, missing = validate_inputs(cfg)
        assert found == 2
        assert missing == []

    def test_system_not_in_variations(self, tmp_path):
        cfg = _make_cfg(tmp_path, systems=['SysA', 'SysB'])
        # SysB not in variations
        _found, missing = validate_inputs(cfg)
        assert any('SysB' in m for m in missing)


# ═══════════════════════════════════════════════════════════════════════════════
#  executor.py — print_config_summary
# ═══════════════════════════════════════════════════════════════════════════════

class TestPrintConfigSummary:
    def test_basic_output(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path)
        print_config_summary(cfg)
        captured = capsys.readouterr().out
        assert 'CONFIGURATION SUMMARY' in captured
        assert 'SysA' in captured

    def test_parallel_mode_displayed(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path, parallel_backend='multiprocessing', n_workers=4)
        print_config_summary(cfg)
        captured = capsys.readouterr().out
        assert 'multiprocessing' in captured
        assert '4' in captured

    def test_serial_mode_displayed(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path, parallel_backend='serial')
        print_config_summary(cfg)
        captured = capsys.readouterr().out
        assert 'serial' in captured

    def test_all_analyses_listed(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path,
                        run_rmsd=True, run_rmsf=True, run_2d_rmsd=True,
                        run_rog=True, run_hbonds=True, run_plots=True)
        print_config_summary(cfg)
        captured = capsys.readouterr().out
        for name in ['RMSD', 'RMSF', '2D-RMSD', 'RoG', 'H-bonds']:
            assert name in captured

    def test_no_analyses_shows_none(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path,
                        run_rmsd=False, run_rmsf=False, run_2d_rmsd=False,
                        run_rog=False, run_hbonds=False)
        print_config_summary(cfg)
        captured = capsys.readouterr().out
        assert '(none)' in captured

    def test_between_pairs_hbonds(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path, run_hbonds=True,
                        between_pairs=[['protein', 'resname URA']])
        print_config_summary(cfg)
        captured = capsys.readouterr().out
        assert 'pair-focused' in captured

    def test_acceptors_sel_hbonds(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path, run_hbonds=True,
                        between_pairs=None, acceptors_sel='resname URA')
        print_config_summary(cfg)
        captured = capsys.readouterr().out
        assert 'atom-focused' in captured

    def test_chain_intervals_displayed(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path, run_rmsf=True,
                        chain_intervals={'SysA': {'A': [1, 100]}})
        print_config_summary(cfg)
        captured = capsys.readouterr().out
        assert 'chains' in captured.lower()

    def test_multiple_ref_selections(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path, run_rmsd=True,
                        ref_selection=['sel_a', 'sel_b'])
        print_config_summary(cfg)
        captured = capsys.readouterr().out
        assert 'sel_a' in captured

    def test_group_selections(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path, run_rmsd=True,
                        group_selections=['protein', 'resname URA'])
        print_config_summary(cfg)
        captured = capsys.readouterr().out
        assert 'sels' in captured.lower()

    def test_plot_groups(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path,
                        plot_groups={'grp1': [('SysA', 'wild')]},
                        replicate_mode='average')
        print_config_summary(cfg)
        captured = capsys.readouterr().out
        assert 'grp1' in captured


# ═══════════════════════════════════════════════════════════════════════════════
#  executor.py — dry_run paths
# ═══════════════════════════════════════════════════════════════════════════════

class TestDryRunPaths:
    """Cover the dry_run=True branch of each analysis runner."""

    def test_run_rmsd_dry_run(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path)
        cfg['_work_dir'] = str(tmp_path)
        run_rmsd(cfg, dry_run=True)
        out = capsys.readouterr().out
        assert 'DRY RUN' in out

    def test_run_rmsf_dry_run(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path)
        cfg['_work_dir'] = str(tmp_path)
        run_rmsf(cfg, dry_run=True)
        out = capsys.readouterr().out
        assert 'DRY RUN' in out

    def test_run_rmsf_dry_run_with_chains(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path, chain_intervals={'SysA': {'A': [1, 100]}})
        cfg['_work_dir'] = str(tmp_path)
        run_rmsf(cfg, dry_run=True)
        out = capsys.readouterr().out
        assert 'Chain intervals' in out

    def test_run_2d_rmsd_dry_run(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path)
        cfg['_work_dir'] = str(tmp_path)
        run_2d_rmsd(cfg, dry_run=True)
        out = capsys.readouterr().out
        assert 'DRY RUN' in out

    def test_run_rog_dry_run(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path)
        cfg['_work_dir'] = str(tmp_path)
        run_rog(cfg, dry_run=True)
        out = capsys.readouterr().out
        assert 'DRY RUN' in out

    def test_run_hbonds_dry_run(self, tmp_path, capsys):
        cfg = _make_cfg(tmp_path)
        cfg['_work_dir'] = str(tmp_path)
        run_hbonds(cfg, dry_run=True)
        out = capsys.readouterr().out
        assert 'DRY RUN' in out


# ═══════════════════════════════════════════════════════════════════════════════
#  ini_to_nf_params.main() — CLI
# ═══════════════════════════════════════════════════════════════════════════════

class TestIniToNfParamsCLI:
    def test_main_stdout(self, tmp_path, monkeypatch, capsys):
        """main() prints YAML to stdout when no -o is given."""
        ini = _write_ini(tmp_path)
        monkeypatch.setattr('sys.argv', ['ini_to_nf_params', ini])
        nf_main()
        out = capsys.readouterr().out
        assert 'num_replicates' in out

    def test_main_output_file(self, tmp_path, monkeypatch):
        """main() writes to file when -o is given."""
        ini = _write_ini(tmp_path)
        out_file = str(tmp_path / 'params.yml')
        monkeypatch.setattr('sys.argv', ['ini_to_nf_params', ini, '-o', out_file])
        nf_main()
        assert os.path.exists(out_file)
        with open(out_file, encoding='utf-8') as f:
            content = f.read()
        assert 'num_replicates' in content

    def test_dump_yaml_simple_string_unquoted(self):
        """Simple strings without special chars are unquoted."""
        yaml = _dump_yaml({'format': 'dcd'})
        assert 'format: dcd' in yaml

    def test_dump_yaml_fallback_type(self):
        """Non-standard types get quoted as strings."""
        yaml = _dump_yaml({'data': (1, 2, 3)})
        assert 'data:' in yaml


# ═══════════════════════════════════════════════════════════════════════════════
#  compare_outputs.main() — CLI
# ═══════════════════════════════════════════════════════════════════════════════

class TestCompareOutputsCLI:
    def test_main_matching(self, tmp_path, monkeypatch):
        """main() prints report and returns normally when all pickles match."""
        py_dir = tmp_path / 'py'
        nf_dir = tmp_path / 'nf'
        for d in (py_dir / 'rmsd', nf_dir / 'rmsd'):
            d.mkdir(parents=True)
        data = {'v': np.array([1.0, 2.0])}
        for d in (py_dir / 'rmsd', nf_dir / 'rmsd'):
            with open(d / 'test.pkl', 'wb') as f:
                pickle.dump(data, f)
        monkeypatch.setattr('sys.argv', [
            'compare_outputs',
            '--dir-python', str(py_dir),
            '--dir-nextflow', str(nf_dir),
        ])
        # No SystemExit when everything matches
        co_main()

    def test_main_mismatch(self, tmp_path, monkeypatch):
        """main() exits 1 when pickles differ."""
        py_dir = tmp_path / 'py' / 'rmsd'
        nf_dir = tmp_path / 'nf' / 'rmsd'
        py_dir.mkdir(parents=True)
        nf_dir.mkdir(parents=True)
        with open(py_dir / 'test.pkl', 'wb') as f:
            pickle.dump({'v': np.array([1.0])}, f)
        with open(nf_dir / 'test.pkl', 'wb') as f:
            pickle.dump({'v': np.array([999.0])}, f)
        monkeypatch.setattr('sys.argv', [
            'compare_outputs',
            '--dir-python', str(tmp_path / 'py'),
            '--dir-nextflow', str(tmp_path / 'nf'),
        ])
        with pytest.raises(SystemExit) as exc_info:
            co_main()
        assert exc_info.value.code == 1

    def test_main_output_file(self, tmp_path, monkeypatch):
        """main() writes report to file with -o flag."""
        py_dir = tmp_path / 'py' / 'rmsd'
        nf_dir = tmp_path / 'nf' / 'rmsd'
        py_dir.mkdir(parents=True)
        nf_dir.mkdir(parents=True)
        data = {'v': np.array([1.0])}
        for d in (py_dir, nf_dir):
            with open(d / 'test.pkl', 'wb') as f:
                pickle.dump(data, f)
        report_file = str(tmp_path / 'report.md')
        monkeypatch.setattr('sys.argv', [
            'compare_outputs',
            '--dir-python', str(tmp_path / 'py'),
            '--dir-nextflow', str(tmp_path / 'nf'),
            '-o', report_file,
        ])
        co_main()  # No SystemExit on match
        assert os.path.exists(report_file)

    def test_main_only_one_side(self, tmp_path, monkeypatch):
        """main() exits 2 when structural mismatch (files only in one dir)."""
        py_dir = tmp_path / 'py' / 'rmsd'
        nf_dir = tmp_path / 'nf'
        py_dir.mkdir(parents=True)
        nf_dir.mkdir(parents=True)
        with open(py_dir / 'test.pkl', 'wb') as f:
            pickle.dump({'v': 1}, f)
        monkeypatch.setattr('sys.argv', [
            'compare_outputs',
            '--dir-python', str(tmp_path / 'py'),
            '--dir-nextflow', str(nf_dir),
        ])
        with pytest.raises(SystemExit) as exc_info:
            co_main()
        assert exc_info.value.code == 2


# ═══════════════════════════════════════════════════════════════════════════════
#  compare_parallel_serial.main() — CLI
# ═══════════════════════════════════════════════════════════════════════════════

class TestCompareParallelSerialCLI:
    def test_main_missing_serial_dir(self, tmp_path, monkeypatch):
        """main() returns 2 when serial dir doesn't exist."""
        monkeypatch.setattr('sys.argv', [
            'compare_parallel_serial',
            '--dir-serial', str(tmp_path / 'nonexistent_serial'),
            '--dir-parallel', str(tmp_path),
        ])
        rc = cps_main()
        assert rc == 2

    def test_main_missing_parallel_dir(self, tmp_path, monkeypatch):
        """main() returns 2 when parallel dir doesn't exist."""
        monkeypatch.setattr('sys.argv', [
            'compare_parallel_serial',
            '--dir-serial', str(tmp_path),
            '--dir-parallel', str(tmp_path / 'nonexistent_parallel'),
        ])
        rc = cps_main()
        assert rc == 2

    def test_main_no_pickles(self, tmp_path, monkeypatch):
        """main() returns 2 when no pkl files exist."""
        serial = tmp_path / 'serial'
        parallel = tmp_path / 'parallel'
        serial.mkdir()
        parallel.mkdir()
        monkeypatch.setattr('sys.argv', [
            'compare_parallel_serial',
            '--dir-serial', str(serial),
            '--dir-parallel', str(parallel),
        ])
        rc = cps_main()
        assert rc == 2

    def test_main_matching_pickles(self, tmp_path, monkeypatch):
        """main() returns 0 when pickles match."""
        serial = tmp_path / 'serial' / 'rmsd'
        parallel = tmp_path / 'parallel' / 'rmsd'
        serial.mkdir(parents=True)
        parallel.mkdir(parents=True)
        data = {'v': np.array([1.0, 2.0, 3.0])}
        for d in (serial, parallel):
            with open(d / 'test.pkl', 'wb') as f:
                pickle.dump(data, f)
        monkeypatch.setattr('sys.argv', [
            'compare_parallel_serial',
            '--dir-serial', str(tmp_path / 'serial'),
            '--dir-parallel', str(tmp_path / 'parallel'),
        ])
        rc = cps_main()
        assert rc == 0

    def test_main_mismatching_pickles(self, tmp_path, monkeypatch):
        """main() returns 1 when pickles mismatch."""
        serial = tmp_path / 'serial' / 'rmsd'
        parallel = tmp_path / 'parallel' / 'rmsd'
        serial.mkdir(parents=True)
        parallel.mkdir(parents=True)
        with open(serial / 'test.pkl', 'wb') as f:
            pickle.dump({'v': np.array([1.0])}, f)
        with open(parallel / 'test.pkl', 'wb') as f:
            pickle.dump({'v': np.array([999.0])}, f)
        monkeypatch.setattr('sys.argv', [
            'compare_parallel_serial',
            '--dir-serial', str(tmp_path / 'serial'),
            '--dir-parallel', str(tmp_path / 'parallel'),
        ])
        rc = cps_main()
        assert rc == 1

    def test_main_output_file(self, tmp_path, monkeypatch):
        """main() writes report when -o is given."""
        serial = tmp_path / 'serial' / 'rmsd'
        parallel = tmp_path / 'parallel' / 'rmsd'
        serial.mkdir(parents=True)
        parallel.mkdir(parents=True)
        data = {'v': np.array([1.0])}
        for d in (serial, parallel):
            with open(d / 'test.pkl', 'wb') as f:
                pickle.dump(data, f)
        report_path = str(tmp_path / 'report.md')
        monkeypatch.setattr('sys.argv', [
            'compare_parallel_serial',
            '--dir-serial', str(tmp_path / 'serial'),
            '--dir-parallel', str(tmp_path / 'parallel'),
            '-o', report_path,
        ])
        rc = cps_main()
        assert rc == 0
        assert os.path.exists(report_path)

    def test_main_structural_mismatch(self, tmp_path, monkeypatch):
        """main() returns 2 when files exist only on one side."""
        serial = tmp_path / 'serial' / 'rmsd'
        parallel = tmp_path / 'parallel'
        serial.mkdir(parents=True)
        parallel.mkdir(parents=True)
        with open(serial / 'test.pkl', 'wb') as f:
            pickle.dump({'v': 1}, f)
        monkeypatch.setattr('sys.argv', [
            'compare_parallel_serial',
            '--dir-serial', str(tmp_path / 'serial'),
            '--dir-parallel', str(parallel),
        ])
        rc = cps_main()
        assert rc == 2


# ═══════════════════════════════════════════════════════════════════════════════
#  compare_parallel_serial — format_report edge cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestFormatReportEdgeCases:
    def test_all_pass_message(self):
        """When all match, report should show the success message."""
        from compare_parallel_serial import ComparisonResult as CPSResult
        results = [CPSResult('test.pkl', 'MATCH', '', None, None)]
        report = format_report(results, '/serial', '/parallel')
        assert 'numerically equivalent' in report

    def test_mismatch_section(self):
        """Mismatches should appear in the report."""
        from compare_parallel_serial import ComparisonResult as CPSResult
        results = [CPSResult('test.pkl', 'MISMATCH', 'values differ',
                             1.5e-3, 2.3e-4)]
        report = format_report(results, '/serial', '/parallel')
        assert 'Mismatches' in report
        assert 'test.pkl' in report


# ═══════════════════════════════════════════════════════════════════════════════
#  compare_outputs — additional _compare_values coverage
# ═══════════════════════════════════════════════════════════════════════════════

class TestCompareValuesAdditional:
    """Cover edge-case branches in _compare_values."""

    def test_list_length_mismatch(self):
        from compare_outputs import _compare_values
        diffs = _compare_values([1, 2], [1, 2, 3])
        assert any('length mismatch' in d for d in diffs)

    def test_float_mismatch(self):
        from compare_outputs import _compare_values
        diffs = _compare_values(1.0, 999.0)
        assert len(diffs) == 1

    def test_scalar_mismatch(self):
        from compare_outputs import _compare_values
        diffs = _compare_values('hello', 'world')
        assert len(diffs) == 1

    def test_dict_extra_key(self):
        from compare_outputs import _compare_values
        diffs = _compare_values({'a': 1}, {'a': 1, 'b': 2})
        assert any('missing in left' in d for d in diffs)

    def test_dict_missing_key(self):
        from compare_outputs import _compare_values
        diffs = _compare_values({'a': 1, 'b': 2}, {'a': 1})
        assert any('missing in right' in d for d in diffs)

    def test_object_with_results_attr(self):
        """Objects with .results attribute compare via results."""
        from compare_outputs import _compare_values

        class FakeResults:
            pass

        class FakeObj:
            def __init__(self, val):
                self.results = FakeResults()
                self.results.data = val

        a = FakeObj(np.array([1.0, 2.0]))
        b = FakeObj(np.array([1.0, 2.0]))
        diffs = _compare_values(a, b)
        assert diffs == []

    def test_generic_object_comparison(self):
        """Objects without .results compare via __dict__."""
        from compare_outputs import _compare_values

        class Simple:
            def __init__(self, x):
                self.x = x

        diffs = _compare_values(Simple(1), Simple(1))
        assert diffs == []

    def test_array_shape_mismatch(self):
        from compare_outputs import _compare_values
        diffs = _compare_values(np.array([1, 2]), np.array([1, 2, 3]))
        assert any('shape mismatch' in d for d in diffs)


# ═══════════════════════════════════════════════════════════════════════════════
#  compare_outputs — _render_markdown with mismatch details
# ═══════════════════════════════════════════════════════════════════════════════

class TestRenderMarkdownDetails:
    def test_mismatch_row(self):
        results = [ComparisonResult('f.pkl', 'MISMATCH', 'values differ')]
        md = _render_markdown(results, '/a', '/b')
        assert 'MISMATCH' in md
        assert 'f.pkl' in md

    def test_only_python_row(self):
        results = [ComparisonResult('f.pkl', 'ONLY_PYTHON', 'Not found')]
        md = _render_markdown(results, '/a', '/b')
        assert 'only-python' in md

    def test_pipe_in_details_escaped(self):
        results = [ComparisonResult('f.pkl', 'MATCH', 'a|b')]
        md = _render_markdown(results, '/a', '/b')
        assert 'a\\|b' in md


# ═══════════════════════════════════════════════════════════════════════════════
#  run_rog_analysis — main() argument parsing
# ═══════════════════════════════════════════════════════════════════════════════

class TestRoGMainArgs:
    def test_invalid_json_exits(self, monkeypatch):
        """main() returns 1 on invalid JSON systems."""
        from run_rog_analysis import main as rog_main
        monkeypatch.setattr('sys.argv', [
            'run_rog_analysis',
            '--systems', '{invalid json',
            '--variations', '{}',
        ])
        rc = rog_main()
        assert rc == 1

    def test_json_file_path(self, tmp_path, monkeypatch):
        """main() reads systems/variations from JSON files."""
        from run_rog_analysis import main as rog_main
        sys_file = tmp_path / 'systems.json'
        var_file = tmp_path / 'variations.json'
        sys_file.write_text('["SysA"]')
        var_file.write_text('{"SysA": ["wild"]}')
        # This will fail at file loading stage since data doesn't exist,
        # but parse should succeed
        monkeypatch.setattr('sys.argv', [
            'run_rog_analysis',
            '--systems', str(sys_file),
            '--variations', str(var_file),
            '--num-replicates', '1',
            '--start-frame', '0',
        ])
        # Will raise because trajectory files don't exist
        with pytest.raises((FileNotFoundError, OSError)):
            rog_main()


# ═══════════════════════════════════════════════════════════════════════════════
#  run_hbonds_analysis — main() argument parsing
# ═══════════════════════════════════════════════════════════════════════════════

class TestHBondsMainArgs:
    def test_invalid_json_exits(self, monkeypatch):
        """main() returns 1 on invalid JSON."""
        from run_hbonds_analysis import main as hbonds_main
        monkeypatch.setattr('sys.argv', [
            'run_hbonds_analysis',
            '--systems', '{bad json',
            '--variations', '{}',
            '--between-pairs', '[["protein", "resname URA"]]',
        ])
        rc = hbonds_main()
        assert rc == 1


# ═══════════════════════════════════════════════════════════════════════════════
#  run_rms_analysis — main() argument parsing
# ═══════════════════════════════════════════════════════════════════════════════

class TestRMSMainArgs:
    def test_invalid_json_exits(self, monkeypatch):
        """main() returns 1 on invalid JSON."""
        from run_rms_analysis import main as rms_main
        monkeypatch.setattr('sys.argv', [
            'run_rms_analysis',
            '--systems', '{bad json',
            '--variations', '{}',
            '--analysis', 'RMSD',
            '--target-selection', 'protein and backbone',
            '--ref-selection', 'protein and backbone',
        ])
        rc = rms_main()
        assert rc == 1

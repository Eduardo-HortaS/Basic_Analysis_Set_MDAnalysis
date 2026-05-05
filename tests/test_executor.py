"""
Tests for executor.py — INI parsing, normalization helpers,
and decoupled analysis/plot toggle validation.

Covers:
  - _normalize_chain_intervals format detection and conversion
  - Per-analysis plot toggle parsing via load_config
  - validate_plot_pickles: error when plot_X=true but run_X=false and no pickles
  - validate_plot_pickles: passes when pickles exist
  - validate_plot_pickles: passes when run_X=true (analysis will produce pickles)
  - validate_plot_pickles: skipped analyses with plot=false produce no errors
"""
import os
import sys
import pickle
import tempfile
import pytest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from executor import (
    _normalize_chain_intervals,
    load_config,
    validate_plot_pickles,
    _collect_pickles,
    setup_workdir,
    _clean_ephemeral_files,
    _ANALYSIS_SUBDIRS,
)


SYSTEMS = ['Ung_G-C_4', 'Mug_G-C_4']


# ─── Helper to write a minimal INI and load it ───────────────────────────────

def _write_ini(tmp_path, analysis_lines, rmsd_lines='time_interval_between_frames = 2.0'):
    """Write a minimal valid INI file and return its path.

    Parameters
    ----------
    rmsd_lines : str
        Extra lines for the [rmsd] section.  Defaults to a valid
        time_interval_between_frames so that the mandatory check passes.
    """
    ini = tmp_path / "test.ini"
    content = (
        "[systems]\n"
        "systems = [\"SysA\"]\n"
        "variations = {\"SysA\": [\"v1\"]}\n"
        "num_replicates = 1\n"
        "\n"
        "[paths]\n"
        f"input_dir = {tmp_path}\n"
        f"output_dir = {tmp_path / 'results'}\n"
        "\n"
        "[analysis]\n"
        f"{analysis_lines}\n"
        "\n"
        "[rmsd]\n"
        f"{rmsd_lines}\n"
    )
    ini.write_text(content)
    return str(ini)


def _make_base_cfg(**overrides):
    """Return a minimal cfg dict suitable for validate_plot_pickles."""
    cfg = {
        'systems': ['SysA'],
        'variations': {'SysA': ['v1']},
        'num_replicates': 1,
        # Defaults: all analyses and plots off
        'run_rmsd': False,    'plot_rmsd': False,
        'run_rmsf': False,    'plot_rmsf': False,
        'run_2d_rmsd': False, 'plot_2d_rmsd': False,
        'run_rog': False,     'plot_rog': False,
        'run_hbonds': False,  'plot_hbonds': False,
        'run_plots': True,
    }
    cfg.update(overrides)
    return cfg


# ─── _normalize_chain_intervals ──────────────────────────────────────────────

class TestNormalizeChainIntervals:
    """Tests for _normalize_chain_intervals."""

    def test_none_returns_none(self):
        assert _normalize_chain_intervals(None, SYSTEMS) is None

    def test_per_system_format_passes_through(self):
        raw = {
            "Ung_G-C_4": {"A": [1, 239]},
            "Mug_G-C_4": {"A": [1, 211]},
        }
        result = _normalize_chain_intervals(raw, SYSTEMS)
        assert result == raw

    def test_legacy_flat_format_broadcasts(self):
        raw = {"A": [1, 120], "B": [121, 239]}
        result = _normalize_chain_intervals(raw, SYSTEMS)
        assert set(result.keys()) == set(SYSTEMS)
        for sys_name in SYSTEMS:
            assert result[sys_name] == {"A": [1, 120], "B": [121, 239]}

    def test_legacy_broadcast_creates_independent_copies(self):
        """Modifying one system's intervals must not affect others."""
        raw = {"A": [1, 100]}
        result = _normalize_chain_intervals(raw, SYSTEMS)
        result['Ung_G-C_4']['A'][1] = 999
        assert result['Mug_G-C_4']['A'] == [1, 100]

    def test_mixed_format_raises(self):
        raw = {"Ung_G-C_4": {"A": [1, 239]}, "B": [1, 100]}
        with pytest.raises(ValueError, match="mixed format"):
            _normalize_chain_intervals(raw, SYSTEMS)

    def test_empty_dict_raises(self):
        with pytest.raises(ValueError, match="JSON dict"):
            _normalize_chain_intervals([], SYSTEMS)

    def test_non_dict_raises(self):
        with pytest.raises(ValueError, match="JSON dict"):
            _normalize_chain_intervals("not a dict", SYSTEMS)

    def test_unrecognized_value_type_raises(self):
        raw = {"A": 42}
        with pytest.raises(ValueError, match="unrecognized format"):
            _normalize_chain_intervals(raw, SYSTEMS)


# ─── Per-analysis plot toggle parsing ─────────────────────────────────────────

class TestPlotToggleParsing:
    """Tests that load_config correctly reads plot_* toggles from INI."""

    def test_plot_toggles_default_to_true(self, tmp_path):
        """When plot_* keys are absent, they should default to true."""
        ini = _write_ini(tmp_path, "run_rmsd = true\nrun_rmsf = false")
        cfg = load_config(ini)
        assert cfg['plot_rmsd'] is True
        assert cfg['plot_rmsf'] is True
        assert cfg['plot_2d_rmsd'] is True
        assert cfg['plot_rog'] is True
        assert cfg['plot_hbonds'] is True

    def test_plot_toggles_explicit_false(self, tmp_path):
        """Explicitly setting plot_rmsd = false should be honoured."""
        ini = _write_ini(tmp_path,
                         "run_rmsd = true\nplot_rmsd = false\n"
                         "run_rmsf = true\nplot_rmsf = false")
        cfg = load_config(ini)
        assert cfg['plot_rmsd'] is False
        assert cfg['plot_rmsf'] is False

    def test_analysis_off_plot_on(self, tmp_path):
        """run_X=false with plot_X=true should both be parsed correctly."""
        ini = _write_ini(tmp_path,
                         "run_rog = false\nplot_rog = true")
        cfg = load_config(ini)
        assert cfg['run_rog'] is False
        assert cfg['plot_rog'] is True

    def test_rmsf_group_selections_none_with_inline_comment(self, tmp_path):
        """Inline comments after 'none' must not trigger JSON decode errors."""
        ini = tmp_path / "test_inline_comment.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "\n"
            "[analysis]\n"
            "run_rmsd = false\n"
            "run_rmsf = true\n"
            "\n"
            "[rmsd]\n"
            "time_interval_between_frames = none\n"
            "\n"
            "[rmsf]\n"
            "chain_intervals = none\n"
            "group_selections = none ; e.g., [\"chainid A\", \"chainid B\"]\n"
        )
        ini.write_text(content)

        cfg = load_config(str(ini))
        assert cfg['rmsf_group_selections'] is None

    def test_hbonds_preset_parsed(self, tmp_path):
        """hbonds_preset should be read from [hbonds] and normalized."""
        ini = tmp_path / "test_hbonds_preset.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "\n"
            "[analysis]\n"
            "run_rmsd = false\n"
            "run_hbonds = true\n"
            "\n"
            "[rmsd]\n"
            "time_interval_between_frames = none\n"
            "\n"
            "[hbonds]\n"
            "hbonds_preset = protein_nucleic\n"
            "between_pairs = none\n"
        )
        ini.write_text(content)

        cfg = load_config(str(ini))
        assert cfg['hbonds_preset'] == 'protein_nucleic'

    def test_hbonds_preset_none_falls_back_to_custom(self, tmp_path):
        """hbonds_preset = none should normalize to the default custom preset."""
        ini = tmp_path / "test_hbonds_preset_none.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "\n"
            "[analysis]\n"
            "run_rmsd = false\n"
            "run_hbonds = true\n"
            "\n"
            "[rmsd]\n"
            "time_interval_between_frames = none\n"
            "\n"
            "[hbonds]\n"
            "hbonds_preset = none\n"
            "between_pairs = none\n"
        )
        ini.write_text(content)

        cfg = load_config(str(ini))
        assert cfg['hbonds_preset'] == 'custom'


# ─── validate_plot_pickles ────────────────────────────────────────────────────

class TestValidatePlotPickles:
    """Tests for validate_plot_pickles — the pickle pre-flight check."""

    def test_no_errors_when_analysis_enabled(self, tmp_path):
        """If run_X is true, no check is needed even without pickles."""
        cfg = _make_base_cfg(run_rmsd=True, plot_rmsd=True)
        errors = validate_plot_pickles(cfg, str(tmp_path))
        assert errors == []

    def test_no_errors_when_plot_disabled(self, tmp_path):
        """If plot_X is false, no check is needed."""
        cfg = _make_base_cfg(run_rmsd=False, plot_rmsd=False)
        errors = validate_plot_pickles(cfg, str(tmp_path))
        assert errors == []

    def test_error_when_plot_on_analysis_off_no_pickles(self, tmp_path):
        """plot_X=true, run_X=false, no pickles → should produce an error."""
        cfg = _make_base_cfg(run_rmsd=False, plot_rmsd=True)
        errors = validate_plot_pickles(cfg, str(tmp_path))
        assert len(errors) == 1
        assert 'RMSD' in errors[0]
        assert 'plot_rmsd = true' in errors[0]
        assert 'run_rmsd = false' in errors[0]
        assert str(tmp_path) in errors[0]

    def test_no_error_when_pickles_exist(self, tmp_path):
        """plot_X=true, run_X=false but pickles present → no error."""
        cfg = _make_base_cfg(run_rmsd=False, plot_rmsd=True)
        # Create a matching pickle file in the per-analysis subdir
        rmsd_dir = tmp_path / 'rmsd'
        rmsd_dir.mkdir()
        pkl = rmsd_dir / 'rmsd_plot_SysA_v1_rep1.pkl'
        pkl.write_bytes(pickle.dumps({'dummy': True}))
        errors = validate_plot_pickles(cfg, str(tmp_path))
        assert errors == []

    def test_multiple_errors_for_multiple_analyses(self, tmp_path):
        """Multiple analyses with missing pickles produce multiple errors."""
        cfg = _make_base_cfg(
            run_rmsd=False, plot_rmsd=True,
            run_rog=False,  plot_rog=True,
        )
        errors = validate_plot_pickles(cfg, str(tmp_path))
        assert len(errors) == 2
        analysis_names = ' '.join(errors)
        assert 'RMSD' in analysis_names
        assert 'Radius of Gyration' in analysis_names

    def test_error_message_contains_expected_pattern(self, tmp_path):
        """The error message should list the expected glob pattern."""
        cfg = _make_base_cfg(run_hbonds=False, plot_hbonds=True)
        errors = validate_plot_pickles(cfg, str(tmp_path))
        assert len(errors) == 1
        assert 'hbonds_plot_SysA_v1_rep1' in errors[0]

    def test_error_message_contains_resolution_hint(self, tmp_path):
        """The error should suggest enabling the analysis or providing pickles."""
        cfg = _make_base_cfg(run_rmsf=False, plot_rmsf=True)
        errors = validate_plot_pickles(cfg, str(tmp_path))
        assert len(errors) == 1
        assert 'Resolution' in errors[0]
        assert 'run_rmsf = true' in errors[0]


# ─── Time-unit and time_interval_between_frames config handling ──────────────

class TestTimeUnitConfig:
    """Tests for global time_unit parsing and conditional time_interval_between_frames requirement."""

    def test_time_unit_from_analysis_section(self, tmp_path):
        """time_unit in [analysis] should be used when present."""
        ini = _write_ini(tmp_path, "run_rmsd = true\ntime_unit = us")
        cfg = load_config(ini)
        assert cfg['time_unit'] == 'us'

    def test_time_unit_fallback_to_rmsd_section(self, tmp_path):
        """When time_unit is absent from [analysis], fall back to [rmsd]."""
        ini = _write_ini(
            tmp_path,
            "run_rmsd = true",
            rmsd_lines="time_interval_between_frames = 2.0\ntime_unit = ms"
        )
        cfg = load_config(ini)
        assert cfg['time_unit'] == 'ms'

    def test_time_unit_defaults_to_ns(self, tmp_path):
        """When time_unit is absent from both sections, default to 'ns'."""
        ini = _write_ini(tmp_path, "run_rmsd = true")
        cfg = load_config(ini)
        assert cfg['time_unit'] == 'ns'

    def test_invalid_time_unit_raises(self, tmp_path):
        """An invalid time_unit in the INI should raise ValueError."""
        ini = _write_ini(tmp_path, "run_rmsd = true\ntime_unit = hours")
        with pytest.raises(ValueError, match="Unsupported time unit"):
            load_config(ini)

    def test_mandatory_time_interval_between_frames(self, tmp_path):
        """time_interval_between_frames = none should cause sys.exit(1)."""
        ini = _write_ini(
            tmp_path,
            "run_rmsd = true",
            rmsd_lines="time_interval_between_frames = none"
        )
        with pytest.raises(SystemExit):
            load_config(ini)

    def test_time_interval_can_be_none_when_rmsd_disabled(self, tmp_path):
        """time_interval_between_frames may be none when RMSD is disabled."""
        ini = _write_ini(
            tmp_path,
            "run_rmsd = false",
            rmsd_lines="time_interval_between_frames = none"
        )
        cfg = load_config(ini)
        assert cfg['run_rmsd'] is False
        assert cfg['time_interval_between_frames'] is None

    def test_valid_time_interval_between_frames(self, tmp_path):
        """A numeric time_interval_between_frames should be parsed correctly."""
        ini = _write_ini(tmp_path, "run_rmsd = true")
        cfg = load_config(ini)
        assert cfg['time_interval_between_frames'] == 2.0

    @pytest.mark.parametrize('rmsd_lines, match', [
        ('time_interval_between_frames = 0', 'greater than zero'),
        ('time_interval_between_frames = -1', 'greater than zero'),
        ('time_interval_between_frames = 1e9', 'sanity range'),
    ])
    def test_invalid_time_interval_between_frames_raises(self, tmp_path, rmsd_lines, match):
        """Out-of-range time_interval_between_frames values should fail fast."""
        ini = _write_ini(tmp_path, "run_rmsd = true", rmsd_lines=rmsd_lines)
        with pytest.raises(ValueError, match=match):
            load_config(ini)

    def test_hbonds_time_unit_removed(self, tmp_path):
        """cfg should NOT contain 'hbonds_time_unit' anymore."""
        ini = _write_ini(tmp_path, "run_rmsd = true")
        cfg = load_config(ini)
        assert 'hbonds_time_unit' not in cfg


# ─── [plot_groups] config parsing ─────────────────────────────────────────────

class TestPlotGroupsParsing:
    """Tests for [plot_groups] section parsing in load_config."""

    def test_plot_groups_parsed_correctly(self, tmp_path):
        """A valid [plot_groups] section should produce cfg['plot_groups'] dict."""
        ini = tmp_path / "test.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\", \"SysB\"]\n"
            "variations = {\"SysA\": [\"v1\"], \"SysB\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "\n"
            "[analysis]\n"
            "run_rmsd = true\n"
            "\n"
            "[rmsd]\n"
            "time_interval_between_frames = 2.0\n"
            "\n"
            "[plot_groups]\n"
            "replicate_mode = separate\n"
            'my_group = [["SysA", "v1"], ["SysB", "v1"]]\n'
        )
        ini.write_text(content)
        cfg = load_config(str(ini))
        assert 'plot_groups' in cfg
        assert 'my_group' in cfg['plot_groups']
        assert cfg['plot_groups']['my_group'] == [('SysA', 'v1'), ('SysB', 'v1')]
        assert cfg['replicate_mode'] == 'separate'

    def test_plot_groups_replicate_mode_average(self, tmp_path):
        """replicate_mode = average should be parsed correctly."""
        ini = tmp_path / "test.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 2\n"
            "\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "\n"
            "[analysis]\n"
            "run_rmsd = true\n"
            "\n"
            "[rmsd]\n"
            "time_interval_between_frames = 2.0\n"
            "\n"
            "[plot_groups]\n"
            "replicate_mode = average\n"
            'grp1 = [["SysA", "v1"]]\n'
        )
        ini.write_text(content)
        cfg = load_config(str(ini))
        assert cfg['replicate_mode'] == 'average'

    def test_plot_groups_invalid_replicate_mode_falls_back(self, tmp_path):
        """Invalid replicate_mode should warn and fall back to 'separate'."""
        ini = tmp_path / "test.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "\n"
            "[analysis]\n"
            "run_rmsd = true\n"
            "\n"
            "[rmsd]\n"
            "time_interval_between_frames = 2.0\n"
            "\n"
            "[plot_groups]\n"
            "replicate_mode = invalid_mode\n"
            'grp1 = [["SysA", "v1"]]\n'
        )
        ini.write_text(content)
        cfg = load_config(str(ini))
        assert cfg['replicate_mode'] == 'separate'

    def test_no_plot_groups_section_gives_empty_dict(self, tmp_path):
        """Without [plot_groups], cfg['plot_groups'] should be empty dict."""
        ini = _write_ini(tmp_path, "run_rmsd = true")
        cfg = load_config(ini)
        assert cfg.get('plot_groups', {}) == {}

    def test_plot_groups_multiple_groups(self, tmp_path):
        """Multiple named groups should all be parsed."""
        ini = tmp_path / "test.ini"
        content = (
            "[systems]\n"
            "systems = [\"A\", \"B\"]\n"
            "variations = {\"A\": [\"v1\", \"v2\"], \"B\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "\n"
            "[analysis]\n"
            "run_rmsd = true\n"
            "\n"
            "[rmsd]\n"
            "time_interval_between_frames = 2.0\n"
            "\n"
            "[plot_groups]\n"
            'group1 = [["A", "v1"], ["A", "v2"]]\n'
            'group2 = [["A", "v1"], ["B", "v1"]]\n'
        )
        ini.write_text(content)
        cfg = load_config(str(ini))
        assert len(cfg['plot_groups']) == 2
        assert 'group1' in cfg['plot_groups']
        assert 'group2' in cfg['plot_groups']


# ─── Executor helper functions ────────────────────────────────────────────────

class TestExecutorHelpers:
    """Tests for the new executor helper functions."""

    def test_find_pickle_for_member_exact(self, tmp_path):
        """_find_pickle_for_member should find an exact match."""
        from executor import _find_pickle_for_member
        pkl = tmp_path / 'rmsd_plot_SysA_v1_rep1.pkl'
        pkl.write_bytes(pickle.dumps({'dummy': True}))
        result = _find_pickle_for_member(str(tmp_path), 'rmsd_plot', 'SysA', 'v1', 1)
        assert result == str(pkl)

    def test_find_pickle_for_member_with_sel_idx(self, tmp_path):
        """_find_pickle_for_member should find _selN pickles."""
        from executor import _find_pickle_for_member
        pkl = tmp_path / 'rmsd_plot_SysA_v1_rep1_sel0.pkl'
        pkl.write_bytes(pickle.dumps({'dummy': True}))
        result = _find_pickle_for_member(str(tmp_path), 'rmsd_plot', 'SysA', 'v1', 1, sel_idx=0)
        assert result == str(pkl)

    def test_find_pickle_for_member_with_ref_idx(self, tmp_path):
        """_find_pickle_for_member should find _refN pickles."""
        from executor import _find_pickle_for_member
        pkl = tmp_path / 'rmsd_plot_SysA_v1_rep1_ref0_sel0.pkl'
        pkl.write_bytes(pickle.dumps({'dummy': True}))
        result = _find_pickle_for_member(str(tmp_path), 'rmsd_plot', 'SysA', 'v1', 1, sel_idx=0, ref_idx=0)
        assert result == str(pkl)

    def test_find_pickle_for_member_missing(self, tmp_path):
        """_find_pickle_for_member should return None when no match found."""
        from executor import _find_pickle_for_member
        result = _find_pickle_for_member(str(tmp_path), 'rmsd_plot', 'SysA', 'v1', 1)
        assert result is None

    def test_detect_selection_indices_with_sels(self, tmp_path):
        """_detect_selection_indices should find _selN indices."""
        from executor import _detect_selection_indices
        (tmp_path / 'rmsd_plot_SysA_v1_rep1_sel0.pkl').touch()
        (tmp_path / 'rmsd_plot_SysA_v1_rep1_sel1.pkl').touch()
        cfg = {'systems': ['SysA'], 'variations': {'SysA': ['v1']}, 'num_replicates': 1}
        result = _detect_selection_indices(str(tmp_path), 'rmsd_plot', cfg)
        assert result == [0, 1]

    def test_detect_selection_indices_no_sels(self, tmp_path):
        """_detect_selection_indices should return [None] when no _selN found."""
        from executor import _detect_selection_indices
        (tmp_path / 'rmsd_plot_SysA_v1_rep1.pkl').touch()
        cfg = {'systems': ['SysA'], 'variations': {'SysA': ['v1']}, 'num_replicates': 1}
        result = _detect_selection_indices(str(tmp_path), 'rmsd_plot', cfg)
        assert result == [None]

    def test_get_selection_label_from_pickle(self, tmp_path):
        """_get_selection_label_from_pickle should read 'selection' and 'ref_selection' keys."""
        from executor import _get_selection_label_from_pickle
        data = {'rmsd_data': None, 'selection': 'protein and backbone', 'ref_selection': 'name CA'}
        pkl = tmp_path / 'test.pkl'
        pkl.write_bytes(pickle.dumps(data))
        sel, ref_sel = _get_selection_label_from_pickle(str(pkl))
        assert sel == 'protein and backbone'
        assert ref_sel == 'name CA'

    def test_get_selection_label_missing_key(self, tmp_path):
        """_get_selection_label_from_pickle should return ('', '') if keys missing."""
        from executor import _get_selection_label_from_pickle
        data = {'rmsd_data': None}
        pkl = tmp_path / 'test.pkl'
        pkl.write_bytes(pickle.dumps(data))
        sel, ref_sel = _get_selection_label_from_pickle(str(pkl))
        assert sel == ''
        assert ref_sel == ''

    def test_detect_ref_indices_with_refs(self, tmp_path):
        """_detect_ref_indices should find _refN indices."""
        from executor import _detect_ref_indices
        (tmp_path / 'rmsd_plot_SysA_v1_rep1_ref0.pkl').touch()
        (tmp_path / 'rmsd_plot_SysA_v1_rep1_ref1.pkl').touch()
        cfg = {'systems': ['SysA'], 'variations': {'SysA': ['v1']}, 'num_replicates': 1}
        result = _detect_ref_indices(str(tmp_path), 'rmsd_plot', cfg)
        assert result == [0, 1]

    def test_detect_ref_indices_no_refs(self, tmp_path):
        """_detect_ref_indices should return [None] when no _refN found."""
        from executor import _detect_ref_indices
        (tmp_path / 'rmsd_plot_SysA_v1_rep1.pkl').touch()
        cfg = {'systems': ['SysA'], 'variations': {'SysA': ['v1']}, 'num_replicates': 1}
        result = _detect_ref_indices(str(tmp_path), 'rmsd_plot', cfg)
        assert result == [None]


# ─── ref_selection parsing ────────────────────────────────────────────────────

class TestRefSelectionParsing:
    """Tests for ref_selection list parsing in load_config."""

    def test_ref_selection_json_list(self, tmp_path):
        """ref_selection as a JSON list should be parsed into a list."""
        ini = _write_ini(tmp_path, "run_rmsd = true")
        # Overwrite the INI with a ref_selection JSON list
        ini_path = tmp_path / "test.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "\n"
            "[analysis]\n"
            "run_rmsd = true\n"
            "\n"
            "[selections]\n"
            'ref_selection = ["protein and backbone", "name CA"]\n'
            "\n"
            "[rmsd]\n"
            "time_interval_between_frames = 2.0\n"
        )
        ini_path.write_text(content)
        cfg = load_config(str(ini_path))
        assert cfg['ref_selection'] == ['protein and backbone', 'name CA']

    def test_ref_selection_plain_string(self, tmp_path):
        """ref_selection as a plain string should be wrapped in a list."""
        ini = _write_ini(tmp_path, "run_rmsd = true")
        cfg = load_config(ini)
        # Default fallback is 'protein and backbone'
        assert cfg['ref_selection'] == ['protein and backbone']

    def test_ref_selection_single_json_string(self, tmp_path):
        """ref_selection as a single JSON string should be wrapped in a list."""
        ini_path = tmp_path / "test.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "\n"
            "[analysis]\n"
            "run_rmsd = true\n"
            "\n"
            "[selections]\n"
            'ref_selection = backbone\n'
            "\n"
            "[rmsd]\n"
            "time_interval_between_frames = 2.0\n"
        )
        ini_path.write_text(content)
        cfg = load_config(str(ini_path))
        assert cfg['ref_selection'] == ['backbone']


# ─── wrap_selection config parsing ────────────────────────────────────────────

class TestWrapSelectionParsing:
    """Tests for wrap_selection parsing in load_config."""

    def test_wrap_selection_defaults_to_auto(self, tmp_path):
        """When wrap_selection is absent, it should default to 'auto'."""
        ini = _write_ini(tmp_path, "run_rmsd = true")
        cfg = load_config(ini)
        assert cfg['wrap_selection'] == 'auto'

    def test_wrap_selection_explicit_auto(self, tmp_path):
        """wrap_selection = auto should be parsed as the string 'auto'."""
        ini_path = tmp_path / "test.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "\n"
            "[analysis]\n"
            "run_rmsd = true\n"
            "\n"
            "[selections]\n"
            "wrap_selection = auto\n"
            "\n"
            "[rmsd]\n"
            "time_interval_between_frames = 2.0\n"
        )
        ini_path.write_text(content)
        cfg = load_config(str(ini_path))
        assert cfg['wrap_selection'] == 'auto'

    def test_wrap_selection_none_gives_python_none(self, tmp_path):
        """wrap_selection = none should be parsed as Python None."""
        ini_path = tmp_path / "test.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "\n"
            "[analysis]\n"
            "run_rmsd = true\n"
            "\n"
            "[selections]\n"
            "wrap_selection = none\n"
            "\n"
            "[rmsd]\n"
            "time_interval_between_frames = 2.0\n"
        )
        ini_path.write_text(content)
        cfg = load_config(str(ini_path))
        assert cfg['wrap_selection'] is None

    def test_wrap_selection_custom_string(self, tmp_path):
        """A custom MDAnalysis selection string should be stored verbatim."""
        ini_path = tmp_path / "test.ini"
        custom = "not (protein or resname GR0 or resname UQ6)"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "\n"
            "[analysis]\n"
            "run_rmsd = true\n"
            "\n"
            "[selections]\n"
            f"wrap_selection = {custom}\n"
            "\n"
            "[rmsd]\n"
            "time_interval_between_frames = 2.0\n"
        )
        ini_path.write_text(content)
        cfg = load_config(str(ini_path))
        assert cfg['wrap_selection'] == custom


# ─── setup_workdir tests ─────────────────────────────────────────────────────

class TestSetupWorkdir:
    def _make_cfg(self, tmp_path):
        """Create a minimal cfg dict with a fake input directory."""
        input_dir = tmp_path / 'input'
        input_dir.mkdir()
        sys_dir = input_dir / 'SysA' / 'v1'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'SysA_system_v1.top').write_bytes(b'fake topology')
        (sys_dir / 'SysA_production_v1_rep_1.dcd').write_bytes(b'fake traj')

        output_dir = tmp_path / 'results'

        return {
            'systems': ['SysA'],
            'variations': {'SysA': ['v1']},
            'num_replicates': 1,
            'input_dir': str(input_dir),
            'output_dir': str(output_dir),
        }

    def test_creates_per_analysis_subdirs(self, tmp_path):
        cfg = self._make_cfg(tmp_path)
        work_dir = setup_workdir(cfg)

        for subdir in _ANALYSIS_SUBDIRS.values():
            assert os.path.isdir(os.path.join(work_dir, subdir)), \
                f"Missing per-analysis subdir: {subdir}"

    def test_creates_data_symlinks(self, tmp_path):
        cfg = self._make_cfg(tmp_path)
        work_dir = setup_workdir(cfg)

        data_dir = os.path.join(work_dir, 'SysA', 'v1')
        assert os.path.isdir(data_dir)
        assert os.path.islink(os.path.join(data_dir, 'SysA_system_v1.top'))

    def test_idempotent(self, tmp_path):
        cfg = self._make_cfg(tmp_path)
        work_dir1 = setup_workdir(cfg)
        work_dir2 = setup_workdir(cfg)
        assert work_dir1 == work_dir2


# ─── _clean_ephemeral_files tests ────────────────────────────────────────────

# ─── [parallelization] config parsing ─────────────────────────────────────────

class TestParallelizationParsing:
    """Tests for [parallelization] section parsing in load_config."""

    def test_defaults_to_serial_no_workers(self, tmp_path):
        """Without [parallelization], backend='serial', n_workers=None."""
        ini = _write_ini(tmp_path, "run_rmsd = true")
        cfg = load_config(ini)
        assert cfg['parallel_backend'] == 'serial'
        assert cfg['n_workers'] is None

    def test_explicit_multiprocessing(self, tmp_path):
        """[parallelization] backend = multiprocessing should be parsed."""
        ini_path = tmp_path / "test.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "[analysis]\n"
            "run_rmsd = true\n"
            "[rmsd]\n"
            "time_interval_between_frames = 2.0\n"
            "[parallelization]\n"
            "backend = multiprocessing\n"
            "n_workers = 4\n"
        )
        ini_path.write_text(content)
        cfg = load_config(str(ini_path))
        assert cfg['parallel_backend'] == 'multiprocessing'
        assert cfg['n_workers'] == 4

    def test_n_workers_none_string(self, tmp_path):
        """n_workers = none should be parsed as Python None."""
        ini_path = tmp_path / "test.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "[analysis]\n"
            "run_rmsd = true\n"
            "[rmsd]\n"
            "time_interval_between_frames = 2.0\n"
            "[parallelization]\n"
            "backend = serial\n"
            "n_workers = none\n"
        )
        ini_path.write_text(content)
        cfg = load_config(str(ini_path))
        assert cfg['n_workers'] is None

    def test_invalid_backend_exits(self, tmp_path):
        """An invalid backend should cause sys.exit(1)."""
        ini_path = tmp_path / "test.ini"
        content = (
            "[systems]\n"
            "systems = [\"SysA\"]\n"
            "variations = {\"SysA\": [\"v1\"]}\n"
            "num_replicates = 1\n"
            "[paths]\n"
            f"input_dir = {tmp_path}\n"
            f"output_dir = {tmp_path / 'results'}\n"
            "[analysis]\n"
            "run_rmsd = true\n"
            "[rmsd]\n"
            "time_interval_between_frames = 2.0\n"
            "[parallelization]\n"
            "backend = mpi\n"
        )
        ini_path.write_text(content)
        with pytest.raises(SystemExit):
            load_config(str(ini_path))


class TestCleanEphemeralFiles:
    def test_removes_rmsfit_files(self, tmp_path):
        (tmp_path / 'rmsfit_SysA_test.dcd').write_bytes(b'aligned')
        (tmp_path / 'rmsfit_SysA_test.xtc').write_bytes(b'aligned')
        (tmp_path / 'rmsd_plot_SysA_v1_rep1.pkl').write_bytes(b'keep')

        removed = _clean_ephemeral_files(str(tmp_path))
        assert len(removed) == 2
        assert not os.path.exists(str(tmp_path / 'rmsfit_SysA_test.dcd'))
        assert os.path.exists(str(tmp_path / 'rmsd_plot_SysA_v1_rep1.pkl'))

    def test_no_rmsfit_no_error(self, tmp_path):
        removed = _clean_ephemeral_files(str(tmp_path))
        assert removed == []


# ─── _collect_pickles with analysis_type ──────────────────────────────────────

class TestCollectPicklesSubdirs:
    def test_finds_in_analysis_subdir(self, tmp_path):
        rmsd_dir = tmp_path / 'rmsd'
        rmsd_dir.mkdir()
        (rmsd_dir / 'rmsd_plot_SysA_v1_rep1.pkl').write_bytes(b'x')

        cfg = _make_base_cfg()
        result = _collect_pickles(str(tmp_path), 'rmsd_plot', cfg, analysis_type='rmsd')
        assert len(result) == 1

    def test_ignores_root_when_analysis_type_set(self, tmp_path):
        # Pickle in root should not be found when analysis_type is specified
        (tmp_path / 'rmsd').mkdir()
        (tmp_path / 'rmsd_plot_SysA_v1_rep1.pkl').write_bytes(b'x')

        cfg = _make_base_cfg()
        result = _collect_pickles(str(tmp_path), 'rmsd_plot', cfg, analysis_type='rmsd')
        assert len(result) == 0

    def test_no_analysis_type_searches_root(self, tmp_path):
        (tmp_path / 'rmsd_plot_SysA_v1_rep1.pkl').write_bytes(b'x')

        cfg = _make_base_cfg()
        result = _collect_pickles(str(tmp_path), 'rmsd_plot', cfg, analysis_type=None)
        assert len(result) == 1

    def test_ignores_legacy_duplicate_selection_pickle(self, tmp_path):
        rmsd_dir = tmp_path / 'rmsd'
        rmsd_dir.mkdir()

        canonical = rmsd_dir / 'rmsd_plot_SysA_v1_rep1.pkl'
        duplicate = rmsd_dir / 'rmsd_plot_SysA_v1_rep1_sel0.pkl'
        extra = rmsd_dir / 'rmsd_plot_SysA_v1_rep1_sel1.pkl'

        canonical.write_bytes(pickle.dumps({'selection': 'protein', 'ref_selection': 'protein and backbone'}))
        duplicate.write_bytes(pickle.dumps({'selection': 'protein', 'ref_selection': 'protein and backbone'}))
        extra.write_bytes(pickle.dumps({'selection': 'resname T44', 'ref_selection': 'protein and backbone'}))

        cfg = _make_base_cfg()
        result = _collect_pickles(
            str(tmp_path),
            'rmsd_plot',
            cfg,
            analysis_type='rmsd',
            prune_stale_duplicates=True,
        )

        assert canonical.as_posix() in result
        assert extra.as_posix() in result
        assert duplicate.as_posix() not in result
        assert len(result) == 2
        assert canonical.exists()
        assert extra.exists()
        assert not duplicate.exists()

    def test_keeps_distinct_system_pickles_with_same_selection(self, tmp_path):
        rmsd_dir = tmp_path / 'rmsd'
        rmsd_dir.mkdir()

        sysa = rmsd_dir / 'rmsd_plot_SysA_v1_rep1.pkl'
        sysb = rmsd_dir / 'rmsd_plot_SysB_v1_rep1.pkl'

        payload = {'selection': 'protein and backbone', 'ref_selection': 'protein and backbone'}
        sysa.write_bytes(pickle.dumps(payload))
        sysb.write_bytes(pickle.dumps(payload))

        cfg = _make_base_cfg(
            systems=['SysA', 'SysB'],
            variations={'SysA': ['v1'], 'SysB': ['v1']},
            num_replicates=1,
        )
        result = _collect_pickles(
            str(tmp_path),
            'rmsd_plot',
            cfg,
            analysis_type='rmsd',
            prune_stale_duplicates=False,
        )

        assert sysa.as_posix() in result
        assert sysb.as_posix() in result
        assert len(result) == 2


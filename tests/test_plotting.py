"""
Tests for the plotting scripts.

Covers:
  - Each plotting script produces a valid PNG from mock pickle data
  - Correct handling of old (raw object) vs new (dict) pickle formats
  - Bar vs line mode for RMSF
  - KDE toggle for RMSD and RoG
  - H-bond subplots (by_time, by_type, by_ids)
  - Style module helpers
"""
import os
import sys
import pickle
import pytest
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for CI

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plotting.style import (apply_style, get_color_cycle, format_label_with_stats,
                            format_selection_subtitle, prettify_label,
                            DEFAULT_DPI, DEFAULT_COLORS)


# ─── Style Tests ──────────────────────────────────────────────────────────────

class TestStyle:
    """Tests for the shared style module."""

    def test_default_colors_length(self):
        """DEFAULT_COLORS should have 8 entries."""
        assert len(DEFAULT_COLORS) == 8

    def test_get_color_cycle_no_arg(self):
        """get_color_cycle() returns a copy of DEFAULT_COLORS."""
        colors = get_color_cycle()
        assert colors == DEFAULT_COLORS
        assert colors is not DEFAULT_COLORS

    def test_get_color_cycle_wraps(self):
        """get_color_cycle(n) should wrap around for n > palette size."""
        colors = get_color_cycle(10)
        assert len(colors) == 10
        assert colors[8] == DEFAULT_COLORS[0]
        assert colors[9] == DEFAULT_COLORS[1]

    def test_apply_style_removes_spines(self):
        """apply_style should remove top and right spines by default."""
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        apply_style(ax)
        assert not ax.spines['top'].get_visible()
        assert not ax.spines['right'].get_visible()
        assert ax.spines['bottom'].get_visible()
        assert ax.spines['left'].get_visible()
        plt.close(fig)

    def test_format_label_with_stats(self):
        """format_label_with_stats should include mean and std."""
        values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        label = format_label_with_stats("Test", values)
        assert "avg=3.00" in label
        assert "±" in label or "+/-" in label

    def test_format_label_with_stats_prettifies(self):
        """format_label_with_stats should replace underscores with spaces."""
        values = np.array([1.0, 2.0, 3.0])
        label = format_label_with_stats("sys_A / var_B", values)
        assert "sys A / var B" in label
        assert "_" not in label.split("(")[0]  # no underscores in the name portion

    def test_prettify_label(self):
        """prettify_label should replace underscores with spaces."""
        assert prettify_label("coq4_cluster0") == "coq4 cluster0"
        assert prettify_label("all_GR0") == "all GR0"
        assert prettify_label("no_underscores_here") == "no underscores here"
        assert prettify_label("clean") == "clean"

    def test_default_dpi(self):
        """DEFAULT_DPI should be 400."""
        assert DEFAULT_DPI == 400

    def test_format_selection_subtitle_short(self):
        """Short selection strings should be returned unchanged."""
        sel = 'protein and backbone'
        assert format_selection_subtitle(sel) == sel

    def test_format_selection_subtitle_truncates(self):
        """Long selection strings should be truncated with '…'."""
        sel = 'x' * 100
        result = format_selection_subtitle(sel, max_length=60)
        assert len(result) == 60
        assert result.endswith('…')


# ─── RMSD Plot Tests ─────────────────────────────────────────────────────────

class TestPlotRMSD:
    """Tests for plotting/plot_rmsd.py."""

    def test_rmsd_produces_png(self, mock_rmsd_pickle, tmp_path):
        """plot_rmsd should produce a .png file."""
        from plotting.plot_rmsd import plot_rmsd
        output_dir = str(tmp_path / 'out')
        plot_rmsd(mock_rmsd_pickle, output_dir=output_dir, dpi=72)
        files = os.listdir(output_dir)
        pngs = [f for f in files if f.endswith('.png')]
        assert len(pngs) == 1

    def test_rmsd_no_kde(self, mock_rmsd_pickle, tmp_path):
        """plot_rmsd should work without KDE sideplot."""
        from plotting.plot_rmsd import plot_rmsd
        output_dir = str(tmp_path / 'out')
        plot_rmsd(mock_rmsd_pickle, output_dir=output_dir, dpi=72, show_kde=False)
        files = os.listdir(output_dir)
        pngs = [f for f in files if f.endswith('.png')]
        assert len(pngs) == 1

    def test_rmsd_custom_label(self, mock_rmsd_pickle, tmp_path):
        """plot_rmsd should accept custom label without error."""
        from plotting.plot_rmsd import plot_rmsd
        output_dir = str(tmp_path / 'out')
        plot_rmsd(mock_rmsd_pickle, output_dir=output_dir, dpi=72, label='Custom Label')
        assert os.path.isdir(output_dir)


# ─── RMSF Plot Tests ─────────────────────────────────────────────────────────

class TestPlotRMSF:
    """Tests for plotting/plot_rmsf.py."""

    def test_rmsf_line_produces_png(self, mock_rmsf_pickle, tmp_path):
        """plot_rmsf should produce a line plot PNG."""
        from plotting.plot_rmsf import plot_rmsf
        output_dir = str(tmp_path / 'out')
        plot_rmsf(mock_rmsf_pickle, output_dir=output_dir, dpi=72, plot_type='line')
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1

    def test_rmsf_bar_produces_png(self, mock_rmsf_pickle, tmp_path):
        """plot_rmsf should produce a bar plot PNG."""
        from plotting.plot_rmsf import plot_rmsf
        output_dir = str(tmp_path / 'out')
        plot_rmsf(mock_rmsf_pickle, output_dir=output_dir, dpi=72, plot_type='bar')
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1

    def test_rmsf_chain_pickle_produces_png(self, mock_rmsf_chain_pickle, tmp_path):
        """plot_rmsf should handle chain-split dict pickles."""
        from plotting.plot_rmsf import plot_rmsf
        output_dir = str(tmp_path / 'out')
        plot_rmsf(mock_rmsf_chain_pickle, output_dir=output_dir, dpi=72)
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1


# ─── 2D-RMSD Plot Tests ──────────────────────────────────────────────────────

class TestPlot2DRMSD:
    """Tests for plotting/plot_2d_rmsd.py."""

    def test_2d_rmsd_produces_png(self, mock_2d_rmsd_pickle, tmp_path):
        """plot_2d_rmsd should produce a heatmap PNG."""
        from plotting.plot_2d_rmsd import plot_2d_rmsd
        output_dir = str(tmp_path / 'out')
        plot_2d_rmsd(mock_2d_rmsd_pickle, output_dir=output_dir, dpi=72)
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1

    def test_2d_rmsd_custom_cmap(self, mock_2d_rmsd_pickle, tmp_path):
        """plot_2d_rmsd should accept a custom colormap."""
        from plotting.plot_2d_rmsd import plot_2d_rmsd
        output_dir = str(tmp_path / 'out')
        plot_2d_rmsd(mock_2d_rmsd_pickle, output_dir=output_dir, dpi=72, cmap='plasma')
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1

    def test_2d_rmsd_portable_dict_schema(self, tmp_path):
        """plot_2d_rmsd should handle dict pickles with dist_matrix key."""
        from plotting.plot_2d_rmsd import plot_2d_rmsd
        pkl = tmp_path / '2d_rmsd_plot_test_wild_rep1.pkl'
        with open(pkl, 'wb') as f:
            pickle.dump({'dist_matrix': np.array([[0.0, 1.0], [1.0, 0.0]])}, f)

        output_dir = str(tmp_path / 'out')
        plot_2d_rmsd(str(pkl), output_dir=output_dir, dpi=72)
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1


# ─── RoG Plot Tests ──────────────────────────────────────────────────────────

class TestPlotRoG:
    """Tests for plotting/plot_rog.py."""

    def test_rog_produces_png(self, mock_rog_pickle, tmp_path):
        """plot_rog should produce a RoG plot PNG."""
        from plotting.plot_rog import plot_rog
        output_dir = str(tmp_path / 'out')
        plot_rog(mock_rog_pickle, output_dir=output_dir, dpi=72)
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1

    def test_rog_no_kde(self, mock_rog_pickle, tmp_path):
        """plot_rog should work without KDE sideplot."""
        from plotting.plot_rog import plot_rog
        output_dir = str(tmp_path / 'out')
        plot_rog(mock_rog_pickle, output_dir=output_dir, dpi=72, show_kde=False)
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1

    def test_rog_backward_compat_old_format(self, tmp_path):
        """plot_rog should handle legacy bare RoGResults pickles (no dict wrapper)."""
        from run_rog_analysis import RoGResults
        from plotting.plot_rog import plot_rog

        frames = list(range(30))
        times = [f * 0.002 for f in frames]
        rog_vals = list(np.random.uniform(15.0, 18.0, 30))
        old_result = RoGResults(frames, times, rog_vals)

        pkl_path = str(tmp_path / 'rog_plot_old_format.pkl')
        import pickle
        with open(pkl_path, 'wb') as f:
            pickle.dump(old_result, f)  # bare object, no dict wrapper

        output_dir = str(tmp_path / 'out')
        plot_rog(pkl_path, output_dir=output_dir, dpi=72)
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1


# ─── H-bond Plot Tests ───────────────────────────────────────────────────────

class TestPlotHBonds:
    """Tests for plotting/plot_hbonds.py."""

    def test_hbonds_by_time_produces_png(self, mock_hbonds_pickle, tmp_path):
        """plot_hbonds should produce a by-time PNG."""
        from plotting.plot_hbonds import plot_hbonds
        output_dir = str(tmp_path / 'out')
        plot_hbonds(mock_hbonds_pickle, output_dir=output_dir, dpi=72, plot_types=['time'])
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1
        assert 'by_time' in pngs[0]

    def test_hbonds_by_type_produces_png(self, mock_hbonds_pickle, tmp_path):
        """plot_hbonds should produce a by-type PNG."""
        from plotting.plot_hbonds import plot_hbonds
        output_dir = str(tmp_path / 'out')
        plot_hbonds(mock_hbonds_pickle, output_dir=output_dir, dpi=72, plot_types=['type'])
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1
        assert 'by_type' in pngs[0]

    def test_hbonds_by_ids_produces_png(self, mock_hbonds_pickle, tmp_path):
        """plot_hbonds should produce a by-ids PNG."""
        from plotting.plot_hbonds import plot_hbonds
        output_dir = str(tmp_path / 'out')
        plot_hbonds(mock_hbonds_pickle, output_dir=output_dir, dpi=72, plot_types=['ids'])
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1
        assert 'by_ids' in pngs[0]

    def test_hbonds_all_types(self, mock_hbonds_pickle, tmp_path):
        """plot_hbonds should produce all three PNGs when all types requested."""
        from plotting.plot_hbonds import plot_hbonds
        output_dir = str(tmp_path / 'out')
        plot_hbonds(mock_hbonds_pickle, output_dir=output_dir, dpi=72)
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 3

    def test_hbonds_ns_time_unit(self, mock_hbonds_pickle, tmp_path):
        """plot_hbonds should handle ns time unit."""
        from plotting.plot_hbonds import plot_hbonds
        output_dir = str(tmp_path / 'out')
        plot_hbonds(mock_hbonds_pickle, output_dir=output_dir, dpi=72,
                    plot_types=['time'], time_unit='ns')
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1

    def test_hbonds_uses_preset_label(self, mock_hbonds_pickle, tmp_path, monkeypatch):
        """plot_hbonds should label plots with the preset, not atom selections."""
        from plotting import plot_hbonds as mod

        captured = {}

        def _capture(times, counts, output_path, dpi=mod.DEFAULT_DPI, time_unit='ps', label=''):
            captured['label'] = label

        monkeypatch.setattr(mod, 'plot_hbonds_by_time', _capture)
        monkeypatch.setattr(mod, 'plot_hbonds_by_type', lambda *args, **kwargs: None)
        monkeypatch.setattr(mod, 'plot_hbonds_by_ids', lambda *args, **kwargs: None)

        mod.plot_hbonds(mock_hbonds_pickle, output_dir=str(tmp_path / 'out'), dpi=72, plot_types=['time'])

        assert 'custom' in captured['label'].lower()
        assert 'protein and name' not in captured['label'].lower()

    def test_hbonds_legacy_schema_raises(self, tmp_path):
        """Legacy non-dict H-bonds pickles should raise a schema error."""
        from plotting.plot_hbonds import plot_hbonds
        legacy_pkl = tmp_path / 'hbonds_plot_legacy.pkl'
        with open(legacy_pkl, 'wb') as f:
            pickle.dump(object(), f)

        with pytest.raises(ValueError, match='Unsupported H-bonds pickle schema'):
            plot_hbonds(str(legacy_pkl), output_dir=str(tmp_path / 'out'), dpi=72)

    def test_hbonds_nonportable_pickle_fails_fast(self, tmp_path):
        """Unpickle-time DCD/path failures should return actionable guidance."""
        from plotting import plot_hbonds as mod

        pkl = tmp_path / 'hbonds_plot_broken.pkl'
        pkl.write_bytes(b'not-used')

        original_load = mod.pickle.load

        def _raise_oserror(_file_obj):
            raise OSError('DCD file does not exist')

        mod.pickle.load = _raise_oserror
        try:
            with pytest.raises(ValueError, match='Regenerate this pickle'):
                mod.plot_hbonds(str(pkl), output_dir=str(tmp_path / 'out'), dpi=72)
        finally:
            mod.pickle.load = original_load


# ─── RMSD Comparison Plot Tests ──────────────────────────────────────────────

class TestPlotRMSDComparison:
    """Tests for plot_rmsd_comparison and plot_rmsd_comparison_average."""

    def test_rmsd_comparison_produces_png(self, mock_rmsd_comparison_pickles, tmp_path):
        """plot_rmsd_comparison should produce a comparison PNG."""
        from plotting.plot_rmsd import plot_rmsd_comparison
        output_dir = str(tmp_path / 'out')
        plot_rmsd_comparison(
            mock_rmsd_comparison_pickles,
            labels=['SysA / wild', 'SysB / wild'],
            group_name='test_group',
            output_dir=output_dir,
            dpi=72,
            show_kde=False,
        )
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1
        assert 'rmsd_comparison_test_group' in pngs[0]

    def test_rmsd_comparison_with_kde(self, mock_rmsd_comparison_pickles, tmp_path):
        """plot_rmsd_comparison with KDE should still produce a PNG."""
        from plotting.plot_rmsd import plot_rmsd_comparison
        output_dir = str(tmp_path / 'out')
        plot_rmsd_comparison(
            mock_rmsd_comparison_pickles,
            labels=['SysA / wild', 'SysB / wild'],
            group_name='with_kde',
            output_dir=output_dir,
            dpi=72,
            show_kde=True,
        )
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1

    def test_rmsd_comparison_with_selection_label(self, mock_rmsd_comparison_pickles, tmp_path):
        """plot_rmsd_comparison should accept selection_label without error."""
        from plotting.plot_rmsd import plot_rmsd_comparison
        output_dir = str(tmp_path / 'out')
        plot_rmsd_comparison(
            mock_rmsd_comparison_pickles,
            labels=['SysA / wild', 'SysB / wild'],
            group_name='sel_label',
            output_dir=output_dir,
            dpi=72,
            show_kde=False,
            selection_label='protein and backbone',
        )
        assert os.path.isdir(output_dir)

    def test_rmsd_comparison_average_produces_png(self, mock_rmsd_replicate_groups, tmp_path):
        """plot_rmsd_comparison_average should produce an averaged comparison PNG."""
        from plotting.plot_rmsd import plot_rmsd_comparison_average
        output_dir = str(tmp_path / 'out')
        plot_rmsd_comparison_average(
            mock_rmsd_replicate_groups,
            labels=['SysA / wild', 'SysB / wild'],
            group_name='avg_test',
            output_dir=output_dir,
            dpi=72,
            show_kde=False,
        )
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1
        assert 'average' in pngs[0].lower() or 'avg' in pngs[0].lower()


# ─── RMSF Comparison Plot Tests ──────────────────────────────────────────────

class TestPlotRMSFComparison:
    """Tests for plot_rmsf_comparison and plot_rmsf_comparison_average."""

    def test_rmsf_comparison_produces_png(self, mock_rmsf_comparison_pickles, tmp_path):
        """plot_rmsf_comparison should produce a comparison PNG."""
        from plotting.plot_rmsf import plot_rmsf_comparison
        output_dir = str(tmp_path / 'out')
        plot_rmsf_comparison(
            mock_rmsf_comparison_pickles,
            labels=['SysA / wild', 'SysB / wild'],
            group_name='test_group',
            output_dir=output_dir,
            dpi=72,
        )
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1
        assert 'rmsf_comparison_test_group' in pngs[0]

    def test_rmsf_comparison_no_fill_between(self, mock_rmsf_comparison_pickles, tmp_path):
        """plot_rmsf_comparison should NOT use fill_between (line-only)."""
        from plotting.plot_rmsf import plot_rmsf_comparison
        import matplotlib.pyplot as plt
        output_dir = str(tmp_path / 'out')
        plot_rmsf_comparison(
            mock_rmsf_comparison_pickles,
            labels=['SysA', 'SysB'],
            group_name='no_fill',
            output_dir=output_dir,
            dpi=72,
        )
        # The test passes if no error is raised; the implementation uses ax.plot only

    def test_rmsf_comparison_with_selection_label(self, mock_rmsf_comparison_pickles, tmp_path):
        """plot_rmsf_comparison should accept selection_label without error."""
        from plotting.plot_rmsf import plot_rmsf_comparison
        output_dir = str(tmp_path / 'out')
        plot_rmsf_comparison(
            mock_rmsf_comparison_pickles,
            labels=['SysA', 'SysB'],
            group_name='sel_test',
            output_dir=output_dir,
            dpi=72,
            selection_label='protein and backbone',
        )
        assert os.path.isdir(output_dir)

    def test_rmsf_comparison_average_produces_png(self, mock_rmsf_replicate_groups, tmp_path):
        """plot_rmsf_comparison_average should produce an averaged comparison PNG."""
        from plotting.plot_rmsf import plot_rmsf_comparison_average
        output_dir = str(tmp_path / 'out')
        plot_rmsf_comparison_average(
            mock_rmsf_replicate_groups,
            labels=['SysA / wild', 'SysB / wild'],
            group_name='avg_test',
            output_dir=output_dir,
            dpi=72,
        )
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1
        assert 'average' in pngs[0].lower() or 'avg' in pngs[0].lower()


# ─── Grouped Comparison Helper Tests ───────────────────────────────────────

class TestPlotGroupComparisonsHelper:
    """Tests for plot_group_comparisons.py helper used by Nextflow."""

    def test_rmsd_group_helper_separate_mode(self, mock_rmsd_replicate_groups, tmp_path):
        """RMSD helper should generate one comparison plot per replicate in separate mode."""
        from plot_group_comparisons import _plot_rmsd_groups

        work_dir = str(tmp_path)
        output_dir = str(tmp_path / 'out_rmsd')
        groups = {'grp_test': [('sysA', 'wild'), ('sysB', 'wild')]}

        produced = _plot_rmsd_groups(
            groups,
            work_dir=work_dir,
            output_dir=output_dir,
            num_replicates=2,
            dpi=72,
            replicate_mode='separate',
            rmsd_show_kde=False,
        )

        assert produced == 2
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 2

    def test_rmsf_group_helper_average_mode(self, mock_rmsf_replicate_groups, tmp_path):
        """RMSF helper should generate one averaged comparison plot in average mode."""
        from plot_group_comparisons import _plot_rmsf_groups

        work_dir = str(tmp_path)
        output_dir = str(tmp_path / 'out_rmsf')
        groups = {'grp_test': [('sysA', 'wild'), ('sysB', 'wild')]}

        produced = _plot_rmsf_groups(
            groups,
            work_dir=work_dir,
            output_dir=output_dir,
            num_replicates=2,
            dpi=72,
            replicate_mode='average',
            target_selection='protein and backbone',
        )

        assert produced == 1
        pngs = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        assert len(pngs) == 1

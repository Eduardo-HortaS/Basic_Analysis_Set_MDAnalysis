"""
Tests for utils.py — transform_trajectory, align_trajectory, and time-unit utilities.

These tests use mocked MDAnalysis objects to verify the utility functions
without requiring actual trajectory data.
"""
import os
import sys
import pytest
import numpy as np
from unittest.mock import MagicMock, patch, call

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import (
    transform_trajectory,
    align_trajectory,
    build_complex_selection,
    TIME_UNITS,
    SUPPORTED_TIME_UNITS,
    TIME_UNIT_LABELS,
    validate_time_unit,
    convert_time_from_ps,
    time_unit_label,
)


class TestTransformTrajectory:
    """Tests for the transform_trajectory function."""

    def test_adds_transformations_to_universe_without_ligand(self):
        """Verify that three transformations (unwrap, center, wrap) are added when no ligand."""
        mock_universe = MagicMock()
        mock_center = MagicMock()
        mock_wrap = MagicMock()

        with patch('utils.trans') as mock_trans:
            mock_trans.unwrap.return_value = 'unwrap_result'
            mock_trans.center_in_box.return_value = 'center_result'
            mock_trans.wrap.return_value = 'wrap_result'

            result = transform_trajectory(mock_universe, mock_center, mock_wrap)

            mock_trans.unwrap.assert_called_once_with(mock_universe.atoms)
            mock_trans.center_in_box.assert_called_once_with(mock_center, wrap=False)
            mock_trans.wrap.assert_called_once_with(mock_wrap, compound='fragments')
            mock_universe.trajectory.add_transformations.assert_called_once_with(
                'unwrap_result', 'center_result', 'wrap_result'
            )
            assert result is mock_universe

    def test_adds_ligand_wrap_step_when_provided(self):
        """Verify that four transformations are added when ligand_selection is given."""
        mock_universe = MagicMock()
        mock_center = MagicMock()
        mock_wrap = MagicMock()
        mock_ligand = MagicMock()

        with patch('utils.trans') as mock_trans:
            mock_trans.unwrap.return_value = 'unwrap_result'
            mock_trans.center_in_box.return_value = 'center_result'
            # wrap is called twice: once for ligand, once for solvent
            mock_trans.wrap.side_effect = ['ligand_wrap_result', 'solvent_wrap_result']

            result = transform_trajectory(mock_universe, mock_center, mock_wrap,
                                          ligand_selection=mock_ligand)

            mock_trans.unwrap.assert_called_once_with(mock_universe.atoms)
            mock_trans.center_in_box.assert_called_once_with(mock_center, wrap=False)
            assert mock_trans.wrap.call_count == 2
            mock_trans.wrap.assert_any_call(mock_ligand, compound='fragments')
            mock_trans.wrap.assert_any_call(mock_wrap, compound='fragments')
            mock_universe.trajectory.add_transformations.assert_called_once_with(
                'unwrap_result', 'center_result', 'ligand_wrap_result', 'solvent_wrap_result'
            )
            assert result is mock_universe

    def test_returns_universe(self):
        """Verify the function returns the universe object."""
        mock_universe = MagicMock()
        with patch('utils.trans'):
            result = transform_trajectory(mock_universe, MagicMock(), MagicMock())
        assert result is mock_universe


class TestBuildComplexSelection:
    """Tests for the build_complex_selection function."""

    def test_auto_selects_protein_and_not_protein(self):
        """In 'auto' mode, center_ag should be 'protein', ligand_ag None, wrap_ag 'not protein'."""
        mock_universe = MagicMock()
        mock_protein = MagicMock(name='protein_ag')
        mock_not_protein = MagicMock(name='not_protein_ag')

        def select_side_effect(sel):
            if sel == 'protein':
                return mock_protein
            elif sel == 'not protein':
                return mock_not_protein
            return MagicMock()

        mock_universe.select_atoms.side_effect = select_side_effect

        center_ag, ligand_ag, wrap_ag = build_complex_selection(mock_universe, wrap_selection='auto')
        assert center_ag is mock_protein
        assert ligand_ag is None
        assert wrap_ag is mock_not_protein

    def test_none_returns_none_none_none(self):
        """When wrap_selection is None, all three return values should be None."""
        mock_universe = MagicMock()
        center_ag, ligand_ag, wrap_ag = build_complex_selection(mock_universe, wrap_selection=None)
        assert center_ag is None
        assert ligand_ag is None
        assert wrap_ag is None
        # Universe should not have been queried
        mock_universe.select_atoms.assert_not_called()

    def test_custom_selection_string(self):
        """A custom wrap_selection: center_ag is protein, ligand_ag is the complement of both."""
        mock_universe = MagicMock()
        mock_protein = MagicMock(name='protein_ag')
        mock_rest = MagicMock(name='rest_ag')
        mock_ligand = MagicMock(name='ligand_ag')
        mock_ligand.__len__ = MagicMock(return_value=50)

        call_count = [0]
        def select_side_effect(sel):
            call_count[0] += 1
            if sel == 'protein':
                return mock_protein
            else:
                return mock_rest

        mock_universe.select_atoms.side_effect = select_side_effect
        # universe.atoms - center_ag - wrap_ag = ligand_ag
        mock_intermediate = MagicMock()
        mock_intermediate.__sub__ = MagicMock(return_value=mock_ligand)
        mock_universe.atoms.__sub__ = MagicMock(return_value=mock_intermediate)

        custom_sel = 'not (protein or resname GR0 or resname UQ6)'
        center_ag, ligand_ag, wrap_ag = build_complex_selection(mock_universe, wrap_selection=custom_sel)

        assert center_ag is mock_protein
        assert wrap_ag is mock_rest
        assert ligand_ag is mock_ligand

    def test_custom_selection_no_ligand(self):
        """When custom selection covers everything except protein, ligand_ag should be None."""
        mock_universe = MagicMock()
        mock_protein = MagicMock(name='protein_ag')
        mock_rest = MagicMock(name='rest_ag')
        mock_empty = MagicMock(name='empty_ag')
        mock_empty.__len__ = MagicMock(return_value=0)

        def select_side_effect(sel):
            if sel == 'protein':
                return mock_protein
            else:
                return mock_rest

        mock_universe.select_atoms.side_effect = select_side_effect
        mock_intermediate = MagicMock()
        mock_intermediate.__sub__ = MagicMock(return_value=mock_empty)
        mock_universe.atoms.__sub__ = MagicMock(return_value=mock_intermediate)

        center_ag, ligand_ag, wrap_ag = build_complex_selection(mock_universe, wrap_selection='not protein')
        assert center_ag is mock_protein
        assert ligand_ag is None
        assert wrap_ag is mock_rest

    def test_default_is_auto(self):
        """The default wrap_selection should be 'auto'."""
        mock_universe = MagicMock()
        mock_protein = MagicMock()
        mock_not_protein = MagicMock()

        def select_side_effect(sel):
            if sel == 'protein':
                return mock_protein
            elif sel == 'not protein':
                return mock_not_protein
            return MagicMock()

        mock_universe.select_atoms.side_effect = select_side_effect

        center_ag, ligand_ag, wrap_ag = build_complex_selection(mock_universe)
        assert center_ag is mock_protein
        assert ligand_ag is None
        assert wrap_ag is mock_not_protein


class TestAlignTrajectory:
    """Tests for the align_trajectory function."""

    def test_rmsf_alignment_creates_average_structure(self):
        """Verify RMSF alignment computes average structure then aligns."""
        mock_universe = MagicMock()
        mock_avg_result = MagicMock()
        mock_ref_universe = MagicMock()
        mock_avg_result.results.universe = mock_ref_universe

        with patch('utils.align') as mock_align:
            mock_align.AverageStructure.return_value.run.return_value = mock_avg_result
            mock_align.AlignTraj.return_value.run.return_value = None

            align_trajectory(mock_universe, 'protein', 'RMSF', 'sys1', 'wild', 1, 'dcd', 100)

            mock_align.AverageStructure.assert_called_once_with(
                mock_universe, mock_universe, select='protein', ref_frame=0
            )
            mock_align.AlignTraj.assert_called_once_with(
                mock_universe, mock_ref_universe, select='protein',
                filename='rmsfit_sys1_production_wild_reduced_rep1.dcd',
                in_memory=False
            )

    def test_2d_rmsd_alignment_self_aligns(self):
        """Verify 2D-RMSD alignment uses self-alignment."""
        mock_universe = MagicMock()

        with patch('utils.align') as mock_align:
            mock_align.AlignTraj.return_value.run.return_value = None

            align_trajectory(mock_universe, 'backbone', '2D-RMSD', 'sys1', 'wild', 1, 'dcd', 50)

            mock_align.AlignTraj.assert_called_once_with(
                mock_universe, mock_universe, select='backbone',
                filename='rmsfit_sys1_production_wild_reduced_rep1.dcd',
                in_memory=False
            )

    def test_filename_format(self):
        """Verify the aligned trajectory filename format is correct."""
        mock_universe = MagicMock()

        with patch('utils.align') as mock_align:
            mock_align.AlignTraj.return_value.run.return_value = None

            align_trajectory(mock_universe, 'sel', '2D-RMSD', 'mysys', 'k84r', 3, 'xtc', 200)

            expected_filename = 'rmsfit_mysys_production_k84r_reduced_rep3.xtc'
            _, kwargs = mock_align.AlignTraj.call_args
            assert kwargs['filename'] == expected_filename


# ─── Time-unit utility tests ─────────────────────────────────────────────────

class TestTimeUnitConstants:
    """Tests for the time-unit constants exposed by utils.py."""

    def test_supported_time_units_tuple(self):
        """SUPPORTED_TIME_UNITS should list the six canonical units."""
        assert SUPPORTED_TIME_UNITS == ('fs', 'ps', 'ns', 'us', 'ms', 's')

    def test_time_units_dict_contains_all_canonical_units(self):
        """TIME_UNITS dict should have entries for every canonical unit."""
        for unit in SUPPORTED_TIME_UNITS:
            assert unit in TIME_UNITS

    def test_time_units_dict_contains_unicode_mu_alias(self):
        """TIME_UNITS should also accept 'μs' as an alias for 'us'."""
        assert 'μs' in TIME_UNITS
        assert TIME_UNITS['μs'] == TIME_UNITS['us']

    def test_time_unit_labels_for_us(self):
        """The pretty-print label for 'us' should be 'µs'."""
        assert TIME_UNIT_LABELS['us'] == 'µs'


class TestValidateTimeUnit:
    """Tests for validate_time_unit()."""

    def test_all_canonical_units_pass(self):
        """Every canonical unit should pass without raising."""
        for unit in SUPPORTED_TIME_UNITS:
            validate_time_unit(unit)  # should not raise

    def test_unicode_mu_passes(self):
        """The 'μs' alias should be accepted."""
        validate_time_unit('μs')

    @pytest.mark.parametrize("bad_unit", ['hours', 'days', '', 'nanoseconds', 'NS', 'Ns'])
    def test_invalid_units_raise(self, bad_unit):
        """Unsupported strings should raise ValueError."""
        with pytest.raises(ValueError, match="Unsupported time unit"):
            validate_time_unit(bad_unit)


class TestConvertTimeFromPs:
    """Tests for convert_time_from_ps()."""

    def test_ps_to_ns(self):
        """1000 ps → 1 ns."""
        assert convert_time_from_ps(1000.0, 'ns') == pytest.approx(1.0)

    def test_ps_to_us(self):
        """1e6 ps → 1 µs."""
        assert convert_time_from_ps(1e6, 'us') == pytest.approx(1.0)

    def test_ps_to_fs(self):
        """1 ps → 1000 fs."""
        assert convert_time_from_ps(1.0, 'fs') == pytest.approx(1000.0)

    def test_ps_identity(self):
        """ps → ps should be identity."""
        assert convert_time_from_ps(42.0, 'ps') == pytest.approx(42.0)

    def test_ps_to_ms(self):
        """1e9 ps → 1 ms."""
        assert convert_time_from_ps(1e9, 'ms') == pytest.approx(1.0)

    def test_ps_to_s(self):
        """1e12 ps → 1 s."""
        assert convert_time_from_ps(1e12, 's') == pytest.approx(1.0)

    def test_numpy_array(self):
        """Should operate element-wise on numpy arrays."""
        arr = np.array([0.0, 1000.0, 2000.0])
        result = convert_time_from_ps(arr, 'ns')
        np.testing.assert_array_almost_equal(result, [0.0, 1.0, 2.0])

    def test_invalid_unit_raises(self):
        """Should raise ValueError for unrecognised unit."""
        with pytest.raises(ValueError):
            convert_time_from_ps(100.0, 'lightyears')


class TestTimeUnitLabel:
    """Tests for time_unit_label()."""

    def test_us_label(self):
        """'us' should display as 'µs'."""
        assert time_unit_label('us') == 'µs'

    def test_ns_label(self):
        """'ns' should display as 'ns'."""
        assert time_unit_label('ns') == 'ns'

    def test_unknown_falls_back(self):
        """An unknown string should be returned as-is."""
        assert time_unit_label('fortnights') == 'fortnights'

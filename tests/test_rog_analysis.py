"""
Tests for run_rog_analysis.py — Radius of Gyration analysis.

Covers:
  - RoGResults class structure
  - RoG pickle output format
  - Skip-if-exists behavior
"""
import os
import sys
import pickle
import pytest
import numpy as np
from unittest.mock import MagicMock, patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from run_rog_analysis import RoGResults, run_rog_analysis


class TestRoGResults:
    """Tests for the RoGResults container class."""

    def test_rog_data_shape(self):
        """RoGResults.rog_data should have shape (n, 3)."""
        frames = [0, 1, 2, 3, 4]
        times = [0.0, 0.1, 0.2, 0.3, 0.4]
        rog_vals = [15.0, 15.5, 16.0, 15.8, 15.2]
        result = RoGResults(frames, times, rog_vals)
        assert result.rog_data.shape == (5, 3)

    def test_rog_data_columns(self):
        """Columns should be [frame, time, rog]."""
        frames = [0, 1, 2]
        times = [0.0, 0.1, 0.2]
        rog_vals = [15.0, 15.5, 16.0]
        result = RoGResults(frames, times, rog_vals)
        np.testing.assert_array_equal(result.rog_data[:, 0], frames)
        np.testing.assert_array_almost_equal(result.rog_data[:, 1], times)
        np.testing.assert_array_almost_equal(result.rog_data[:, 2], rog_vals)

    def test_rog_data_pickleable(self, tmp_path):
        """RoGResults should be pickleable and unpickleable."""
        result = RoGResults([0, 1], [0.0, 0.1], [15.0, 15.5])
        pkl = tmp_path / 'test.pkl'
        with open(pkl, 'wb') as f:
            pickle.dump(result, f)
        with open(pkl, 'rb') as f:
            loaded = pickle.load(f)
        assert hasattr(loaded, 'rog_data')
        assert loaded.rog_data.shape == (2, 3)


class TestRunRoGAnalysis:
    """Tests for the run_rog_analysis function."""

    @patch('run_rog_analysis.mda')
    def test_rog_pickle_created(self, mock_mda, tmp_path):
        """Should create a RoG pickle file."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()

        # Mock universe and trajectory
        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe

        # Mock select_atoms for protein/not protein
        mock_protein = MagicMock()
        mock_not_protein = MagicMock()
        mock_selected = MagicMock()
        mock_selected.radius_of_gyration.return_value = 16.0

        def select_side_effect(sel):
            if sel == 'protein':
                return mock_protein
            elif sel == 'not protein':
                return mock_not_protein
            else:
                return mock_selected
        mock_universe.select_atoms.side_effect = select_side_effect

        # Mock trajectory iteration
        mock_ts1 = MagicMock()
        mock_ts1.frame = 0
        mock_ts1.time = 0.0
        mock_ts2 = MagicMock()
        mock_ts2.frame = 1
        mock_ts2.time = 2000.0
        mock_universe.trajectory.__getitem__ = MagicMock(return_value=[mock_ts1, mock_ts2])

        run_rog_analysis(
            systems=['sys1'],
            variations={'sys1': ['wild']},
            num_replicates=1,
            start_frame=0,
            traj_format='dcd',
            selection='protein and backbone'
        )

        pkl = tmp_path / 'rog_plot_sys1_wild_rep1.pkl'
        assert pkl.exists()

        with open(pkl, 'rb') as f:
            data = pickle.load(f)
        assert isinstance(data, dict)
        assert 'rog_obj' in data
        assert hasattr(data['rog_obj'], 'rog_data')
        assert data['rog_obj'].rog_data.shape[1] == 3
        assert data.get('selection') == 'protein and backbone'

    @patch('run_rog_analysis.mda')
    def test_rog_skip_if_exists(self, mock_mda, tmp_path, capsys):
        """Should skip RoG analysis if pickle already exists."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rog_plot_sys1_wild_rep1.pkl').touch()

        run_rog_analysis(
            systems=['sys1'],
            variations={'sys1': ['wild']},
            num_replicates=1,
            start_frame=0,
            traj_format='dcd'
        )

        captured = capsys.readouterr()
        assert "Skipping" in captured.out
        mock_mda.Universe.assert_not_called()

    @patch('run_rog_analysis.mda')
    def test_rog_pickle_is_dict_with_metadata(self, mock_mda, tmp_path):
        """RoG pickle should be a dict with 'rog_obj' and 'time_unit' keys."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe

        mock_protein = MagicMock()
        mock_not_protein = MagicMock()
        mock_selected = MagicMock()
        mock_selected.radius_of_gyration.return_value = 16.0

        def select_side_effect(sel):
            if sel == 'protein':
                return mock_protein
            elif sel == 'not protein':
                return mock_not_protein
            else:
                return mock_selected
        mock_universe.select_atoms.side_effect = select_side_effect

        mock_ts1 = MagicMock()
        mock_ts1.frame = 0
        mock_ts1.time = 0.0
        mock_ts2 = MagicMock()
        mock_ts2.frame = 1
        mock_ts2.time = 2000.0
        mock_universe.trajectory.__getitem__ = MagicMock(return_value=[mock_ts1, mock_ts2])

        run_rog_analysis(
            systems=['sys1'],
            variations={'sys1': ['wild']},
            num_replicates=1,
            start_frame=0,
            traj_format='dcd',
            time_unit='ns',
        )

        pkl = tmp_path / 'rog_plot_sys1_wild_rep1.pkl'
        assert pkl.exists()

        with open(pkl, 'rb') as f:
            data = pickle.load(f)
        assert isinstance(data, dict)
        assert 'rog_obj' in data
        assert data['time_unit'] == 'ns'
        assert hasattr(data['rog_obj'], 'rog_data')

    @patch('run_rog_analysis.mda')
    def test_rog_time_unit_us(self, mock_mda, tmp_path):
        """RoG with time_unit='us' should convert ps → µs."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe

        mock_protein = MagicMock()
        mock_not_protein = MagicMock()
        mock_selected = MagicMock()
        mock_selected.radius_of_gyration.return_value = 16.0

        def select_side_effect(sel):
            if sel == 'protein':
                return mock_protein
            elif sel == 'not protein':
                return mock_not_protein
            else:
                return mock_selected
        mock_universe.select_atoms.side_effect = select_side_effect

        # ts.time is in ps; 1e6 ps = 1 µs
        mock_ts = MagicMock()
        mock_ts.frame = 0
        mock_ts.time = 1e6  # 1 µs in ps
        mock_universe.trajectory.__getitem__ = MagicMock(return_value=[mock_ts])

        run_rog_analysis(
            systems=['sys1'],
            variations={'sys1': ['wild']},
            num_replicates=1,
            start_frame=0,
            traj_format='dcd',
            time_unit='us',
        )

        pkl = tmp_path / 'rog_plot_sys1_wild_rep1.pkl'
        with open(pkl, 'rb') as f:
            data = pickle.load(f)
        # time column should be 1.0 µs
        np.testing.assert_almost_equal(data['rog_obj'].rog_data[0, 1], 1.0)

    def test_rog_invalid_time_unit_raises(self):
        """An invalid time_unit should raise ValueError."""
        with pytest.raises(ValueError, match="Unsupported time unit"):
            run_rog_analysis(
                systems=['sys1'],
                variations={'sys1': ['wild']},
                num_replicates=1,
                start_frame=0,
                traj_format='dcd',
                time_unit='hours',
            )

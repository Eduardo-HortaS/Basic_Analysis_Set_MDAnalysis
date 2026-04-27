"""
Tests for run_rms_analysis.py — RMSD, RMSF, 2D-RMSD analysis functions.

Covers:
  - validate_chain_intervals logic
  - RMSD with DCD time-axis correction
  - RMSF with chain split
  - Skip-if-exists behavior
"""
import os
import sys
import pickle
import pytest
import numpy as np
from unittest.mock import MagicMock, patch, ANY

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from run_rms_analysis import validate_chain_intervals, run_rms_analysis


class TestValidateChainIntervals:
    """Tests for the validate_chain_intervals function."""

    def _make_universe(self, resids):
        """Create a mock universe with given resids in the target selection."""
        u = MagicMock()
        mock_atoms = MagicMock()
        mock_residues = MagicMock()
        mock_residues.resids = np.array(resids)
        mock_atoms.residues = mock_residues
        u.select_atoms.return_value = mock_atoms
        return u

    def test_valid_intervals(self):
        """Valid non-overlapping intervals should pass without errors."""
        u = self._make_universe(range(1, 101))
        chain_intervals = {"A": [1, 50], "B": [51, 100]}
        # Should not raise
        validate_chain_intervals(u, chain_intervals, 'protein and backbone')

    def test_invalid_start_resid_zero(self):
        """Should raise ValueError if start_resid is 0 (not 1-based)."""
        u = self._make_universe(range(1, 101))
        chain_intervals = {"A": [0, 50]}
        with pytest.raises(ValueError, match="less than 1"):
            validate_chain_intervals(u, chain_intervals, 'protein and backbone')

    def test_start_greater_than_end(self):
        """Should raise ValueError if start > end."""
        u = self._make_universe(range(1, 101))
        chain_intervals = {"A": [60, 50]}
        with pytest.raises(ValueError, match="start_resid=60 > end_resid=50"):
            validate_chain_intervals(u, chain_intervals, 'protein and backbone')

    def test_missing_resids(self):
        """Should raise ValueError if interval resids are not in trajectory."""
        u = self._make_universe(range(1, 51))
        chain_intervals = {"A": [1, 50], "B": [51, 100]}
        with pytest.raises(ValueError, match="not found in target selection"):
            validate_chain_intervals(u, chain_intervals, 'protein and backbone')

    def test_overlapping_intervals_warns(self, capsys):
        """Overlapping intervals should print a warning."""
        u = self._make_universe(range(1, 101))
        chain_intervals = {"A": [1, 60], "B": [50, 100]}
        validate_chain_intervals(u, chain_intervals, 'protein and backbone')
        captured = capsys.readouterr()
        assert "overlapping" in captured.out.lower()

    def test_gap_in_coverage_warns(self, capsys):
        """Gaps between intervals should print a warning."""
        u = self._make_universe(range(1, 101))
        chain_intervals = {"A": [1, 40], "B": [61, 100]}
        validate_chain_intervals(u, chain_intervals, 'protein and backbone')
        captured = capsys.readouterr()
        assert "not covered" in captured.out.lower()


class TestRunRmsAnalysisRMSD:
    """Tests for RMSD analysis with DCD time correction."""

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_rmsd_pickle_contains_metadata(self, mock_rms, mock_mda, mock_transform, tmp_path):
        """RMSD pickle should contain dict with time metadata."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()
        # Pre-create aligned trajectory so the new align_trajectory call is skipped
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe
        # select_atoms must return something with non-zero len
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=100)
        mock_universe.select_atoms.return_value = mock_atoms

        mock_rmsd_result = MagicMock()
        mock_rmsd_result.results.rmsd = np.column_stack([
            np.arange(10, dtype=float),
            np.arange(10, dtype=float),
            np.random.uniform(1.0, 3.0, 10)
        ])
        mock_rms.RMSD.return_value = mock_rmsd_result
        mock_rmsd_result.run.return_value = mock_rmsd_result

        captured_data = {}

        def capture_pickle_dump(obj, f):
            """Intercept pickle.dump to capture the data without actually pickling MagicMock."""
            captured_data.update(obj if isinstance(obj, dict) else {'obj': obj})
            # Write a placeholder so the file is non-empty
            import pickle as real_pickle
            real_pickle.dumps(b'placeholder')

        with patch('run_rms_analysis.pickle') as mock_pickle:
            mock_pickle.dump.side_effect = capture_pickle_dump
            run_rms_analysis(
                systems=['sys1'],
                variations={'sys1': ['wild']},
                num_replicates=1,
                analysis='RMSD',
                target_selection='protein and backbone',
                ref_selection='protein and backbone',
                start_frame=0,
                traj_format='dcd',
                time_interval_between_frames=2.0,
                time_unit='ns'
            )

        assert 'rmsd_data' in captured_data
        assert captured_data['time_corrected'] is True
        assert captured_data['time_unit'] == 'ns'
        assert captured_data['time_interval_between_frames'] == 2.0
        assert captured_data['time_length'] == pytest.approx(0.018)

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_rmsd_time_correction_ns(self, mock_rms, mock_mda, mock_transform, tmp_path):
        """Time column should be corrected to nanoseconds when time_unit='ns'."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=100)
        mock_universe.select_atoms.return_value = mock_atoms

        frames = np.arange(5, dtype=float)
        rmsd_array = np.column_stack([
            frames,
            np.zeros(5),  # placeholder time
            np.ones(5)    # placeholder RMSD
        ])
        mock_rmsd_result = MagicMock()
        mock_rmsd_result.results.rmsd = rmsd_array
        mock_rms.RMSD.return_value = mock_rmsd_result
        mock_rmsd_result.run.return_value = mock_rmsd_result

        with patch('run_rms_analysis.pickle'):
            run_rms_analysis(
                systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
                analysis='RMSD', target_selection='sel', ref_selection='sel',
                start_frame=0, traj_format='dcd',
                time_interval_between_frames=2.0, time_unit='ns'
            )

        # frame * 2.0 ps / 1000 = ns
        expected_times = frames * 2.0 / 1000.0
        np.testing.assert_array_almost_equal(rmsd_array[:, 1], expected_times)
        assert mock_rmsd_result.results.rmsd is rmsd_array

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_rmsd_time_correction_ps(self, mock_rms, mock_mda, mock_transform, tmp_path):
        """Time column should remain in picoseconds when time_unit='ps'."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=100)
        mock_universe.select_atoms.return_value = mock_atoms

        frames = np.arange(5, dtype=float)
        rmsd_array = np.column_stack([
            frames, np.zeros(5), np.ones(5)
        ])
        mock_rmsd_result = MagicMock()
        mock_rmsd_result.results.rmsd = rmsd_array
        mock_rms.RMSD.return_value = mock_rmsd_result
        mock_rmsd_result.run.return_value = mock_rmsd_result

        with patch('run_rms_analysis.pickle'):
            run_rms_analysis(
                systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
                analysis='RMSD', target_selection='sel', ref_selection='sel',
                start_frame=0, traj_format='dcd',
                time_interval_between_frames=2.0, time_unit='ps'
            )

        expected_times = frames * 2.0
        np.testing.assert_array_almost_equal(rmsd_array[:, 1], expected_times)

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_rmsd_time_correction_us(self, mock_rms, mock_mda, mock_transform, tmp_path):
        """Time column should be corrected to microseconds when time_unit='us'."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=100)
        mock_universe.select_atoms.return_value = mock_atoms

        frames = np.arange(5, dtype=float)
        rmsd_array = np.column_stack([
            frames, np.zeros(5), np.ones(5)
        ])
        mock_rmsd_result = MagicMock()
        mock_rmsd_result.results.rmsd = rmsd_array
        mock_rms.RMSD.return_value = mock_rmsd_result
        mock_rmsd_result.run.return_value = mock_rmsd_result

        with patch('run_rms_analysis.pickle'):
            run_rms_analysis(
                systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
                analysis='RMSD', target_selection='sel', ref_selection='sel',
                start_frame=0, traj_format='dcd',
                time_interval_between_frames=2.0, time_unit='us'
            )

        # frame * 2.0 ps / 1e6 = µs
        expected_times = frames * 2.0 / 1e6
        np.testing.assert_array_almost_equal(rmsd_array[:, 1], expected_times)

    def test_rmsd_invalid_time_unit_raises(self, tmp_path):
        """An invalid time_unit should raise ValueError immediately."""
        with pytest.raises(ValueError, match="Unsupported time unit"):
            run_rms_analysis(
                systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
                analysis='RMSD', target_selection='sel', ref_selection='sel',
                start_frame=0, traj_format='dcd',
                time_unit='hours',
            )


class TestRunRmsAnalysisRMSF:
    """Tests for RMSF analysis with chain split."""

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.align_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_rmsf_chain_split_produces_chain_pickles(self, mock_rms, mock_mda,
                                                       mock_align, mock_transform, tmp_path):
        """RMSF with chain_intervals should produce per-chain pickle files."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()

        # Create fake aligned trajectory file so it uses the pre-aligned branch
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe

        # --- Build a mock AtomGroup that faithfully reproduces real MDAnalysis
        # AtomGroup iteration behavior. In production, the code does:
        #   atom_resids = np.array([atom.resid for atom in rmsf_atoms])
        # This iterates over the AtomGroup atom-by-atom. Real backbone
        # selections yield ~4 atoms per residue (N, CA, C, O). We simulate
        # that here: 4 atoms per residue for 100 residues = 400 atoms total.
        atoms_per_residue = 4
        num_residues = 100
        num_atoms = atoms_per_residue * num_residues

        # Create mock Atom objects — one per atom, each with the correct .resid
        mock_atom_list = []
        for resid in range(1, num_residues + 1):
            for _ in range(atoms_per_residue):
                atom = MagicMock()
                atom.resid = resid
                mock_atom_list.append(atom)

        # Build the mock AtomGroup: iterable over individual atoms
        mock_atoms = MagicMock()
        mock_atoms.__iter__ = MagicMock(return_value=iter(mock_atom_list))
        mock_residues = MagicMock()
        mock_residues.resids = np.arange(1, num_residues + 1)
        mock_atoms.residues = mock_residues
        mock_universe.select_atoms.return_value = mock_atoms

        # Mock RMSF results — one value per atom (not per residue).
        # The production code averages atom-level RMSF to per-residue values.
        rng = np.random.default_rng(42)
        mock_rmsf = MagicMock()
        mock_rmsf.results.rmsf = rng.uniform(0.5, 3.0, num_atoms)
        mock_rms.RMSF.return_value.run.return_value = mock_rmsf
        mock_rms.RMSF.return_value = mock_rmsf
        mock_rmsf.run.return_value = mock_rmsf

        chain_intervals = {"sys1": {"A": [1, 50], "B": [51, 100]}}

        run_rms_analysis(
            systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
            analysis='RMSF', target_selection='protein and backbone',
            ref_selection='protein and backbone', start_frame=0, traj_format='dcd',
            chain_intervals=chain_intervals
        )

        # Check that chain pickle files exist
        chain_a = tmp_path / 'rmsf_plot_sys1_wild_rep1_chainA.pkl'
        chain_b = tmp_path / 'rmsf_plot_sys1_wild_rep1_chainB.pkl'
        assert chain_a.exists(), "Chain A pickle should exist"
        assert chain_b.exists(), "Chain B pickle should exist"

        # Verify chain A content
        with open(chain_a, 'rb') as f:
            data_a = pickle.load(f)
        assert data_a['chain_id'] == 'A'
        assert len(data_a['rmsf']) == 50
        assert data_a['resids'][0] == 1  # renumbered
        assert data_a['resids'][-1] == 50

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.align_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_rmsf_no_chain_split_produces_single_pickle(self, mock_rms, mock_mda,
                                                          mock_align, mock_transform, tmp_path):
        """RMSF without chain_intervals should produce a per-residue dict pickle."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        # Build realistic atom-level mock (4 atoms per residue, 25 residues)
        atoms_per_residue = 4
        num_residues = 25
        num_atoms = atoms_per_residue * num_residues

        mock_atom_list = []
        for resid in range(1, num_residues + 1):
            for _ in range(atoms_per_residue):
                atom = MagicMock()
                atom.resid = resid
                mock_atom_list.append(atom)

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe
        mock_atoms = MagicMock()
        mock_atoms.__iter__ = MagicMock(return_value=iter(mock_atom_list))
        mock_universe.select_atoms.return_value = mock_atoms

        rng = np.random.default_rng(99)
        mock_rmsf = MagicMock()
        mock_rmsf.results.rmsf = rng.uniform(0.5, 3.0, num_atoms)
        mock_rms.RMSF.return_value = mock_rmsf
        mock_rmsf.run.return_value = mock_rmsf

        run_rms_analysis(
            systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
            analysis='RMSF', target_selection='protein and backbone',
            ref_selection='protein and backbone', start_frame=0, traj_format='dcd',
            chain_intervals=None
        )

        pkl_path = tmp_path / 'rmsf_plot_sys1_wild_rep1.pkl'
        assert pkl_path.exists(), "Single RMSF pickle should exist"

        with open(pkl_path, 'rb') as f:
            data = pickle.load(f)

        # Should be a dict with per-residue RMSF (25 residues, not 100 atoms)
        assert isinstance(data, dict)
        assert len(data['rmsf']) == num_residues
        assert len(data['resids']) == num_residues
        assert data['chain_id'] == ''


class TestPerSystemChainIntervals:
    """Tests for per-system chain_intervals with different systems."""

    def _setup_mock_system(self, tmp_path, system_name, variation, num_residues, atoms_per_residue=4):
        """Helper: create dirs and return mock atoms for a system."""
        sys_dir = tmp_path / system_name / variation
        sys_dir.mkdir(parents=True)
        (sys_dir / f'{system_name}_system_{variation}.top').touch()
        (sys_dir / f'{system_name}_production_{variation}_rep_1.dcd').touch()
        (tmp_path / f'rmsfit_{system_name}_production_{variation}_reduced_rep1.dcd').touch()

        num_atoms = atoms_per_residue * num_residues
        mock_atom_list = []
        for resid in range(1, num_residues + 1):
            for _ in range(atoms_per_residue):
                atom = MagicMock()
                atom.resid = resid
                mock_atom_list.append(atom)

        return mock_atom_list, num_atoms

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.align_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_per_system_different_intervals(self, mock_rms, mock_mda,
                                             mock_align, mock_transform, tmp_path):
        """Per-system chain_intervals should apply different ranges to different systems."""
        os.chdir(tmp_path)

        # System 1: 50 residues → chain A [1,50]
        atoms_sys1, n_atoms_sys1 = self._setup_mock_system(tmp_path, 'sys1', 'wild', 50)
        # System 2: 80 residues → chain A [1,80]
        atoms_sys2, n_atoms_sys2 = self._setup_mock_system(tmp_path, 'sys2', 'wild', 80)

        rng = np.random.default_rng(42)

        # Track which system is being loaded to return correct mock
        call_count = [0]

        def universe_side_effect(*args, **kwargs):
            mock_u = MagicMock()
            call_count[0] += 1
            # Odd calls → sys1 (50 residues), Even calls → sys2 (80 residues)
            if call_count[0] <= 1:
                atoms = atoms_sys1
                num_res = 50
                n_atoms = n_atoms_sys1
            else:
                atoms = atoms_sys2
                num_res = 80
                n_atoms = n_atoms_sys2

            mock_atoms = MagicMock()
            mock_atoms.__iter__ = MagicMock(return_value=iter(list(atoms)))
            mock_residues = MagicMock()
            mock_residues.resids = np.arange(1, num_res + 1)
            mock_atoms.residues = mock_residues
            mock_u.select_atoms.return_value = mock_atoms

            mock_rmsf_obj = MagicMock()
            mock_rmsf_obj.results.rmsf = rng.uniform(0.5, 3.0, n_atoms)
            mock_rmsf_obj.run.return_value = mock_rmsf_obj
            mock_rms.RMSF.return_value = mock_rmsf_obj

            return mock_u

        mock_mda.Universe.side_effect = universe_side_effect

        chain_intervals = {
            "sys1": {"A": [1, 50]},
            "sys2": {"A": [1, 80]}
        }

        run_rms_analysis(
            systems=['sys1', 'sys2'],
            variations={'sys1': ['wild'], 'sys2': ['wild']},
            num_replicates=1,
            analysis='RMSF',
            target_selection='protein and backbone',
            ref_selection='protein and backbone',
            start_frame=0,
            traj_format='dcd',
            chain_intervals=chain_intervals
        )

        # Verify per-system pickle files
        chain_sys1 = tmp_path / 'rmsf_plot_sys1_wild_rep1_chainA.pkl'
        chain_sys2 = tmp_path / 'rmsf_plot_sys2_wild_rep1_chainA.pkl'
        assert chain_sys1.exists(), "sys1 chain A pickle should exist"
        assert chain_sys2.exists(), "sys2 chain A pickle should exist"

        with open(chain_sys1, 'rb') as f:
            data_sys1 = pickle.load(f)
        with open(chain_sys2, 'rb') as f:
            data_sys2 = pickle.load(f)

        assert len(data_sys1['rmsf']) == 50, "sys1 should have 50 residues"
        assert len(data_sys2['rmsf']) == 80, "sys2 should have 80 residues"

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.align_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_legacy_flat_format_still_works(self, mock_rms, mock_mda,
                                             mock_align, mock_transform, tmp_path):
        """Legacy flat chain_intervals dict should still work via fallback."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        atoms_per_residue = 4
        num_residues = 100
        num_atoms = atoms_per_residue * num_residues

        mock_atom_list = []
        for resid in range(1, num_residues + 1):
            for _ in range(atoms_per_residue):
                atom = MagicMock()
                atom.resid = resid
                mock_atom_list.append(atom)

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe
        mock_atoms = MagicMock()
        mock_atoms.__iter__ = MagicMock(return_value=iter(mock_atom_list))
        mock_residues = MagicMock()
        mock_residues.resids = np.arange(1, num_residues + 1)
        mock_atoms.residues = mock_residues
        mock_universe.select_atoms.return_value = mock_atoms

        rng = np.random.default_rng(42)
        mock_rmsf = MagicMock()
        mock_rmsf.results.rmsf = rng.uniform(0.5, 3.0, num_atoms)
        mock_rms.RMSF.return_value = mock_rmsf
        mock_rmsf.run.return_value = mock_rmsf

        # Legacy flat format — not wrapped in system name
        chain_intervals = {"A": [1, 50], "B": [51, 100]}

        run_rms_analysis(
            systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
            analysis='RMSF', target_selection='protein and backbone',
            ref_selection='protein and backbone', start_frame=0, traj_format='dcd',
            chain_intervals=chain_intervals
        )

        chain_a = tmp_path / 'rmsf_plot_sys1_wild_rep1_chainA.pkl'
        chain_b = tmp_path / 'rmsf_plot_sys1_wild_rep1_chainB.pkl'
        assert chain_a.exists(), "Chain A pickle should exist with legacy format"
        assert chain_b.exists(), "Chain B pickle should exist with legacy format"


class TestSkipIfExists:
    """Tests for skip-if-exists behavior."""

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_skip_rmsd_if_pickle_exists(self, mock_rms, mock_mda, mock_transform, tmp_path, capsys):
        """Should skip RMSD analysis if pickle file already exists."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        # Pre-create the pickle file
        (tmp_path / 'rmsd_plot_sys1_wild_rep1.pkl').touch()

        run_rms_analysis(
            systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
            analysis='RMSD', target_selection='sel', ref_selection='sel',
            start_frame=0, traj_format='dcd'
        )

        captured = capsys.readouterr()
        assert "Skipping" in captured.out
        # RMSD no longer skips Universe creation (it checks per-selection inside)
        # but rms.RMSD should NOT have been called
        mock_rms.RMSD.assert_not_called()


class TestGroupSelections:
    """Tests for group_selections: separate RMSD runs per selection."""

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_group_selections_produce_per_selection_pickles(self, mock_rms, mock_mda, mock_transform, tmp_path):
        """Each group_selection should produce a _selN.pkl file."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe
        # select_atoms returns non-empty AtomGroup for valid selections
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=100)
        mock_universe.select_atoms.return_value = mock_atoms

        # Real MDAnalysis RMSD with groupselections produces 4 columns:
        # [frame, time, alignment_rmsd, group_rmsd]
        mock_rmsd_result = MagicMock()
        mock_rmsd_result.results.rmsd = np.column_stack([
            np.arange(10, dtype=float),
            np.arange(10, dtype=float) * 2.0,
            np.random.uniform(1.0, 3.0, 10),
            np.random.uniform(1.0, 3.0, 10)
        ])
        mock_rms.RMSD.return_value = mock_rmsd_result
        mock_rmsd_result.run.return_value = mock_rmsd_result

        # Capture pickle data — MagicMock can't be pickled, so intercept dump
        captured = {}
        def capture_dump(obj, f):
            """Write a pickleable subset and record the full dict for assertions."""
            import pickle as real_pickle
            key = os.path.basename(f.name)  # Use basename for reliable lookup
            captured[key] = dict(obj) if isinstance(obj, dict) else obj
            # Write a pickleable version (ensure arrays are plain numpy data)
            safe_obj = dict(obj) if isinstance(obj, dict) else obj
            real_pickle.dump(safe_obj, f)

        with patch('run_rms_analysis.pickle') as mock_pickle:
            mock_pickle.dump.side_effect = capture_dump
            run_rms_analysis(
                systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
                analysis='RMSD', target_selection='protein and backbone',
                ref_selection='protein and backbone',
                start_frame=0, traj_format='dcd',
                group_selections=['protein and backbone', 'nucleic'],
                time_interval_between_frames=2.0, time_unit='ns'
            )

        base = tmp_path / 'rmsd_plot_sys1_wild_rep1.pkl'
        sel0 = tmp_path / 'rmsd_plot_sys1_wild_rep1_sel0.pkl'
        sel1 = tmp_path / 'rmsd_plot_sys1_wild_rep1_sel1.pkl'
        assert base.exists(), "Canonical target-selection pickle should exist"
        assert sel0.exists(), "Selection 0 pickle should exist"
        assert not sel1.exists(), "Duplicate target selection should be deduplicated"

        # Verify selection key from captured data
        base_data = captured['rmsd_plot_sys1_wild_rep1.pkl']
        sel0_data = captured['rmsd_plot_sys1_wild_rep1_sel0.pkl']
        assert base_data['selection'] == 'protein and backbone'
        assert sel0_data['selection'] == 'nucleic'

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_invalid_selection_warns_and_skips(self, mock_rms, mock_mda, mock_transform, tmp_path, capsys):
        """An invalid selection should print a WARNING and skip, not crash."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe

        # First selection is valid, second raises error
        call_count = [0]
        def select_side_effect(sel_str):
            call_count[0] += 1
            if 'nucleic' in sel_str:
                raise Exception("Selection not found")
            mock_atoms = MagicMock()
            mock_atoms.__len__ = MagicMock(return_value=100)
            return mock_atoms
        mock_universe.select_atoms.side_effect = select_side_effect

        # 4 columns to match real MDAnalysis groupselections output
        mock_rmsd_result = MagicMock()
        mock_rmsd_result.results.rmsd = np.column_stack([
            np.arange(10, dtype=float),
            np.arange(10, dtype=float) * 2.0,
            np.random.uniform(1.0, 3.0, 10),
            np.random.uniform(1.0, 3.0, 10)
        ])
        mock_rms.RMSD.return_value = mock_rmsd_result
        mock_rmsd_result.run.return_value = mock_rmsd_result

        # Should NOT raise
        with patch('run_rms_analysis.pickle'):
            run_rms_analysis(
                systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
                analysis='RMSD', target_selection='protein',
                ref_selection='protein',
                start_frame=0, traj_format='dcd',
                group_selections=['protein and backbone', 'nucleic'],
                time_interval_between_frames=2.0, time_unit='ns'
            )

        captured = capsys.readouterr()
        assert "WARNING" in captured.out
        # sel0 file created by open() (even if pickle.dump is mocked)
        assert (tmp_path / 'rmsd_plot_sys1_wild_rep1_sel0.pkl').exists()
        # sel1 should NOT exist because the selection raised an error
        assert not (tmp_path / 'rmsd_plot_sys1_wild_rep1_sel1.pkl').exists()

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_empty_selection_warns_and_skips(self, mock_rms, mock_mda, mock_transform, tmp_path, capsys):
        """A selection matching 0 atoms should print WARNING and skip."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe

        # Return an empty AtomGroup
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=0)
        mock_universe.select_atoms.return_value = mock_atoms

        run_rms_analysis(
            systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
            analysis='RMSD', target_selection='nonexistent',
            ref_selection='nonexistent',
            start_frame=0, traj_format='dcd',
            group_selections=['nonexistent'],
            time_interval_between_frames=2.0, time_unit='ns'
        )

        captured = capsys.readouterr()
        assert "WARNING" in captured.out
        assert "0 atoms" in captured.out

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_no_group_selections_uses_target_selection(self, mock_rms, mock_mda, mock_transform, tmp_path):
        """Without group_selections, target_selection is used and no _selN suffix."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=100)
        mock_universe.select_atoms.return_value = mock_atoms

        mock_rmsd_result = MagicMock()
        mock_rmsd_result.results.rmsd = np.column_stack([
            np.arange(10, dtype=float),
            np.arange(10, dtype=float) * 2.0,
            np.random.uniform(1.0, 3.0, 10)
        ])
        mock_rms.RMSD.return_value = mock_rmsd_result
        mock_rmsd_result.run.return_value = mock_rmsd_result

        # Capture pickle data — MagicMock can't be pickled
        captured = {}
        def capture_dump(obj, f):
            import pickle as real_pickle
            key = os.path.basename(f.name)
            captured[key] = dict(obj) if isinstance(obj, dict) else obj
            safe_obj = dict(obj) if isinstance(obj, dict) else obj
            real_pickle.dump(safe_obj, f)

        with patch('run_rms_analysis.pickle') as mock_pickle:
            mock_pickle.dump.side_effect = capture_dump
            run_rms_analysis(
                systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
                analysis='RMSD', target_selection='protein and backbone',
                ref_selection='protein and backbone',
                start_frame=0, traj_format='dcd',
                group_selections=None,
                time_interval_between_frames=2.0, time_unit='ns'
            )

        # Standard naming (no _sel suffix) should be used
        std_pkl = tmp_path / 'rmsd_plot_sys1_wild_rep1.pkl'
        assert std_pkl.exists(), "Standard (non-sel) pickle should exist"

        std_data = captured['rmsd_plot_sys1_wild_rep1.pkl']
        assert std_data['selection'] == 'protein and backbone'


class TestRunRmsAnalysis2DRMSD:
    """Tests for 2D-RMSD output schema."""

    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.build_complex_selection', return_value=(MagicMock(), MagicMock(), MagicMock()))
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.diffusionmap')
    def test_2d_rmsd_pickle_contains_portable_matrix(self, mock_diffusionmap, mock_mda,
                                                     _mock_build, _mock_transform, tmp_path):
        """2D-RMSD pickle should contain dict with `dist_matrix` key."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe

        matrix_obj = MagicMock()
        matrix_obj.results.dist_matrix = np.array([[0.0, 1.0], [1.0, 0.0]])
        matrix_obj.run.return_value = matrix_obj
        mock_diffusionmap.DistanceMatrix.return_value = matrix_obj

        run_rms_analysis(
            systems=['sys1'],
            variations={'sys1': ['wild']},
            num_replicates=1,
            analysis='2D-RMSD',
            target_selection='protein and backbone',
            ref_selection='protein and backbone',
            start_frame=0,
            traj_format='dcd',
            time_interval_between_frames=2.0,
            time_unit='ns'
        )

        pkl = tmp_path / '2d_rmsd_plot_sys1_wild_rep1.pkl'
        assert pkl.exists()
        with open(pkl, 'rb') as f:
            data = pickle.load(f)
        assert isinstance(data, dict)
        assert 'dist_matrix' in data
        np.testing.assert_array_equal(data['dist_matrix'], np.array([[0.0, 1.0], [1.0, 0.0]]))


class TestRmsdCreatesAlignedTrajectory:
    """Tests that RMSD analysis creates an aligned trajectory for reuse by RMSF/2D-RMSD."""

    @patch('run_rms_analysis.align_trajectory')
    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_rmsd_creates_aligned_traj_when_missing(self, mock_rms, mock_mda,
                                                      mock_transform, mock_align, tmp_path):
        """RMSD should call align_trajectory when aligned trajectory file doesn't exist."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        # Do NOT create rmsfit_sys1_production_wild_reduced_rep1.dcd

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=100)
        mock_universe.select_atoms.return_value = mock_atoms

        mock_rmsd_result = MagicMock()
        mock_rmsd_result.results.rmsd = np.column_stack([
            np.arange(5, dtype=float),
            np.arange(5, dtype=float),
            np.ones(5)
        ])
        mock_rms.RMSD.return_value = mock_rmsd_result
        mock_rmsd_result.run.return_value = mock_rmsd_result

        with patch('run_rms_analysis.pickle'):
            run_rms_analysis(
                systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
                analysis='RMSD', target_selection='protein', ref_selection='protein',
                start_frame=0, traj_format='dcd'
            )

        # align_trajectory should have been called with 2D-RMSD strategy
        mock_align.assert_called_once_with(
            mock_universe, 'protein', '2D-RMSD', 'sys1', 'wild', 1, 'dcd', 0
        )

    @patch('run_rms_analysis.align_trajectory')
    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_rmsd_skips_aligned_traj_when_exists(self, mock_rms, mock_mda,
                                                    mock_transform, mock_align, tmp_path):
        """RMSD should NOT call align_trajectory when aligned trajectory already exists."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        # Pre-create the aligned trajectory
        (tmp_path / 'rmsfit_sys1_production_wild_reduced_rep1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=100)
        mock_universe.select_atoms.return_value = mock_atoms

        mock_rmsd_result = MagicMock()
        mock_rmsd_result.results.rmsd = np.column_stack([
            np.arange(5, dtype=float),
            np.arange(5, dtype=float),
            np.ones(5)
        ])
        mock_rms.RMSD.return_value = mock_rmsd_result
        mock_rmsd_result.run.return_value = mock_rmsd_result

        with patch('run_rms_analysis.pickle'):
            run_rms_analysis(
                systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
                analysis='RMSD', target_selection='protein', ref_selection='protein',
                start_frame=0, traj_format='dcd'
            )

        # align_trajectory should NOT have been called
        mock_align.assert_not_called()

    @patch('run_rms_analysis.align_trajectory')
    @patch('run_rms_analysis.transform_trajectory')
    @patch('run_rms_analysis.mda')
    @patch('run_rms_analysis.rms')
    def test_rmsd_align_uses_ref_selection(self, mock_rms, mock_mda,
                                             mock_transform, mock_align, tmp_path):
        """align_trajectory should be called with ref_selection, not target_selection."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=100)
        mock_universe.select_atoms.return_value = mock_atoms

        mock_rmsd_result = MagicMock()
        mock_rmsd_result.results.rmsd = np.column_stack([
            np.arange(5, dtype=float),
            np.arange(5, dtype=float),
            np.ones(5),
            np.ones(5)
        ])
        mock_rms.RMSD.return_value = mock_rmsd_result
        mock_rmsd_result.run.return_value = mock_rmsd_result

        with patch('run_rms_analysis.pickle'):
            run_rms_analysis(
                systems=['sys1'], variations={'sys1': ['wild']}, num_replicates=1,
                analysis='RMSD', target_selection='nucleic',
                ref_selection='protein and backbone',
                start_frame=10, traj_format='dcd'
            )

        # The ref_selection ('protein and backbone') should be used for alignment, not target
        mock_align.assert_called_once()
        call_args = mock_align.call_args
        assert call_args[0][1] == 'protein and backbone'  # selection argument
        assert call_args[0][2] == '2D-RMSD'  # analysis type for alignment

"""
Tests for run_hbonds_analysis.py — Hydrogen Bonds analysis.

Covers:
  - ValueError when no selections provided
  - Skip-if-exists behavior
  - between_pairs vs atom-focused selection paths
"""
import os
import sys
import numpy as np
import pytest
from unittest.mock import MagicMock, patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from run_hbonds_analysis import (
    run_hbonds_analysis,
    _build_hbonds_payload,
    _build_atom_labels_map,
    _resolve_hbond_selection_config,
)


class _DummyHBonds:
    def __init__(self, by_time=None):
        self.times = np.array([0.0, 1.0])
        self._by_time = np.array([2, 1]) if by_time is None else by_time

    def count_by_time(self):
        return self._by_time

    def count_by_type(self):
        return np.array([['DON', 'ACC', '2']], dtype=object)

    def count_by_ids(self):
        return np.array([[1, 2, 3, 2]])


class _NonPortable:
    pass


class TestPortablePayload:
    def test_build_payload_is_portable(self):
        hb = _DummyHBonds()
        payload = _build_hbonds_payload(
            hb,
            d_a_cutoff=3.5,
            d_h_a_angle_cutoff=150.0,
            start_frame=0,
            update_selections=True,
            acceptors_sel='protein and name O*',
            hydrogens_sel='protein and name H*',
            between_pairs=None,
            parallel_backend='serial',
            n_workers=1,
        )

        assert isinstance(payload['times'], np.ndarray)
        assert isinstance(payload['count_by_time'], np.ndarray)
        assert isinstance(payload['count_by_type'], np.ndarray)
        assert isinstance(payload['count_by_ids'], np.ndarray)
        assert payload['time_unit'] == 'ps'
        assert payload['time_corrected'] is False
        assert payload['hbonds_preset'] == 'custom'
        assert payload['parallel_backend'] == 'serial'
        assert payload['n_workers'] == 1

    def test_build_payload_rejects_nonportable_objects(self):
        bad = np.array([_NonPortable()], dtype=object)
        hb = _DummyHBonds(by_time=bad)

        with pytest.raises(TypeError, match='non-portable object'):
            _build_hbonds_payload(
                hb,
                d_a_cutoff=3.5,
                d_h_a_angle_cutoff=150.0,
                start_frame=0,
                update_selections=True,
                acceptors_sel=None,
                hydrogens_sel=None,
                between_pairs=[["protein", "nucleic"]],
                parallel_backend='multiprocessing',
                n_workers=4,
            )


class TestHBondPresets:
    def test_explicit_between_pairs_overrides_preset(self):
        explicit_pairs = [['protein', 'resname LIG']]
        _preset, _acc, _hyd, pairs = _resolve_hbond_selection_config(
            'protein_nucleic',
            acceptors_sel=None,
            hydrogens_sel=None,
            between_pairs=explicit_pairs,
        )
        assert pairs == explicit_pairs

    def test_invalid_preset_raises(self):
        with pytest.raises(ValueError, match='Unknown hbonds_preset'):
            _resolve_hbond_selection_config(
                'unknown_preset',
                acceptors_sel=None,
                hydrogens_sel=None,
                between_pairs=None,
            )


class TestHBondAtomLabelMap:
    def test_atom_labels_use_residue_only_labels(self):
        class _Atom:
            def __init__(self, index, atom_id, resname, resid, name, chain_id='', segid=''):
                self.index = index
                self.id = atom_id
                self.resname = resname
                self.resid = resid
                self.name = name
                self.chainID = chain_id
                self.segid = segid

        class _Universe:
            def __init__(self, atoms):
                self.atoms = atoms

        universe = _Universe([
            _Atom(10, 10, 'ASP', 10, 'OD1', chain_id='PROA'),
            _Atom(11, 11, 'ASP', 10, 'HD1', chain_id='PROA'),
            _Atom(20, 20, 'GLU', 20, 'OE2', chain_id='PROB'),
        ])
        count_by_ids = np.array([[10, 11, 20, 7]])

        labels = _build_atom_labels_map(universe, count_by_ids)

        assert labels[10] == 'ASP10'
        assert labels[11] == 'ASP10'
        assert labels[20] == 'GLU20'

    def test_atom_labels_use_residue_only_labels_without_chainid(self):
        class _Atom:
            def __init__(self, index, atom_id, resname, resid, name, segid=''):
                self.index = index
                self.id = atom_id
                self.resname = resname
                self.resid = resid
                self.name = name
                self.segid = segid

        class _Universe:
            def __init__(self, atoms):
                self.atoms = atoms

        universe = _Universe([
            _Atom(30, 30, 'SER', 30, 'OG', segid='PROC'),
        ])
        count_by_ids = np.array([[30, 30, 30, 1]])

        labels = _build_atom_labels_map(universe, count_by_ids)

        assert labels[30] == 'SER30'

    def test_atom_labels_use_residue_only_labels_without_name_suffix(self):
        class _Atom:
            def __init__(self, index, atom_id, resname, resid, name):
                self.index = index
                self.id = atom_id
                self.resname = resname
                self.resid = resid
                self.name = name
                self.chainID = ''
                self.segid = ''

        class _Universe:
            def __init__(self, atoms):
                self.atoms = atoms

        universe = _Universe([_Atom(30, 30, 'SER', 30, 'OG')])
        count_by_ids = np.array([[30, 30, 30, 1]])

        labels = _build_atom_labels_map(universe, count_by_ids)

        assert labels[30] == 'SER30'


class TestRunHBondsAnalysis:
    """Tests for the run_hbonds_analysis function."""

    def test_raises_without_selections(self, tmp_path):
        """Should raise ValueError when no selection parameters are provided."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()

        with patch('run_hbonds_analysis.mda') as mock_mda, \
             patch('run_hbonds_analysis.build_complex_selection', return_value=(MagicMock(), MagicMock(), MagicMock())), \
             patch('run_hbonds_analysis.transform_trajectory'):
            mock_mda.Universe.return_value = MagicMock()

            with pytest.raises(ValueError, match="must provide"):
                run_hbonds_analysis(
                    systems=['sys1'],
                    variations={'sys1': ['wild']},
                    num_replicates=1,
                    d_a_cutoff=3.5,
                    d_h_a_angle_cutoff=150.0,
                    start_frame=0,
                    traj_format='dcd',
                    acceptors_sel=None,
                    hydrogens_sel=None,
                    between_pairs=None
                )

    @patch('run_hbonds_analysis.mda')
    @patch('run_hbonds_analysis.hydrogenbonds')
    def test_skip_if_exists(self, mock_hb_module, mock_mda, tmp_path, capsys):
        """Should skip H-bond analysis if pickle already exists."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()
        (tmp_path / 'hbonds_plot_sys1_wild_rep1.pkl').touch()

        run_hbonds_analysis(
            systems=['sys1'],
            variations={'sys1': ['wild']},
            num_replicates=1,
            d_a_cutoff=3.5,
            d_h_a_angle_cutoff=150.0,
            start_frame=0,
            traj_format='dcd',
            acceptors_sel='protein and name O*',
            hydrogens_sel='nucleic and name H*'
        )

        captured = capsys.readouterr()
        assert "Skipping" in captured.out
        mock_mda.Universe.assert_not_called()

    @patch('run_hbonds_analysis.transform_trajectory')
    @patch('run_hbonds_analysis.build_complex_selection', return_value=(MagicMock(), MagicMock(), MagicMock()))
    @patch('run_hbonds_analysis.mda')
    @patch('run_hbonds_analysis.hydrogenbonds')
    def test_between_pairs_path(self, mock_hb_module, mock_mda, mock_build, mock_transform, tmp_path):
        """Should use between_pairs when provided."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe

        mock_hba = MagicMock()
        mock_hba.times = np.array([0.0, 1.0])
        mock_hba.count_by_time.return_value = np.array([2, 1])
        mock_hba.count_by_type.return_value = np.array([['A', 'B', '3']])
        mock_hba.count_by_ids.return_value = np.array([[1, 2, 3, 3]])
        mock_hb_module.HydrogenBondAnalysis.return_value = mock_hba

        between = [['protein and resnum 73', 'nucleic and name NH*']]

        with patch('run_hbonds_analysis.pickle') as mock_pickle:
            run_hbonds_analysis(
                systems=['sys1'],
                variations={'sys1': ['wild']},
                num_replicates=1,
                d_a_cutoff=3.5,
                d_h_a_angle_cutoff=150.0,
                start_frame=0,
                traj_format='dcd',
                between_pairs=between
            )

        mock_hb_module.HydrogenBondAnalysis.assert_called_once()
        dumped = mock_pickle.dump.call_args[0][0]
        assert isinstance(dumped, dict)
        assert {'times', 'count_by_time', 'count_by_type', 'count_by_ids'} <= set(dumped)
        assert dumped['hbonds_preset'] == 'custom'

    @patch('run_hbonds_analysis.transform_trajectory')
    @patch('run_hbonds_analysis.build_complex_selection', return_value=(MagicMock(), MagicMock(), MagicMock()))
    @patch('run_hbonds_analysis.mda')
    @patch('run_hbonds_analysis.hydrogenbonds')
    def test_atom_focused_path(self, mock_hb_module, mock_mda, mock_build, mock_transform, tmp_path):
        """Should use acceptors_sel/hydrogens_sel when provided."""
        os.chdir(tmp_path)
        sys_dir = tmp_path / 'sys1' / 'wild'
        sys_dir.mkdir(parents=True)
        (sys_dir / 'sys1_system_wild.top').touch()
        (sys_dir / 'sys1_production_wild_rep_1.dcd').touch()

        mock_universe = MagicMock()
        mock_mda.Universe.return_value = mock_universe

        mock_hba = MagicMock()
        mock_hba.times = np.array([0.0, 1.0])
        mock_hba.count_by_time.return_value = np.array([1, 0])
        mock_hba.count_by_type.return_value = np.array([['A', 'B', '1']])
        mock_hba.count_by_ids.return_value = np.array([[1, 2, 3, 1]])
        mock_hb_module.HydrogenBondAnalysis.return_value = mock_hba

        with patch('run_hbonds_analysis.pickle') as mock_pickle:
            run_hbonds_analysis(
                systems=['sys1'],
                variations={'sys1': ['wild']},
                num_replicates=1,
                d_a_cutoff=3.5,
                d_h_a_angle_cutoff=150.0,
                start_frame=0,
                traj_format='dcd',
                acceptors_sel='protein and name O*',
                hydrogens_sel='nucleic and name H*'
            )

        mock_hb_module.HydrogenBondAnalysis.assert_called_once()
        dumped = mock_pickle.dump.call_args[0][0]
        assert isinstance(dumped, dict)
        assert dumped['acceptors_sel'] == 'protein and name O*'
        assert dumped['hydrogens_sel'] == 'nucleic and name H*'
        assert dumped['hbonds_preset'] == 'custom'

"""
Shared pytest fixtures for the MD Analysis test suite.
"""
import os
import sys
import pickle
import pytest
import numpy as np
import tempfile
import shutil

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ─── Path Fixtures ────────────────────────────────────────────────────────────

@pytest.fixture
def project_root():
    """Returns the project root directory."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture
def test_data_dir(project_root):
    """Returns the .test_data directory path."""
    return os.path.join(project_root, '.test_data')


@pytest.fixture
def test_system():
    """Returns the test system name."""
    return 'Ung_G-C_4'


@pytest.fixture
def test_variation():
    """Returns the test variation name."""
    return 'wild'


@pytest.fixture
def test_top_file(test_data_dir, test_system, test_variation):
    """Returns the path to the test topology file."""
    return os.path.join(test_data_dir, test_system, test_variation,
                        f'{test_system}_system_{test_variation}.top')


@pytest.fixture
def test_traj_file(test_data_dir, test_system, test_variation):
    """Returns the path to the test trajectory file."""
    return os.path.join(test_data_dir, test_system, test_variation,
                        f'{test_system}_production_{test_variation}_rep_1.dcd')


@pytest.fixture
def test_data_available(test_top_file, test_traj_file):
    """Check if test data files exist."""
    return os.path.exists(test_top_file) and os.path.exists(test_traj_file)


@pytest.fixture
def tmp_output_dir(tmp_path):
    """Provides a temporary directory for test outputs."""
    return str(tmp_path)


@pytest.fixture
def tmp_workdir(tmp_path, test_data_dir, test_system, test_variation):
    """
    Creates a temporary working directory with symlinks to test data,
    mimicking the {system}/{variation}/ structure the scripts expect.
    """
    workdir = str(tmp_path / 'work')
    data_dir = os.path.join(workdir, test_system, test_variation)
    os.makedirs(data_dir, exist_ok=True)

    src_dir = os.path.join(test_data_dir, test_system, test_variation)
    if os.path.exists(src_dir):
        for fname in os.listdir(src_dir):
            src = os.path.join(src_dir, fname)
            dst = os.path.join(data_dir, fname)
            if os.path.isfile(src):
                os.symlink(src, dst)

    return workdir


# ─── Mock Data Classes (module-level for pickle compatibility) ────────────────

class _RMSDInnerResults:
    def __init__(self):
        self.rmsd = np.column_stack([
            np.arange(10, dtype=float),
            np.arange(10, dtype=float) * 2.0,
            np.random.uniform(1.0, 3.0, 10)
        ])


class MockRMSDResults:
    """Mock RMSD results object (pickleable)."""
    def __init__(self):
        self.results = _RMSDInnerResults()


class _RMSFInnerResults:
    def __init__(self, n=50):
        self.rmsf = np.random.uniform(0.5, 4.0, n)


class MockRMSFResults:
    """Mock RMSF results object (pickleable)."""
    def __init__(self, n=50):
        self.results = _RMSFInnerResults(n)


class _2DRMSDInnerResults:
    def __init__(self, n=20):
        self.dist_matrix = np.random.uniform(0, 5.0, (n, n))


class Mock2DRMSDResults:
    """Mock 2D-RMSD results object (pickleable)."""
    def __init__(self, n=20):
        self.results = _2DRMSDInnerResults(n)


class _HBondsInnerResults:
    def __init__(self):
        self.hbonds = np.array([
            [0, 10, 11, 20, 2.8, 155.0],
            [0, 30, 31, 40, 3.0, 160.0],
            [10, 10, 11, 20, 2.9, 152.0],
            [10, 50, 51, 60, 2.7, 165.0],
            [20, 10, 11, 20, 2.85, 158.0],
        ])


class MockHBondsResults:
    """Mock H-bonds results object (pickleable)."""
    def __init__(self):
        self.times = np.arange(0, 100, 10, dtype=float)
        self.n_frames = 10
        self.results = _HBondsInnerResults()

    def count_by_time(self):
        counts = np.zeros(len(self.times))
        for i, t in enumerate(self.times):
            counts[i] = np.sum(self.results.hbonds[:, 0] == t)
        return counts

    def count_by_type(self):
        return np.array([['TIP3:OT', 'TIP3:OT', '5']], dtype='<U21')

    def count_by_ids(self):
        return np.array([
            [10, 11, 20, 3],
            [30, 31, 40, 1],
            [50, 51, 60, 1],
        ])


# ─── Mock Data Fixtures ──────────────────────────────────────────────────────

@pytest.fixture
def mock_rmsd_pickle(tmp_path):
    """Creates a mock RMSD pickle file (new dict format)."""
    data = {
        'rmsd_data': MockRMSDResults().results.rmsd,
        'time_corrected': True,
        'time_unit': 'ns',
        'time_interval_between_frames': 2.0,
        'ref_selection': 'protein and backbone',
    }
    pkl_path = str(tmp_path / 'rmsd_plot_test_wild_rep1.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(data, f)
    return pkl_path


@pytest.fixture
def mock_rmsf_pickle(tmp_path):
    """Creates a mock RMSF pickle file (raw rms.RMSF-like object)."""
    pkl_path = str(tmp_path / 'rmsf_plot_test_wild_rep1.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(MockRMSFResults(), f)
    return pkl_path


@pytest.fixture
def mock_rmsf_chain_pickle(tmp_path):
    """Creates a mock chain-split RMSF pickle file (dict format)."""
    data = {
        'rmsf': np.random.uniform(0.5, 4.0, 25),
        'resids': np.arange(1, 26),
        'chain_id': 'A',
        'original_resids': np.arange(1, 26)
    }
    pkl_path = str(tmp_path / 'rmsf_plot_test_wild_rep1_chainA.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(data, f)
    return pkl_path


@pytest.fixture
def mock_2d_rmsd_pickle(tmp_path):
    """Creates a mock 2D-RMSD pickle file."""
    pkl_path = str(tmp_path / '2d_rmsd_plot_test_wild_rep1.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(Mock2DRMSDResults(), f)
    return pkl_path


@pytest.fixture
def mock_rog_pickle(tmp_path):
    """Creates a mock RoG pickle file (new dict format with metadata)."""
    from run_rog_analysis import RoGResults

    frames = list(range(50))
    times = [f * 0.002 for f in frames]  # ns
    rog_vals = list(np.random.uniform(15.0, 18.0, 50))

    result = RoGResults(frames, times, rog_vals)
    rog_data = {
        'rog_obj': result,
        'time_unit': 'ns',
        'selection': 'protein and backbone',
    }
    pkl_path = str(tmp_path / 'rog_plot_test_wild_rep1.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(rog_data, f)
    return pkl_path


@pytest.fixture
def mock_hbonds_pickle(tmp_path):
    """Creates a mock H-bonds pickle file using portable dict schema."""
    mock = MockHBondsResults()
    data = {
        'times': mock.times,
        'count_by_time': mock.count_by_time(),
        'count_by_type': mock.count_by_type(),
        'count_by_ids': mock.count_by_ids(),
        'd_a_cutoff': 3.5,
        'd_h_a_angle_cutoff': 150.0,
        'start_frame': 0,
        'update_selections': True,
        'acceptors_sel': 'protein and name O*',
        'hydrogens_sel': 'nucleic and name H*',
        'between_pairs': None,
        'parallel_backend': 'serial',
        'n_workers': 1,
    }
    pkl_path = str(tmp_path / 'hbonds_plot_test_wild_rep1.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(data, f)
    return pkl_path


# ─── Comparison-plot fixtures ─────────────────────────────────────────────────

@pytest.fixture
def mock_rmsd_pickle_with_selection(tmp_path):
    """Creates an RMSD pickle containing the 'selection' metadata key."""
    data = {
        'rmsd_data': MockRMSDResults().results.rmsd,
        'time_corrected': True,
        'time_unit': 'ns',
        'time_interval_between_frames': 2.0,
        'selection': 'protein and backbone',
        'ref_selection': 'protein and backbone',
    }
    pkl_path = str(tmp_path / 'rmsd_plot_test_wild_rep1_sel0.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(data, f)
    return pkl_path


@pytest.fixture
def mock_rmsd_comparison_pickles(tmp_path):
    """Creates two RMSD pickles suitable for plot_rmsd_comparison tests."""
    paths = []
    for name in ['sysA_wild_rep1', 'sysB_wild_rep1']:
        data = {
            'rmsd_data': MockRMSDResults().results.rmsd,
            'time_corrected': True,
            'time_unit': 'ns',
            'time_interval_between_frames': 2.0,
            'selection': 'protein and backbone',
            'ref_selection': 'protein and backbone',
        }
        pkl_path = str(tmp_path / f'rmsd_plot_{name}.pkl')
        with open(pkl_path, 'wb') as f:
            pickle.dump(data, f)
        paths.append(pkl_path)
    return paths


@pytest.fixture
def mock_rmsf_comparison_pickles(tmp_path):
    """Creates two RMSF chain-dict pickles suitable for plot_rmsf_comparison tests."""
    paths = []
    for name in ['sysA_wild_rep1_chainA', 'sysB_wild_rep1_chainA']:
        data = {
            'rmsf': np.random.uniform(0.5, 4.0, 25),
            'resids': np.arange(1, 26),
            'chain_id': 'A',
            'original_resids': np.arange(1, 26),
        }
        pkl_path = str(tmp_path / f'rmsf_plot_{name}.pkl')
        with open(pkl_path, 'wb') as f:
            pickle.dump(data, f)
        paths.append(pkl_path)
    return paths


@pytest.fixture
def mock_rmsd_replicate_groups(tmp_path):
    """Creates two groups of 2 replicate RMSD pickles for averaging tests."""
    groups = []
    for sys_name in ['sysA_wild', 'sysB_wild']:
        group = []
        for rep in [1, 2]:
            data = {
                'rmsd_data': MockRMSDResults().results.rmsd,
                'time_corrected': True,
                'time_unit': 'ns',
                'time_interval_between_frames': 2.0,
                'selection': 'protein and backbone',
                'ref_selection': 'protein and backbone',
            }
            pkl_path = str(tmp_path / f'rmsd_plot_{sys_name}_rep{rep}.pkl')
            with open(pkl_path, 'wb') as f:
                pickle.dump(data, f)
            group.append(pkl_path)
        groups.append(group)
    return groups


@pytest.fixture
def mock_rmsf_replicate_groups(tmp_path):
    """Creates two groups of 2 replicate RMSF pickles for averaging tests."""
    groups = []
    for sys_name in ['sysA_wild', 'sysB_wild']:
        group = []
        for rep in [1, 2]:
            data = {
                'rmsf': np.random.uniform(0.5, 4.0, 25),
                'resids': np.arange(1, 26),
                'chain_id': 'A',
                'original_resids': np.arange(1, 26),
            }
            pkl_path = str(tmp_path / f'rmsf_plot_{sys_name}_rep{rep}_chainA.pkl')
            with open(pkl_path, 'wb') as f:
                pickle.dump(data, f)
            group.append(pkl_path)
        groups.append(group)
    return groups

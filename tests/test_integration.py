"""
Integration tests that use the real .test_data/ trajectory files.

These tests load actual AMBER topology + DCD trajectory files and run
abbreviated analyses. They are slower and require test data to be present.

Mark: @pytest.mark.integration
Skip: automatically skipped when test data is not available.
"""
import os
import sys
import pickle
import pytest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

pytestmark = pytest.mark.integration


@pytest.fixture
def integration_workdir(tmp_path, test_data_dir, test_system, test_variation):
    """
    Creates a temporary working directory with symlinks to real test data,
    mimicking the {system}/{variation}/ structure the scripts expect.
    Changes into that directory for the duration of the test.
    """
    workdir = tmp_path / 'work'
    data_dir = workdir / test_system / test_variation
    data_dir.mkdir(parents=True)

    src_dir = os.path.join(test_data_dir, test_system, test_variation)
    if not os.path.exists(src_dir):
        pytest.skip(f"Test data not available at {src_dir}")

    for fname in os.listdir(src_dir):
        src = os.path.join(src_dir, fname)
        dst = str(data_dir / fname)
        if os.path.isfile(src):
            os.symlink(src, dst)

    original_dir = os.getcwd()
    os.chdir(str(workdir))
    yield str(workdir)
    os.chdir(original_dir)


class TestIntegrationRMSD:
    """Integration test for RMSD analysis with real data."""

    def test_rmsd_real_data(self, integration_workdir, test_system, test_variation, test_data_available):
        """Run actual RMSD analysis on test trajectory."""
        if not test_data_available:
            pytest.skip("Test data not available")

        from run_rms_analysis import run_rms_analysis

        run_rms_analysis(
            systems=[test_system],
            variations={test_system: [test_variation]},
            num_replicates=1,
            analysis='RMSD',
            target_selection='protein and backbone',
            ref_selection='protein and backbone',
            start_frame=0,
            traj_format='dcd',
            time_interval_between_frames=2.0,
            time_unit='ns'
        )

        pkl = os.path.join(integration_workdir, f'rmsd_plot_{test_system}_{test_variation}_rep1.pkl')
        assert os.path.exists(pkl)

        with open(pkl, 'rb') as f:
            data = pickle.load(f)

        assert isinstance(data, dict)
        assert data['time_corrected'] is True
        rmsd_data = data['rmsd_data']
        assert rmsd_data.shape[0] > 0
        assert rmsd_data.shape[1] >= 3
        # Time should be in ns (values should be < 1 for a 251-frame traj at 2ps)
        assert rmsd_data[-1, 1] < 1.0  # 250 * 2 / 1000 = 0.5 ns


class TestIntegrationRoG:
    """Integration test for RoG analysis with real data."""

    def test_rog_real_data(self, integration_workdir, test_system, test_variation, test_data_available):
        """Run actual RoG analysis on test trajectory."""
        if not test_data_available:
            pytest.skip("Test data not available")

        from run_rog_analysis import run_rog_analysis

        run_rog_analysis(
            systems=[test_system],
            variations={test_system: [test_variation]},
            num_replicates=1,
            start_frame=0,
            traj_format='dcd',
            selection='protein and backbone'
        )

        pkl = os.path.join(integration_workdir, f'rog_plot_{test_system}_{test_variation}_rep1.pkl')
        assert os.path.exists(pkl)

        with open(pkl, 'rb') as f:
            data = pickle.load(f)

        # New dict format: {'rog_obj': RoGResults, 'time_unit': str}
        assert isinstance(data, dict)
        assert 'rog_obj' in data
        assert 'time_unit' in data
        rog_obj = data['rog_obj']
        assert hasattr(rog_obj, 'rog_data')
        assert rog_obj.rog_data.shape[0] > 0
        assert rog_obj.rog_data.shape[1] == 3
        # RoG values should be positive and reasonable (10-50 Å for proteins)
        rog_vals = rog_obj.rog_data[:, 2]
        assert np.all(rog_vals > 0)
        assert np.all(rog_vals < 100)


class TestIntegrationRMSF:
    """Integration test for RMSF analysis with real data."""

    def test_rmsf_real_data(self, integration_workdir, test_system, test_variation, test_data_available):
        """Run actual RMSF analysis on test trajectory."""
        if not test_data_available:
            pytest.skip("Test data not available")

        from run_rms_analysis import run_rms_analysis

        run_rms_analysis(
            systems=[test_system],
            variations={test_system: [test_variation]},
            num_replicates=1,
            analysis='RMSF',
            target_selection='protein and backbone',
            ref_selection='protein and backbone',
            start_frame=0,
            traj_format='dcd',
            chain_intervals={"A": [1, 239]}
        )

        chain_pkl = os.path.join(integration_workdir, f'rmsf_plot_{test_system}_{test_variation}_rep1_chainA.pkl')
        assert os.path.exists(chain_pkl)

        with open(chain_pkl, 'rb') as f:
            data = pickle.load(f)

        assert data['chain_id'] == 'A'
        assert len(data['rmsf']) > 0
        assert np.all(data['rmsf'] >= 0)

    def test_rmsf_chain_split_real_data(self, integration_workdir, test_system, test_variation, test_data_available):
        """Run RMSF with chain split on test trajectory."""
        if not test_data_available:
            pytest.skip("Test data not available")

        from run_rms_analysis import run_rms_analysis

        # Protein resids 1-239
        chain_intervals = {"A": [1, 120], "B": [121, 239]}

        run_rms_analysis(
            systems=[test_system],
            variations={test_system: [test_variation]},
            num_replicates=1,
            analysis='RMSF',
            target_selection='protein and backbone',
            ref_selection='protein and backbone',
            start_frame=0,
            traj_format='dcd',
            chain_intervals=chain_intervals
        )

        chain_a = os.path.join(integration_workdir, f'rmsf_plot_{test_system}_{test_variation}_rep1_chainA.pkl')
        chain_b = os.path.join(integration_workdir, f'rmsf_plot_{test_system}_{test_variation}_rep1_chainB.pkl')
        assert os.path.exists(chain_a)
        assert os.path.exists(chain_b)

        with open(chain_a, 'rb') as f:
            data_a = pickle.load(f)
        with open(chain_b, 'rb') as f:
            data_b = pickle.load(f)

        assert data_a['chain_id'] == 'A'
        assert data_b['chain_id'] == 'B'
        # Chain A: resids 1-120, Chain B: resids 121-239
        assert len(data_a['rmsf']) > 0
        assert len(data_b['rmsf']) > 0
        assert data_a['resids'][0] == 1
        assert data_b['resids'][0] == 1  # renumbered


class TestIntegrationPlotting:
    """Integration tests that produce plots from real analysis results."""

    def test_rmsd_plot_from_real_data(self, integration_workdir, test_system, test_variation, test_data_available):
        """Run RMSD analysis then produce plot from real data."""
        if not test_data_available:
            pytest.skip("Test data not available")

        from run_rms_analysis import run_rms_analysis
        from plotting.plot_rmsd import plot_rmsd

        run_rms_analysis(
            systems=[test_system],
            variations={test_system: [test_variation]},
            num_replicates=1,
            analysis='RMSD',
            target_selection='protein and backbone',
            ref_selection='protein and backbone',
            start_frame=0,
            traj_format='dcd',
            time_interval_between_frames=2.0,
            time_unit='ns'
        )

        pkl = os.path.join(integration_workdir, f'rmsd_plot_{test_system}_{test_variation}_rep1.pkl')
        plot_dir = os.path.join(integration_workdir, 'plots')
        plot_rmsd(pkl, output_dir=plot_dir, dpi=72)

        pngs = [f for f in os.listdir(plot_dir) if f.endswith('.png')]
        assert len(pngs) == 1
        # File should be non-trivially sized (at least a few KB)
        png_path = os.path.join(plot_dir, pngs[0])
        assert os.path.getsize(png_path) > 1000

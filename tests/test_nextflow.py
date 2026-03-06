"""Nextflow pipeline tests.

These tests require ``nextflow`` to be installed and available in ``$PATH``.
They are gated behind the ``@pytest.mark.nextflow`` marker and will be
automatically skipped when Nextflow is not available.

Run with::

    uv run pytest tests/test_nextflow.py -v -m nextflow
"""
import os
import pickle
import shutil
import subprocess
import sys
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from run_rog_analysis import RoGResults  # noqa: F401, E402  — needed for pickle

# ─── Skip if Nextflow is not installed ────────────────────────────────────────
NF_AVAILABLE = shutil.which('nextflow') is not None

pytestmark = [
    pytest.mark.nextflow,
    pytest.mark.skipif(not NF_AVAILABLE, reason='Nextflow not installed'),
]


@pytest.fixture
def project_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture
def nf_outdir(tmp_path):
    """Temporary output directory for Nextflow results."""
    return str(tmp_path / 'nf_results')


# ─── Tests ────────────────────────────────────────────────────────────────────

class TestNextflowSmoke:
    """Smoke tests that run the pipeline with -profile test and verify output structure."""

    def test_nextflow_version(self):
        """Nextflow is callable and reports a version."""
        result = subprocess.run(['nextflow', '-version'],
                                capture_output=True, text=True, timeout=30)
        assert result.returncode == 0
        assert 'nextflow' in result.stdout.lower() or 'nextflow' in result.stderr.lower()

    def test_pipeline_help(self, project_root):
        """Pipeline can at least be parsed by Nextflow (syntax check)."""
        result = subprocess.run(
            ['nextflow', 'run', os.path.join(project_root, 'main.nf'), '--help'],
            capture_output=True, text=True, timeout=60,
            cwd=project_root)
        # --help is not implemented, but NF should at least parse the file.
        # Returncode may be 0 or 1 depending on NF version, but we shouldn't
        # get a syntax/compile error.
        assert 'Compilation error' not in result.stderr
        assert 'Unable to parse' not in result.stderr


class TestNextflowTestProfile:
    """Run the pipeline with -profile test (requires .test_data/)."""

    @pytest.fixture
    def test_data_exists(self, project_root):
        test_dir = os.path.join(project_root, '.test_data', 'Ung_G-C_4', 'wild')
        if not os.path.isdir(test_dir):
            pytest.skip('Test data not available at .test_data/')
        return True

    def test_rmsd_produces_pickle(self, project_root, nf_outdir, test_data_exists):
        """Run RMSD analysis via Nextflow and verify pickle output."""
        result = subprocess.run(
            ['nextflow', 'run', os.path.join(project_root, 'main.nf'),
             '-profile', 'test',
             '--run_rmsf', 'false',
             '--run_rog', 'false',
             '--plot_rmsd', 'false',
             '--plot_rmsf', 'false',
             '--plot_rog', 'false',
             '--outdir', nf_outdir],
            capture_output=True, text=True, timeout=1800,
            cwd=project_root)

        assert result.returncode == 0, f"Nextflow failed:\n{result.stderr}"

        # Check RMSD pickle was published
        rmsd_dir = os.path.join(nf_outdir, 'rmsd')
        assert os.path.isdir(rmsd_dir), f"Missing rmsd output dir: {rmsd_dir}"
        pickles = [f for f in os.listdir(rmsd_dir) if f.endswith('.pkl')]
        assert len(pickles) >= 1, f"No RMSD pickles found in {rmsd_dir}"

        # Validate portable pickle content.
        pkl_path = os.path.join(rmsd_dir, pickles[0])
        with open(pkl_path, 'rb') as f:
            data = pickle.load(f)
        assert isinstance(data, dict)
        assert 'rmsd_data' in data

    def test_rog_produces_pickle(self, project_root, nf_outdir, test_data_exists):
        """Run RoG analysis via Nextflow and verify pickle output."""
        result = subprocess.run(
            ['nextflow', 'run', os.path.join(project_root, 'main.nf'),
             '-profile', 'test',
             '--run_rmsd', 'false',
             '--run_rmsf', 'false',
             '--plot_rmsd', 'false',
             '--plot_rmsf', 'false',
             '--plot_rog', 'false',
             '--outdir', nf_outdir],
            capture_output=True, text=True, timeout=1800,
            cwd=project_root)

        assert result.returncode == 0, f"Nextflow failed:\n{result.stderr}"

        rog_dir = os.path.join(nf_outdir, 'rog')
        assert os.path.isdir(rog_dir)
        pickles = [f for f in os.listdir(rog_dir) if f.endswith('.pkl')]
        assert len(pickles) >= 1

        pkl_path = os.path.join(rog_dir, pickles[0])
        with open(pkl_path, 'rb') as f:
            data = pickle.load(f)
        assert isinstance(data, dict)
        assert 'rog_obj' in data

    def test_hbonds_plotting_path(self, project_root, nf_outdir, test_data_exists):
        """Run H-bonds + plotting via Nextflow and verify portable outputs/plots."""
        result = subprocess.run(
            ['nextflow', 'run', os.path.join(project_root, 'main.nf'),
             '-profile', 'test',
             '--run_rmsd', 'false',
             '--run_rmsf', 'false',
             '--run_rog', 'false',
             '--run_hbonds', 'true',
             '--plot_rmsd', 'false',
             '--plot_rmsf', 'false',
             '--plot_rog', 'false',
             '--plot_hbonds', 'true',
             '--between_pairs', '[["protein", "nucleic"]]',
             '--outdir', nf_outdir],
            capture_output=True, text=True, timeout=1800,
            cwd=project_root)

        assert result.returncode == 0, f"Nextflow failed:\n{result.stderr}"

        hb_dir = os.path.join(nf_outdir, 'hbonds')
        assert os.path.isdir(hb_dir)
        pickles = [f for f in os.listdir(hb_dir) if f.endswith('.pkl')]
        assert len(pickles) >= 1

        with open(os.path.join(hb_dir, pickles[0]), 'rb') as f:
            data = pickle.load(f)
        assert isinstance(data, dict)
        assert 'count_by_time' in data
        assert 'times' in data

        plots_dir = os.path.join(nf_outdir, 'plots', 'hbonds')
        assert os.path.isdir(plots_dir)
        pngs = [f for f in os.listdir(plots_dir) if f.endswith('.png')]
        assert len(pngs) >= 1


class TestNextflowResume:
    """Test that -resume works without re-running completed tasks."""

    @pytest.fixture
    def test_data_exists(self, project_root):
        test_dir = os.path.join(project_root, '.test_data', 'Ung_G-C_4', 'wild')
        if not os.path.isdir(test_dir):
            pytest.skip('Test data not available at .test_data/')
        return True

    def test_resume_uses_cache(self, project_root, nf_outdir, test_data_exists):
        """Second run with -resume should use cached results."""
        base_cmd = [
            'nextflow', 'run', os.path.join(project_root, 'main.nf'),
            '-profile', 'test',
            '--run_rmsf', 'false',
            '--run_rog', 'false',
            '--plot_rmsd', 'false',
            '--outdir', nf_outdir,
        ]

        # First run
        r1 = subprocess.run(base_cmd, capture_output=True, text=True,
                            timeout=1800, cwd=project_root)
        assert r1.returncode == 0, f"First run failed:\n{r1.stderr}"

        # Second run with -resume
        r2 = subprocess.run(base_cmd + ['-resume'], capture_output=True,
                            text=True, timeout=300, cwd=project_root)
        assert r2.returncode == 0
        # Look for "cached" in the output (case-insensitive; Nextflow
        # ≥25.x uses lowercase "cached: N" in its progress line).
        assert 'cached' in r2.stdout.lower(), \
            "Expected cached processes on resume"

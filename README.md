# Basic_Analysis_Set_MDAnalysis

A toolkit for standard MD trajectory analyses using [MDAnalysis](https://www.mdanalysis.org/).
Supports **RMSD**, **2D-RMSD**, **RMSF** (with per-chain split), **Radius of Gyration**,
and **Hydrogen Bonds**, with professional publication-quality plotting and an optional
Nextflow pipeline for batch execution.

## Features

- **RMSD** with DCD time-axis correction and KDE density sideplot
- **RMSF** with per-chain splitting and 1-based PDB-style residue validation
- **Named comparison plot groups** for RMSD/RMSF (`[plot_groups]`, `replicate_mode = separate|average`)
- **2D-RMSD** pairwise distance matrix heatmaps
- **Radius of Gyration** over time with KDE density sideplot
- **Hydrogen Bonds** — count by time, type, and residue-only H-bond labels, with replicate averaging support
- Consistent, accessible color palette and publication-quality styling
- **Native parallelization** via MDAnalysis 2.8+ split-apply-combine (multiprocessing / Dask) for RMSD and H-bonds
- Nextflow DSL2 pipeline with configurable per-analysis toggles and MDAnalysis parallelization support
- INI → Nextflow YAML converter for cross-runner parity
- Semantic output comparison tool (`compare_outputs.py`)
- Parallel vs serial comparison tool (`compare_parallel_serial.py`)
- Comprehensive pytest test suite (395 unit tests + 5 integration tests + 6 Nextflow tests)

## Installation

> **Recommended**: Use [Astral uv](https://docs.astral.sh/uv/) for
> dependency management.  It is fast, reproducible, and manages the
> virtual environment automatically.

```bash
# Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Set up the project (creates .venv, installs all deps including dev tools)
uv sync --all-extras
```

To enable the **Dask** parallel backend (optional):

```bash
uv sync --extra parallel
# or: pip install 'dask[distributed]'
```

<details>
<summary>Alternative: pip (not recommended)</summary>

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install pytest pytest-html pytest-cov   # dev deps
```

</details>

## Directory Structure

The scripts accept trajectory data organized as either of these layouts:
```
{system}/{system}_system_{variation}.top
{system}/{system}_production_{variation}_rep_{N}.dcd

# legacy nested layout still supported
{system}/{variation}/{system}_system_{variation}.top
{system}/{variation}/{system}_production_{variation}_rep_{N}.dcd
```

## Project Layout

```
├── executor.py              # INI-driven Python orchestrator
├── run_rms_analysis.py      # RMSD, RMSF, 2D-RMSD analyses
├── run_rog_analysis.py      # Radius of Gyration analysis
├── run_hbonds_analysis.py   # Hydrogen Bond analysis
├── parallelization.py       # MDAnalysis 2.8+ parallel backend helpers
├── utils.py                 # Shared trajectory transformations
├── plotting/
│   ├── style.py             # Shared colors, styling, helpers
│   ├── plot_rmsd.py         # RMSD + KDE plot
│   ├── plot_rmsf.py         # Per-residue RMSF (line/bar)
│   ├── plot_2d_rmsd.py      # 2D-RMSD heatmap
│   ├── plot_rog.py          # RoG + KDE plot
│   └── plot_hbonds.py       # H-bond plots (by time/type/IDs)
├── plot_group_comparisons.py # Grouped RMSD/RMSF comparison plotting helper
├── ini_to_nf_params.py      # INI → Nextflow YAML converter
├── compare_outputs.py       # Semantic pickle comparison (Python vs Nextflow)
├── compare_parallel_serial.py  # Parallel vs serial result comparison
├── main.nf                  # Nextflow DSL2 pipeline
├── nextflow.config          # Nextflow configuration
├── conf/                    # Nextflow profiles (base, test, local16)
├── tests/                   # pytest test suite
├── pyproject.toml           # uv / PEP 621 project metadata
├── uv.lock                  # Reproducible dependency lock
└── docs/
    ├── ANALYZING_PROTOCOL.md
    ├── MANUAL_TESTING_PYTHON.md
    ├── MANUAL_TESTING_NEXTFLOW.md
    ├── REPRODUCIBILITY_GUIDE.md
    └── TEST_DETAILS.md
```

## Quick Start

### Python orchestrator (serial)

```bash
# Full pipeline driven by INI config
uv run python executor.py example_pipeline.ini

# Individual analysis
uv run python run_rms_analysis.py \
  --systems '["Ung_G-C_4"]' \
  --variations '{"Ung_G-C_4": ["wild"]}' \
  --analysis RMSD \
  --target-selection 'protein and backbone' \
  --ref-selection 'protein and backbone' \
  --num-replicates 3 --start-frame 500 \
  --time-interval-between-frames 2.0 --time-unit ns
```

### Nextflow pipeline (parallel — recommended for heavy workloads)

```bash
# Run with test data
nextflow run main.nf -profile test

# Run with test data + local 16-core/16-GB tuning profile
nextflow run main.nf -profile test,local16

# Run with params converted from INI
uv run python ini_to_nf_params.py my_config.ini -o params.yml
nextflow run main.nf -params-file params.yml

# Grouped comparison plotting (RMSD/RMSF only) is propagated from INI
# [plot_groups] into Nextflow params via ini_to_nf_params.py

# With MDAnalysis-level parallelization (RMSD + H-bonds)
nextflow run main.nf -profile test \
  --parallel_backend multiprocessing --n_workers 4

# local16 profile already enables multiprocessing defaults
# (override with CLI flags if needed)
nextflow run main.nf -profile test,local16
```

> **Note**: 2D-RMSD and H-bonds are compute-intensive.  For production
> systems with multiple replicates, prefer the Nextflow pipeline which
> parallelises each (system, variation, replicate) independently.

## Plot Groups (RMSD/RMSF)

Named comparison plot groups let you overlay systems/variations in shared plots.
This is supported in both `executor.py` and the Nextflow pipeline.

Scope and behavior:
- Supported analyses: RMSD and RMSF only.
- Unsupported grouped comparisons: 2D-RMSD, RoG, H-bonds.
- `replicate_mode = separate` creates one comparison plot per replicate.
- `replicate_mode = average` averages replicates before plotting.
- Missing member pickles are warned and skipped (run continues).

Example INI section:

```ini
[plot_groups]
replicate_mode = # average or separate
systemX_variationX_vs_systemX_variationX = [["systemX", "variationX"], "]]
```

## Parallelization (MDAnalysis 2.8+)

RMSD and Hydrogen Bonds analyses support native MDAnalysis parallelization
via the split-apply-combine framework.  Both Python (`executor.py`) and
Nextflow orchestration modes support these settings identically.

### Python executor (INI config)

```ini
[parallelization]
backend = multiprocessing   ; or 'dask' (requires dask[distributed])
n_workers = 4               ; or 'none' for auto-detect (min(cpu_count, 4))
```

### Nextflow

```bash
nextflow run main.nf -profile test \
  --parallel_backend multiprocessing --n_workers 4
```

When a parallel backend is active, Nextflow dynamically sets `cpus`
for `RUN_RMSD` and `RUN_HBONDS` to match `n_workers`.

### Supported analyses

| Analysis | Supports parallel? |
|----------|--------------------|
| RMSD | ✅ multiprocessing / dask |
| H-bonds | ✅ multiprocessing / dask |
| RMSF | ❌ serial only |
| 2D-RMSD | ❌ serial only |
| RoG | ❌ serial only |

Parallel arguments are **not** passed to serial-only analyses in either
orchestration mode.

**Safety**: If trajectory transformations are not parallelizable, the runner
automatically falls back to serial execution with a warning.

**Validate** parallel results match serial:

```bash
# Automatic — runs serial & parallel pipelines, then compares:
uv run python compare_parallel_serial.py \
    --config test_comparison.ini \
    --backend multiprocessing --n-workers 4

# Or reuse existing serial results for speed:
uv run python compare_parallel_serial.py \
    --config test_comparison.ini \
    --backend multiprocessing --n-workers 4 \
    --dir-serial results_comparison/work
```

## Testing

```bash
# Unit tests only (fast, no data needed)
uv run pytest tests/ --ignore=tests/test_integration.py --ignore=tests/test_nextflow.py -q

# All tests (requires .test_data/)
uv run pytest tests/ -v

# Nextflow tests (requires Nextflow + .test_data/)
uv run pytest tests/ -m nextflow -v
```

## Analysis Requirements

1. **Basic Configuration**
   The analysis will at least consist of x systems, y variations, z replicates, topology format and trajectory format, and start frame (ex. 500).

   ### Hydrogen Bonds Analysis

   Add one of two pairs of atom selections:

   - **Option 1:** Atom-focused analysis
     - `acceptors_sel` (e.g. `'protein and name O*'`)
     - `hydrogens_sel` (e.g. `'nucleic and name H*'`)

   - **Option 2:** Pair-focused analysis
     - `between_pairs` list of lists (e.g. `[['protein and resnum 73', 'nucleic and name NH*'], ['protein and resnum 73', 'nucleic and name N*'], ...]`)

   **Additional parameters:**
   - `d_a_cutoff` (e.g. 3.5)
   - `d_h_a_angle_cutoff` (e.g. 150)
   - `update_selections` may be passed as `False` to improve performance, but this will keep selection from updating over time, which is generally not recommended.

   ### RMS* Analysis (RMSF, RMSD or 2D-RMSD)

   - `reference_selection` (e.g. `'protein and backbone'`)
   - `target_selection` (e.g. `''`)
   - `chain_intervals` (optional, for RMSF chain split, e.g. `'{"A": [1, 120], "B": [121, 239]}'`)
  - `time_interval_between_frames` (required when RMSD is enabled; in ps)

2. **Directory Structure**
   The script will be run on the parent folder for both systems, where each has its own folder:

   ```text
   ./{system}/{variation}/
   ```

## Configuration

Configuration is delegated to CLI arguments (JSON strings or file paths).

Old/Pure Python Version:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16389354.svg)](https://doi.org/10.5281/zenodo.16389354)

Current Version Deposition Pending, due to requiring further testing.

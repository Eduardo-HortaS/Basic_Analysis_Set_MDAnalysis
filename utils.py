import os
from MDAnalysis.analysis import align
import MDAnalysis.transformations as trans


# ─── Time-unit conversion utilities ──────────────────────────────────────────
#
# MDAnalysis stores trajectory timestamps in picoseconds (ps) internally.
# These helpers convert from ps to any supported unit and provide pretty
# labels for plot axes.
#
# Canonical supported units: fs, ps, ns, us, ms, s
# The alias 'μs' (Unicode mu) is accepted in look-ups but the canonical
# spelling is 'us'.

TIME_UNITS: dict[str, float] = {
    'fs':  1e-3,   # 1 ps = 1000 fs   → divide ps by 1e-3 to get fs
    'ps':  1.0,    # identity
    'ns':  1e3,    # 1 ns = 1000 ps
    'us':  1e6,    # 1 µs = 1e6 ps
    'μs':  1e6,    # Unicode alias
    'ms':  1e9,    # 1 ms = 1e9 ps
    's':   1e12,   # 1 s  = 1e12 ps
}
"""Mapping from unit name to the number of picoseconds in one of that unit."""

SUPPORTED_TIME_UNITS: tuple[str, ...] = ('fs', 'ps', 'ns', 'us', 'ms', 's')
"""Canonical set of accepted time-unit strings (no Unicode aliases)."""

TIME_UNIT_LABELS: dict[str, str] = {
    'fs': 'fs',
    'ps': 'ps',
    'ns': 'ns',
    'us': 'µs',
    'ms': 'ms',
    's':  's',
}
"""Pretty-print labels for axis / title display (e.g. 'us' → 'µs')."""


def validate_time_unit(unit: str) -> None:
    """
    Raise ``ValueError`` if *unit* is not a recognised time-unit string.

    Accepted values: any key in :data:`TIME_UNITS` (includes the ``μs``
    Unicode alias).
    """
    if unit not in TIME_UNITS:
        raise ValueError(
            f"Unsupported time unit '{unit}'. "
            f"Accepted values: {', '.join(SUPPORTED_TIME_UNITS)}"
        )


def convert_time_from_ps(time_ps, target_unit: str):
    """
    Convert a time value (or array of values) **from picoseconds** to
    *target_unit*.

    Parameters
    ----------
    time_ps : float or numpy.ndarray
        Time value(s) in picoseconds.
    target_unit : str
        One of the keys in :data:`TIME_UNITS`.

    Returns
    -------
    float or numpy.ndarray
        Converted value(s).

    Raises
    ------
    ValueError
        If *target_unit* is not recognised.
    """
    validate_time_unit(target_unit)
    return time_ps / TIME_UNITS[target_unit]


def time_unit_label(unit: str) -> str:
    """
    Return a human-friendly label for *unit* (e.g. ``'us'`` → ``'µs'``).

    Falls back to the raw string when no pretty-print mapping exists.
    """
    return TIME_UNIT_LABELS.get(unit, unit)


# ─── Trajectory utilities ────────────────────────────────────────────────────

# Recognised special values for the *wrap_selection* parameter.
_WRAP_AUTO = 'auto'


def build_complex_selection(universe, wrap_selection='auto'):
    """
    Build atom groups for PBC centering and wrapping.

    Returns a 3-tuple ``(center_ag, ligand_ag, wrap_ag)``:

    * *center_ag* — atoms to centre on (the "anchor", typically protein
      only).  The system is shifted so this group's centre of mass sits
      at the box centre.
    * *ligand_ag* — non-protein molecules that should remain near the
      protein (e.g. ligands, cofactors).  These are wrapped **by
      fragment** *after* centering so they land in the nearest periodic
      image.  ``None`` when there are no such molecules.
    * *wrap_ag* — remaining molecules (solvent, ions) wrapped by fragment
      into the primary image.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Loaded Universe.
    wrap_selection : str or None
        ``'auto'`` (default)
            Centre on ``'protein'``; wrap ``'not protein'`` (ligands,
            solvent, ions) by fragment.  In this mode *ligand_ag* is
            ``None`` because ligands are part of *wrap_ag* and will be
            placed in the nearest image automatically.
        ``None``
            Disable all wrapping transformations.  The function returns
            ``(None, None, None)``.
        Any other string
            Interpreted as an MDAnalysis selection for atoms to wrap as
            **solvent/ions** (i.e. *wrap_ag*).  The centering anchor is
            always ``'protein'``, and any atoms that are neither protein,
            nor matched by this selection, become *ligand_ag* — they are
            wrapped by fragment near the protein in a dedicated step.

    Returns
    -------
    center_ag : AtomGroup or None
    ligand_ag : AtomGroup or None
    wrap_ag : AtomGroup or None
    """
    if wrap_selection is None:
        return None, None, None

    if wrap_selection == _WRAP_AUTO:
        center_ag = universe.select_atoms('protein')
        wrap_ag = universe.select_atoms('not protein')
        ligand_ag = None
    else:
        center_ag = universe.select_atoms('protein')
        wrap_ag = universe.select_atoms(wrap_selection)
        # Atoms that are neither protein nor in wrap_ag are ligands
        # that should be kept near the protein.
        ligand_ag = universe.atoms - center_ag - wrap_ag
        if len(ligand_ag) == 0:
            ligand_ag = None

    return center_ag, ligand_ag, wrap_ag


def transform_trajectory(universe, center_selection, wrap_selection,
                         ligand_selection=None, strict_wrapping=False):
    """
    Unwrap, center, and wrap trajectory data.

    The pipeline is:

    1. **Unwrap** every atom so that all molecules are whole (no fragment
       split across periodic boundaries).
    2. **Centre** so that *center_selection* (protein) sits at the box
       centre.  No atoms are wrapped during this step.
    3. *(optional)* **Wrap** *ligand_selection* back into the primary
       image **by fragment**.  Because the protein is already centred,
       this places each ligand molecule in the nearest periodic image
       — i.e. right next to the protein.  Skipped when
       *ligand_selection* is ``None``.
    4. **Wrap** *wrap_selection* (solvent, ions) back into the primary
       image **by fragment**.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe whose trajectory will be transformed **in place**.
    center_selection : AtomGroup
        Atoms to centre on (typically the protein).
    wrap_selection : AtomGroup
        Atoms to wrap back into the primary cell by fragment (solvent/ions).
    ligand_selection : AtomGroup or None, optional
        Ligand atoms to wrap near the protein before wrapping solvent.
        When ``None``, no separate ligand-wrapping step is performed
        (appropriate for ``'auto'`` mode where ligands are already in
        *wrap_selection*).
    strict_wrapping : bool, optional
        If True, raise an error if fragment data is not available.
        If False (default), warn and fall back to residue-based wrapping.
    """

    # Prefer fragment-based transforms when bonds are available.
    # GRO/topologies without bond data have no fragments; in that case,
    # try to guess bonds once before falling back to residue-based wrapping.
    try:
        _ = universe.atoms.fragments
        has_fragments = True
    except Exception:
        has_fragments = False

    guessed_bonds = False
    if not has_fragments:
        try:
            universe.atoms.guess_bonds()
            _ = universe.atoms.fragments
            has_fragments = True
            guessed_bonds = True
            print("INFO: Inferred bonds from coordinates for fragment-aware wrapping.")
        except Exception:
            has_fragments = False

    if not has_fragments and strict_wrapping:
        raise RuntimeError(
            "ERROR: Topology lacks bond/fragment information required for proper PBC wrapping.\n"
            "  This typically happens when using GRO topology files without bond data.\n"
            "  Solutions:\n"
            f"  1. Use 'gmx trjconv -pbc whole' to pre-wrap the trajectory before analysis.\n"
            f"  2. Disable strict wrapping by removing the --strict-wrapping flag.\n"
            f"  3. Ensure your topology file includes explicit bond data (e.g., PSF for NAMD).\n"
        )

    transforms = []
    if has_fragments:
        transforms.append(trans.unwrap(universe.atoms))
    else:
        print(
            "NEXTFLOW_WRAPPING_WARNING: No bond/fragments data found; "
            "using residue-based wrapping and skipping unwrap."
        )

    transforms.append(trans.center_in_box(center_selection, wrap=False))

    wrap_compound = 'fragments' if has_fragments else 'residues'

    if ligand_selection is not None:
        transforms.append(trans.wrap(ligand_selection, compound=wrap_compound))

    transforms.append(trans.wrap(wrap_selection, compound=wrap_compound))

    universe.trajectory.add_transformations(*transforms)
    return universe

def align_trajectory(universe, selection, analysis, system, variation, rep, traj_format, start_frame):
    """
    Aligns trajectory data based on the given selection.
    """
    if analysis == "RMSF":
        average = align.AverageStructure(universe, universe, select=selection, ref_frame=0).run()
        ref_rmsf = average.results.universe
        _ = align.AlignTraj(universe, ref_rmsf, select=selection, filename=f"rmsfit_{system}_production_{variation}_reduced_rep{rep}.{traj_format}", in_memory=False).run(start=start_frame, stop=None, step=1)
    elif analysis == "2D-RMSD":
        _ = align.AlignTraj(universe, universe, select=selection, filename=f"rmsfit_{system}_production_{variation}_reduced_rep{rep}.{traj_format}", in_memory=False).run(start=start_frame, stop=None, step=1)


def resolve_trajectory_file(system, variation, rep, traj_format, base_dir=None):
    """Return a trajectory path supporting both ``rep_1`` and ``rep1`` naming.

    Parameters
    ----------
    system, variation : str
        System and variation labels.
    rep : int
        Replicate number.
    traj_format : str
        Trajectory file extension without leading dot.
    base_dir : str or None, optional
        Directory under which candidates are checked.  When None,
        candidates are checked relative to the current working directory.

    Returns
    -------
    tuple[str, list[str]]
        ``(resolved_path, candidates)`` where ``resolved_path`` is the first
        existing candidate (or the legacy default if none exists), and
        ``candidates`` contains both attempted paths in priority order.
    """
    candidates = [
        f'{system}/{variation}/{system}_production_{variation}_rep_{rep}.{traj_format}',
        f'{system}/{variation}/{system}_production_{variation}_rep{rep}.{traj_format}',
    ]

    for rel_path in candidates:
        path = os.path.join(base_dir, rel_path) if base_dir else rel_path
        if os.path.isfile(path):
            return path, candidates

    default_path = os.path.join(base_dir, candidates[0]) if base_dir else candidates[0]
    return default_path, candidates


def resolve_topology_file(system, variation, top_format, base_dir=None):
    """Return the expected topology path for a system/variation."""
    rel_path = f'{system}/{variation}/{system}_system_{variation}.{top_format}'
    path = os.path.join(base_dir, rel_path) if base_dir else rel_path
    return path, rel_path


def resolve_reference_pdb_file(system, variation, base_dir=None):
    """Return required reference PDB path for chain metadata."""
    rel_path = f'{system}/{variation}/{system}_system_{variation}.pdb'
    path = os.path.join(base_dir, rel_path) if base_dir else rel_path
    return path, rel_path


def cleanup_work_directory(work_dir, preserve_extensions=None, verbose=True):
    """
    Clean up intermediate files in work directory while preserving analysis outputs.
    
    Removes:
    - work/ directory (Nextflow intermediates)
    - results_* directories (old analysis outputs)
    - __pycache__/ directories
    - .pytest_cache/ directories
    - logs/ directory
    - *.log files
    - Old analysis scripts (if present)
    
    Preserves:
    - All *.pkl files (analysis pickles for plotting)
    - All *.png files (plots and figures)
    - Source code (*.py, *.nf, *.sh, *.yml, *.ini files)
    - Configuration and documentation
    
    Parameters
    ----------
    work_dir : str
        Path to the working directory to clean
    preserve_extensions : list, optional
        File extensions to always preserve. Default: ['.pkl', '.png']
    verbose : bool, default True
        If True, print cleanup summary
    
    Returns
    -------
    dict
        Statistics: {
            'files_removed': int,
            'dirs_removed': int,
            'space_freed_mb': float,
            'removed_items': list of removed paths
        }
    """
    if preserve_extensions is None:
        preserve_extensions = ['.pkl', '.png']
    
    removed_items = []
    space_freed = 0
    
    # Define directories to remove
    dirs_to_remove = [
        'work',
        '__pycache__',
        '.pytest_cache',
        'logs'
    ]
    
    # Define file patterns to remove
    patterns_to_remove = [
        'results_*',
        '*.log',
        '.coverage',
        '*.pyc'
    ]
    
    # Target directories in work_dir
    for dirname in dirs_to_remove:
        dir_path = os.path.join(work_dir, dirname)
        if os.path.isdir(dir_path):
            try:
                for root, dirs, files in os.walk(dir_path, topdown=False):
                    # Remove files
                    for file in files:
                        file_path = os.path.join(root, file)
                        try:
                            size = os.path.getsize(file_path)
                            os.remove(file_path)
                            space_freed += size
                            removed_items.append(file_path)
                        except Exception as e:
                            if verbose:
                                print(f"Warning: Could not remove {file_path}: {e}")
                    # Remove subdirectories
                    for dir_name in dirs:
                        dir_to_remove = os.path.join(root, dir_name)
                        try:
                            os.rmdir(dir_to_remove)
                            removed_items.append(dir_to_remove)
                        except Exception as e:
                            if verbose:
                                print(f"Warning: Could not remove {dir_to_remove}: {e}")
                # Remove main directory
                try:
                    os.rmdir(dir_path)
                    removed_items.append(dir_path)
                except Exception as e:
                    if verbose:
                        print(f"Warning: Could not remove {dir_path}: {e}")
            except Exception as e:
                if verbose:
                    print(f"Warning: Error processing {dir_path}: {e}")
    
    # Remove files matching patterns (at root level)
    import glob as glob_module
    for pattern in patterns_to_remove:
        pattern_path = os.path.join(work_dir, pattern)
        for file_path in glob_module.glob(pattern_path):
            if os.path.isfile(file_path):
                # Don't remove files with preserve extensions
                if any(file_path.endswith(ext) for ext in preserve_extensions):
                    continue
                try:
                    size = os.path.getsize(file_path)
                    os.remove(file_path)
                    space_freed += size
                    removed_items.append(file_path)
                except Exception as e:
                    if verbose:
                        print(f"Warning: Could not remove {file_path}: {e}")
    
    stats = {
        'files_removed': sum(1 for item in removed_items if os.path.isfile(item)),
        'dirs_removed': sum(1 for item in removed_items if os.path.isdir(item)),
        'space_freed_mb': round(space_freed / (1024 * 1024), 2),
        'removed_items': removed_items
    }
    
    if verbose:
        print(f"\n{'='*60}")
        print("Cleanup Summary")
        print(f"{'='*60}")
        print(f"Files removed:    {stats['files_removed']}")
        print(f"Directories removed: {stats['dirs_removed']}")
        print(f"Space freed:      {stats['space_freed_mb']} MB")
        if stats['files_removed'] + stats['dirs_removed'] > 0:
            print(f"✓ Cleanup successful")
        print(f"{'='*60}\n")
    
    return stats

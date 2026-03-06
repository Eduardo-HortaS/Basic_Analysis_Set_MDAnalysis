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
                         ligand_selection=None):
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
    """

    transforms = [
        trans.unwrap(universe.atoms),
        trans.center_in_box(center_selection, wrap=False),
    ]

    if ligand_selection is not None:
        transforms.append(trans.wrap(ligand_selection, compound='fragments'))

    transforms.append(trans.wrap(wrap_selection, compound='fragments'))

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

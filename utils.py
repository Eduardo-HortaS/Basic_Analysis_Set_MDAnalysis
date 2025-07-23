from MDAnalysis.analysis import align
import MDAnalysis.transformations as trans

def transform_trajectory(universe, prewrap_selection, wrap_selection):
    """
    Unwraps and wraps trajectory data based on the given selections.
    """

    transforms = [trans.unwrap(prewrap_selection),
                trans.center_in_box(prewrap_selection, wrap=True),
                trans.wrap(wrap_selection)]

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

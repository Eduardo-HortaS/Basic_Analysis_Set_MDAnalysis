#!/usr/bin/env nextflow

/*
 * Trajectory Preparation Nextflow Pipeline
 * =======================================
 * Nextflow DSL2 pipeline to run the python trajectory prep tool.
 * Unwraps, centers, wraps, strips, and renumbers MD trajectories.
 */

nextflow.enable.dsl = 2

// Default Parameters
params.topology      = null
params.trajectory    = null
params.outdir        = 'results'
params.anchor        = 'protein'
params.strip_sel     = 'not water and not (type CLA SOD POT)'
params.guess_bonds   = false
params.save_pdb      = true
params.save_psf      = true
params.save_traj     = true


process PREPARE_TRAJECTORY {
    tag "${trajectory.baseName}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path topology
    path trajectory

    output:
    path "*_stripped.pdb", emit: pdb, optional: true
    path "*_stripped.psf", emit: psf, optional: true
    path "*_stripped.dcd", emit: traj, optional: true

    script:
    def guess_bonds_flag = params.guess_bonds ? "--guess-bonds" : ""
    def out_pdb_arg = params.save_pdb ? "--out-pdb ${trajectory.baseName}_stripped.pdb" : ""
    def out_psf_arg = params.save_psf ? "--out-psf ${trajectory.baseName}_stripped.psf" : ""
    def out_traj_arg = params.save_traj ? "--out-traj ${trajectory.baseName}_stripped.dcd" : ""
    """
    python ${projectDir}/trajectory_prep.py \
        --topology ${topology} \
        --trajectory ${trajectory} \
        ${out_pdb_arg} \
        ${out_psf_arg} \
        ${out_traj_arg} \
        --anchor "${params.anchor}" \
        --strip-sel "${params.strip_sel}" \
        ${guess_bonds_flag}
    """
}

workflow {
    if (!params.topology || !params.trajectory) {
        error "ERROR: --topology and --trajectory parameters must be specified."
    }

    // Create channel for topology (value channel using first/collect)
    topology_ch = Channel.fromPath(params.topology).first()

    // Create channel for trajectories (queue channel, can contain multiple matching files)
    trajectory_ch = Channel.fromPath(params.trajectory)

    // Run preparation process
    PREPARE_TRAJECTORY(topology_ch, trajectory_ch)
}

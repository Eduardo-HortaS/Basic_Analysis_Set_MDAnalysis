#!/usr/bin/env nextflow

/*
 * MD Analysis Pipeline
 * ====================
 * A Nextflow DSL2 pipeline for running molecular dynamics trajectory analyses
 * (RMSD, RMSF, 2D-RMSD, RoG, Hydrogen Bonds) and generating publication-quality plots.
 *
 * Uses scripts from the Basic_Analysis_Set_MDAnalysis toolkit.
 *
 * nf-core compliance notes:
 *   - Uses DSL2 with modular processes
 *   - Configurable via nextflow.config with profiles
 *   - publishDir for results output
 *   - Params for all user-configurable options
 *   Deviations from full nf-core template:
 *   - No MultiQC (not applicable to MD analysis)
 *   - No nf-validation/tower integration (overkill for this specialized pipeline)
 *   - No nf-core create boilerplate (this is a domain-specific, not community, pipeline)
 */

nextflow.enable.dsl = 2

// ─── Default Parameters ──────────────────────────────────────────────────────
params.systems          = null    // JSON string or path: e.g., '["Ung_G-C_4"]'
params.variations       = null    // JSON string or path: e.g., '{"Ung_G-C_4": ["wild"]}'
params.num_replicates   = 1
params.traj_format      = 'dcd'
params.top_format       = 'top'
params.start_frame      = 0
params.input_dir        = '.'     // Base directory containing {system}/{variation}/ structure
params.outdir           = './results'

// Analysis toggles
params.run_rmsd         = true
params.run_rmsf         = true
params.run_2d_rmsd      = false
params.run_rog          = true
params.run_hbonds       = false

// Per-analysis plot toggles (independent of analysis toggles)
params.plot_rmsd        = true
params.plot_rmsf        = true
params.plot_2d_rmsd     = true
params.plot_rog         = true
params.plot_hbonds      = true

// RMSD-specific
params.target_selection = 'protein and backbone'
params.ref_selection    = 'protein and backbone'
params.group_selections = null
params.time_interval_between_frames = null
params.time_unit        = 'ns'

// RMSF-specific
params.chain_intervals  = null

// RoG-specific
params.rog_selection    = 'protein and backbone'

// H-bonds specific
params.acceptors_sel    = null
params.hydrogens_sel    = null
params.between_pairs    = null
params.d_a_cutoff       = 3.5
params.d_h_a_angle_cutoff = 150.0
params.update_selections = true

// PBC wrapping
params.wrap_selection   = 'auto'   // 'auto', 'none', or an MDAnalysis selection string

// Plotting
params.dpi              = 400
params.hbonds_top_n     = 20

// Parallelization — MDAnalysis 2.8+ native parallel backends
// 'serial' (default), 'multiprocessing', or 'dask'
params.parallel_backend = 'serial'
params.n_workers        = null    // null = auto-detect

// ─── Helper: Resolve script paths ────────────────────────────────────────────
def scriptDir = "${projectDir}"

// ─── Input Channel ───────────────────────────────────────────────────────────

workflow {
    // Parse systems and variations
    if (!params.systems || !params.variations) {
        error "ERROR: --systems and --variations must be provided."
    }

    // time_interval_between_frames is mandatory — it cannot be inferred from
    // trajectory files and must always be provided explicitly.
    if (params.run_rmsd && !params.time_interval_between_frames) {
        error "ERROR: 'time_interval_between_frames' is not set.\n" +
              "  This value CANNOT be inferred from trajectory files and must be provided explicitly (in picoseconds).\n" +
              "  Please set  --time_interval_between_frames <value>  and re-run.\n" +
              "  See docs/MANUAL_TESTING_NEXTFLOW.md for examples."
    }

    // Create a channel of [system, variation, rep] tuples
    Channel
        .fromList(generate_combinations())
        .set { analysis_inputs }

    // Run analyses conditionally
    if (params.run_rmsd) {
        RUN_RMSD(analysis_inputs)
        if (params.plot_rmsd) {
            PLOT_RMSD(RUN_RMSD.out.pickle)
        }
    }

    if (params.run_rmsf) {
        RUN_RMSF(analysis_inputs)
        if (params.plot_rmsf) {
            PLOT_RMSF(RUN_RMSF.out.pickle.flatten())
        }
    }

    if (params.run_2d_rmsd) {
        RUN_2D_RMSD(analysis_inputs)
        if (params.plot_2d_rmsd) {
            PLOT_2D_RMSD(RUN_2D_RMSD.out.pickle)
        }
    }

    if (params.run_rog) {
        RUN_ROG(analysis_inputs)
        if (params.plot_rog) {
            PLOT_ROG(RUN_ROG.out.pickle)
        }
    }

    if (params.run_hbonds) {
        RUN_HBONDS(analysis_inputs)
        if (params.plot_hbonds) {
            PLOT_HBONDS(RUN_HBONDS.out.pickle)
        }
    }
}

// ─── Helper function ─────────────────────────────────────────────────────────

import groovy.json.JsonSlurper

def generate_combinations() {
    def jsonSlurper = new JsonSlurper()
    def systems
    def variations

    // Parse systems
    def systems_file = new File(params.systems.toString())
    if (systems_file.exists()) {
        systems = jsonSlurper.parse(systems_file)
    } else {
        systems = jsonSlurper.parseText(params.systems.toString())
    }

    // Parse variations
    def variations_file = new File(params.variations.toString())
    if (variations_file.exists()) {
        variations = jsonSlurper.parse(variations_file)
    } else {
        variations = jsonSlurper.parseText(params.variations.toString())
    }

    def combos = []
    systems.each { system ->
        variations[system].each { variation ->
            (1..params.num_replicates).each { rep ->
                def traj = file("${params.input_dir}/${system}/${variation}/${system}_production_${variation}_rep_${rep}.${params.traj_format}")
                def top  = file("${params.input_dir}/${system}/${variation}/${system}_system_${variation}.${params.top_format}")
                combos.add(tuple(system, variation, rep, top, traj))
            }
        }
    }
    return combos
}

// ─── Processes ───────────────────────────────────────────────────────────────

process RUN_RMSD {
    tag "${system}_${variation}_rep${rep}"
    label 'process_medium'
    cpus { params.parallel_backend != 'serial' && params.n_workers ? params.n_workers : (params.parallel_backend != 'serial' ? 4 : 1) }
    publishDir "${params.outdir}/rmsd", mode: 'copy'

    input:
    tuple val(system), val(variation), val(rep), path(top_file), path(traj_file)

    output:
    path("rmsd_plot_${system}_${variation}_rep${rep}.pkl"), emit: pickle

    script:
    def time_arg = params.time_interval_between_frames ? "--time-interval-between-frames ${params.time_interval_between_frames} --time-unit ${params.time_unit}" : ""
    def group_arg = params.group_selections ? "--group-selections '${params.group_selections}'" : ""
    def wrap_arg = params.wrap_selection == 'none' ? "--wrap-selection none" : (params.wrap_selection != 'auto' ? "--wrap-selection '${params.wrap_selection}'" : "")
    def parallel_arg = params.parallel_backend != 'serial' ? "--parallel-backend ${params.parallel_backend}" : ""
    def nworkers_arg = params.n_workers ? "--n-workers ${params.n_workers}" : ""
    """
    mkdir -p ${system}/${variation}
    ln -sf \$(readlink -f ${top_file}) ${system}/${variation}/${system}_system_${variation}.${params.top_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep_${rep}.${params.traj_format}

    python ${scriptDir}/run_rms_analysis.py \\
        --systems '["${system}"]' \\
        --variations '{"${system}": ["${variation}"]}' \\
        --num-replicates 1 \\
        --analysis RMSD \\
        --target-selection '${params.target_selection}' \\
        --ref-selection '${params.ref_selection}' \\
        --start-frame ${params.start_frame} \\
        --traj-format ${params.traj_format} \\
        --top-format ${params.top_format} \\
        ${time_arg} ${group_arg} ${wrap_arg} ${parallel_arg} ${nworkers_arg}

    mv rmsd_plot_${system}_${variation}_rep1.pkl rmsd_plot_${system}_${variation}_rep${rep}.pkl || true
    """
}

process RUN_RMSF {
    tag "${system}_${variation}_rep${rep}"
    label 'process_medium'
    publishDir "${params.outdir}/rmsf", mode: 'copy'

    input:
    tuple val(system), val(variation), val(rep), path(top_file), path(traj_file)

    output:
    path("rmsf_plot_${system}_${variation}_rep${rep}*.pkl"), emit: pickle

    script:
    def chain_arg = params.chain_intervals ? "--chain-intervals '${params.chain_intervals}'" : ""
    def wrap_arg = params.wrap_selection == 'none' ? "--wrap-selection none" : (params.wrap_selection != 'auto' ? "--wrap-selection '${params.wrap_selection}'" : "")
    """
    mkdir -p ${system}/${variation}
    ln -sf \$(readlink -f ${top_file}) ${system}/${variation}/${system}_system_${variation}.${params.top_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep_${rep}.${params.traj_format}

    python ${scriptDir}/run_rms_analysis.py \\
        --systems '["${system}"]' \\
        --variations '{"${system}": ["${variation}"]}' \\
        --num-replicates 1 \\
        --analysis RMSF \\
        --target-selection '${params.target_selection}' \\
        --ref-selection '${params.ref_selection}' \\
        --start-frame ${params.start_frame} \\
        --traj-format ${params.traj_format} \\
        --top-format ${params.top_format} \\
        ${chain_arg} ${wrap_arg}

    for f in rmsf_plot_${system}_${variation}_rep1*.pkl; do
        newname=\$(echo "\$f" | sed "s/rep1/rep${rep}/")
        [ "\$f" != "\$newname" ] && mv "\$f" "\$newname" || true
    done
    """
}

process RUN_2D_RMSD {
    tag "${system}_${variation}_rep${rep}"
    label 'process_high'
    publishDir "${params.outdir}/2d_rmsd", mode: 'copy'

    input:
    tuple val(system), val(variation), val(rep), path(top_file), path(traj_file)

    output:
    path("2d_rmsd_plot_${system}_${variation}_rep${rep}.pkl"), emit: pickle

    script:
    def wrap_arg = params.wrap_selection == 'none' ? "--wrap-selection none" : (params.wrap_selection != 'auto' ? "--wrap-selection '${params.wrap_selection}'" : "")
    """
    mkdir -p ${system}/${variation}
    ln -sf \$(readlink -f ${top_file}) ${system}/${variation}/${system}_system_${variation}.${params.top_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep_${rep}.${params.traj_format}

    python ${scriptDir}/run_rms_analysis.py \\
        --systems '["${system}"]' \\
        --variations '{"${system}": ["${variation}"]}' \\
        --num-replicates 1 \\
        --analysis 2D-RMSD \\
        --target-selection '${params.target_selection}' \\
        --ref-selection '${params.ref_selection}' \\
        --start-frame ${params.start_frame} \\
        --traj-format ${params.traj_format} \\
        --top-format ${params.top_format} \\
        ${wrap_arg}

    mv 2d_rmsd_plot_${system}_${variation}_rep1.pkl 2d_rmsd_plot_${system}_${variation}_rep${rep}.pkl || true
    """
}

process RUN_ROG {
    tag "${system}_${variation}_rep${rep}"
    label 'process_medium'
    publishDir "${params.outdir}/rog", mode: 'copy'

    input:
    tuple val(system), val(variation), val(rep), path(top_file), path(traj_file)

    output:
    path("rog_plot_${system}_${variation}_rep${rep}.pkl"), emit: pickle

    script:
    def wrap_arg = params.wrap_selection == 'none' ? "--wrap-selection none" : (params.wrap_selection != 'auto' ? "--wrap-selection '${params.wrap_selection}'" : "")
    """
    mkdir -p ${system}/${variation}
    ln -sf \$(readlink -f ${top_file}) ${system}/${variation}/${system}_system_${variation}.${params.top_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep_${rep}.${params.traj_format}

    python ${scriptDir}/run_rog_analysis.py \\
        --systems '["${system}"]' \\
        --variations '{"${system}": ["${variation}"]}' \\
        --num-replicates 1 \\
        --start-frame ${params.start_frame} \\
        --traj-format ${params.traj_format} \\
        --top-format ${params.top_format} \\
        --selection '${params.rog_selection}' \\
        --time-unit ${params.time_unit} \\
        ${wrap_arg}

    mv rog_plot_${system}_${variation}_rep1.pkl rog_plot_${system}_${variation}_rep${rep}.pkl || true
    """
}

process RUN_HBONDS {
    tag "${system}_${variation}_rep${rep}"
    label 'process_high'
    cpus { params.parallel_backend != 'serial' && params.n_workers ? params.n_workers : (params.parallel_backend != 'serial' ? 4 : 1) }
    publishDir "${params.outdir}/hbonds", mode: 'copy'

    input:
    tuple val(system), val(variation), val(rep), path(top_file), path(traj_file)

    output:
    path("hbonds_plot_${system}_${variation}_rep${rep}.pkl"), emit: pickle

    script:
    def sel_arg = ""
    if (params.between_pairs) {
        sel_arg = "--between-pairs '${params.between_pairs}'"
    } else if (params.acceptors_sel && params.hydrogens_sel) {
        sel_arg = "--atom-selections '${params.acceptors_sel}' '${params.hydrogens_sel}'"
    } else {
        error "H-bonds analysis requires either --acceptors_sel + --hydrogens_sel or --between_pairs"
    }
    def update_arg = params.update_selections ? "" : "--no-update-selections"
    def wrap_arg = params.wrap_selection == 'none' ? "--wrap-selection none" : (params.wrap_selection != 'auto' ? "--wrap-selection '${params.wrap_selection}'" : "")
    def parallel_arg = params.parallel_backend != 'serial' ? "--parallel-backend ${params.parallel_backend}" : ""
    def nworkers_arg = params.n_workers ? "--n-workers ${params.n_workers}" : ""
    """
    mkdir -p ${system}/${variation}
    ln -sf \$(readlink -f ${top_file}) ${system}/${variation}/${system}_system_${variation}.${params.top_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep_${rep}.${params.traj_format}

    python ${scriptDir}/run_hbonds_analysis.py \\
        --systems '["${system}"]' \\
        --variations '{"${system}": ["${variation}"]}' \\
        --num-replicates 1 \\
        --start-frame ${params.start_frame} \\
        --traj-format ${params.traj_format} \\
        --top-format ${params.top_format} \\
        --d-a-cutoff ${params.d_a_cutoff} \\
        --d-h-a-angle-cutoff ${params.d_h_a_angle_cutoff} \\
        ${update_arg} \\
        ${sel_arg} ${wrap_arg} ${parallel_arg} ${nworkers_arg}

    mv hbonds_plot_${system}_${variation}_rep1.pkl hbonds_plot_${system}_${variation}_rep${rep}.pkl || true
    """
}

// ─── Plotting Processes ──────────────────────────────────────────────────────

process PLOT_RMSD {
    tag "${pickle.baseName}"
    label 'process_low'
    publishDir "${params.outdir}/plots/rmsd", mode: 'copy'

    input:
    path(pickle)

    output:
    path("*.png"), emit: plot

    script:
    """
    python ${scriptDir}/plotting/plot_rmsd.py \\
        --pickle-file ${pickle} \\
        --output-dir . \\
        --dpi ${params.dpi}
    """
}

process PLOT_RMSF {
    tag "${pickle.baseName}"
    label 'process_low'
    publishDir "${params.outdir}/plots/rmsf", mode: 'copy'

    input:
    path(pickle)

    output:
    path("*.png"), emit: plot

    script:
    """
    python ${scriptDir}/plotting/plot_rmsf.py \\
        --pickle-file ${pickle} \\
        --output-dir . \\
        --dpi ${params.dpi}
    """
}

process PLOT_2D_RMSD {
    tag "${pickle.baseName}"
    label 'process_low'
    publishDir "${params.outdir}/plots/2d_rmsd", mode: 'copy'

    input:
    path(pickle)

    output:
    path("*.png"), emit: plot

    script:
    """
    python ${scriptDir}/plotting/plot_2d_rmsd.py \\
        --pickle-file ${pickle} \\
        --output-dir . \\
        --dpi ${params.dpi}
    """
}

process PLOT_ROG {
    tag "${pickle.baseName}"
    label 'process_low'
    publishDir "${params.outdir}/plots/rog", mode: 'copy'

    input:
    path(pickle)

    output:
    path("*.png"), emit: plot

    script:
    """
    python ${scriptDir}/plotting/plot_rog.py \\
        --pickle-file ${pickle} \\
        --output-dir . \\
        --dpi ${params.dpi}
    """
}

process PLOT_HBONDS {
    tag "${pickle.baseName}"
    label 'process_low'
    publishDir "${params.outdir}/plots/hbonds", mode: 'copy'

    input:
    path(pickle)

    output:
    path("*.png"), emit: plot

    script:
    """
    python ${scriptDir}/plotting/plot_hbonds.py \\
        --pickle-file ${pickle} \\
        --output-dir . \\
        --dpi ${params.dpi} \\
        --time-unit ${params.time_unit} \\
        --top-n ${params.hbonds_top_n}
    """
}

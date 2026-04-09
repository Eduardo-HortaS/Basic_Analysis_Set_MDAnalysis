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
params.outdir           = 'results_nf'

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
params.strict_wrapping  = false

// Plotting
params.dpi              = 400
params.hbonds_top_n     = 20
params.rmsd_show_kde    = true

// Comparison plot groups (RMSD/RMSF)
params.plot_groups      = null  // JSON: {"group_name": [["SysA", "wild"], ...]}
params.replicate_mode   = 'separate'  // 'separate' or 'average'

// Parallelization — MDAnalysis 2.8+ native parallel backends
// 'serial' (default), 'multiprocessing', or 'dask'
params.parallel_backend = 'serial'
params.n_workers        = null    // null = auto-detect

// ─── Helper: Resolve script paths ────────────────────────────────────────────
def scriptDir = "${projectDir}"
def analysisDirAbs = new File(params.input_dir.toString()).absolutePath
def outdirAbs = resolve_outdir_abs(params.outdir, analysisDirAbs)

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

    def _plot_groups_map = parse_plot_groups(params.plot_groups)
    def _has_plot_groups = !_plot_groups_map.isEmpty()

    // Wrap each trajectory once and fan out wrapped inputs to all analyses.
    RUN_WRAP(analysis_inputs)
    def wrapped_inputs = RUN_WRAP.out.wrapped

    // Run analyses conditionally
    if (params.run_rmsd) {
        RUN_RMSD(wrapped_inputs)

        if (params.plot_rmsd) {
            PLOT_RMSD(RUN_RMSD.out.pickle)
            if (_has_plot_groups) {
                PLOT_RMSD_GROUPS(PLOT_RMSD.out.plot.collect())
            }
        }
    }

    if (params.run_rmsf) {
        RUN_RMSF(wrapped_inputs)

        if (params.plot_rmsf) {
            PLOT_RMSF(RUN_RMSF.out.pickle.flatten())
            if (_has_plot_groups) {
                PLOT_RMSF_GROUPS(PLOT_RMSF.out.plot.collect())
            }
        }
    }

    if (params.run_2d_rmsd) {
        RUN_2D_RMSD(wrapped_inputs)
        if (params.plot_2d_rmsd) {
            PLOT_2D_RMSD(RUN_2D_RMSD.out.pickle)
        }
    }

    if (params.run_rog) {
        RUN_ROG(wrapped_inputs)
        if (params.plot_rog) {
            PLOT_ROG(RUN_ROG.out.pickle)
        }
    }

    if (params.run_hbonds) {
        RUN_HBONDS(wrapped_inputs)
        if (params.plot_hbonds) {
            PLOT_HBONDS(RUN_HBONDS.out.pickle)
        }
    }
}

// ─── Helper function ─────────────────────────────────────────────────────────

import groovy.json.JsonOutput
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
                def traj = resolve_traj_file(system, variation, rep)
                def top  = file("${params.input_dir}/${system}/${variation}/${system}_system_${variation}.${params.top_format}")
                combos.add(tuple(system, variation, rep, top, traj))
            }
        }
    }
    return combos
}

def resolve_traj_file(system, variation, rep) {
    def withUnderscore = file("${params.input_dir}/${system}/${variation}/${system}_production_${variation}_rep_${rep}.${params.traj_format}")
    if (withUnderscore.exists()) {
        return withUnderscore
    }

    def noUnderscore = file("${params.input_dir}/${system}/${variation}/${system}_production_${variation}_rep${rep}.${params.traj_format}")
    if (noUnderscore.exists()) {
        return noUnderscore
    }

    // Keep legacy path as default so downstream errors show the historical pattern.
    return withUnderscore
}

def resolve_outdir_abs(rawOutdir, analysisDirAbs) {
    def configured = rawOutdir?.toString()?.trim()

    // Keep legacy defaults mapped to a deterministic analysis-local folder.
    if (!configured || configured in ['results', './results', 'results_nf', './results_nf']) {
        return new File(analysisDirAbs, 'results_nf').absolutePath
    }

    def configuredFile = new File(configured)
    if (configuredFile.isAbsolute()) {
        return configuredFile.absolutePath
    }

    // Resolve relative outdir from the analysis input folder, not launchDir.
    return new File(analysisDirAbs, configured).absolutePath
}

def parse_plot_groups(rawValue) {
    if (!rawValue) {
        return [:]
    }

    if (rawValue instanceof Map) {
        return rawValue as Map
    }

    def jsonSlurper = new JsonSlurper()
    def value = rawValue.toString()
    def maybeFile = new File(value)

    def parsed
    if (maybeFile.exists()) {
        parsed = jsonSlurper.parse(maybeFile)
    } else {
        parsed = jsonSlurper.parseText(value)
    }

    if (!(parsed instanceof Map)) {
        error "ERROR: plot_groups must be a JSON object mapping group names to member lists"
    }
    return parsed as Map
}

// ─── Processes ───────────────────────────────────────────────────────────────

process RUN_WRAP {
    tag "${system}_${variation}_rep${rep}"
    label 'process_medium'
    cache 'lenient'
    publishDir "${outdirAbs}/wrapped", mode: 'copy'

    input:
    tuple val(system), val(variation), val(rep), path(top_file), path(traj_file)

    output:
    tuple val(system), val(variation), val(rep), path(top_file), path("wrapped_${system}_production_${variation}_rep${rep}.${params.traj_format}"), emit: wrapped

    script:
    def wrappedName = "wrapped_${system}_production_${variation}_rep${rep}.${params.traj_format}"
    """
    persisted_dir="${outdirAbs}/wrapped"
    persisted_file="\$persisted_dir/${wrappedName}"

    if [ -s "\$persisted_file" ]; then
        cp -f "\$persisted_file" "${wrappedName}"
    else
        python - "${top_file}" "${traj_file}" "${wrappedName}" "${params.wrap_selection}" "${params.strict_wrapping}" <<'PY'
import shutil
import sys

import MDAnalysis as mda

from utils import build_complex_selection, transform_trajectory

top_file, traj_file, wrapped_out, wrap_selection, strict_raw = sys.argv[1:6]
strict_wrapping = strict_raw.lower() in {'true', '1', 'yes'}

if wrap_selection.lower() == 'none':
    shutil.copy2(traj_file, wrapped_out)
    sys.exit(0)

u = mda.Universe(top_file, traj_file)
complex_ag, ligand_ag, rest_ag = build_complex_selection(u, wrap_selection=wrap_selection)

if complex_ag is None:
    shutil.copy2(traj_file, wrapped_out)
else:
    transform_trajectory(
        u,
        complex_ag,
        rest_ag,
        ligand_selection=ligand_ag,
        strict_wrapping=strict_wrapping,
    )
    with mda.Writer(wrapped_out, n_atoms=u.atoms.n_atoms) as writer:
        for _ts in u.trajectory:
            writer.write(u.atoms)
PY

        cp -f "${wrappedName}" "\$persisted_file" 2>/dev/null || true
    fi
    """
}

process RUN_RMSD {
    tag "${system}_${variation}_rep${rep}"
    label 'process_medium'
    cpus { params.parallel_backend != 'serial' && params.n_workers ? params.n_workers : (params.parallel_backend != 'serial' ? 4 : 1) }
    publishDir "${outdirAbs}/rmsd", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(system), val(variation), val(rep), path(top_file), path(traj_file)

    output:
    path("rmsd_plot_${system}_${variation}_rep${rep}*.pkl"), emit: pickle

    script:
    def time_arg = params.time_interval_between_frames ? "--time-interval-between-frames ${params.time_interval_between_frames} --time-unit ${params.time_unit}" : ""
    def group_selections_list = []
    if (params.group_selections) {
        if (params.group_selections instanceof Collection) {
            group_selections_list = params.group_selections as List
        } else {
            def raw = params.group_selections.toString().trim()
            if (raw.startsWith('[') && raw.endsWith(']')) {
                try {
                    group_selections_list = new JsonSlurper().parseText(raw) as List
                } catch (Exception ignored) {
                    // Support legacy non-JSON list syntax like:
                    // [chainid A, chainid B, chainid C]
                    def inner = raw.substring(1, raw.length() - 1).trim()
                    if (inner) {
                        group_selections_list = inner
                            .split(/\s*,\s*/)
                            .collect { it.toString().replaceAll(/^['\"]|['\"]$/, '') }
                            .findAll { it }
                    }
                }
            } else {
                group_selections_list = [raw]
            }
        }
    }
    def group_arg = group_selections_list ? "--group-selections ${group_selections_list.collect { "'${it.toString().replace("'", "'\\''")}'" }.join(' ')}" : ""
    def wrap_arg = "--wrap-selection none"
    def strict_arg = ""
    def parallel_arg = params.parallel_backend != 'serial' ? "--parallel-backend ${params.parallel_backend}" : ""
    def nworkers_arg = params.n_workers ? "--n-workers ${params.n_workers}" : ""
    """
    mkdir -p ${system}/${variation}
    ln -sf \$(readlink -f ${top_file}) ${system}/${variation}/${system}_system_${variation}.${params.top_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep_1.${params.traj_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep1.${params.traj_format}

    cp -f ${outdirAbs}/rmsd/rmsd_plot_${system}_${variation}_rep${rep}*.pkl . 2>/dev/null || true

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
        ${time_arg} ${group_arg} ${wrap_arg} ${strict_arg} ${parallel_arg} ${nworkers_arg}

    for f in rmsd_plot_${system}_${variation}_rep1*.pkl; do
        if [ -e "\$f" ]; then
            newname=\$(echo "\$f" | sed "s/rep1/rep${rep}/")
            [ "\$f" != "\$newname" ] && mv "\$f" "\$newname" || true
        fi
    done
    """
}

process RUN_RMSF {
    tag "${system}_${variation}_rep${rep}"
    label 'process_medium'
    publishDir "${outdirAbs}/rmsf", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(system), val(variation), val(rep), path(top_file), path(traj_file)

    output:
    path("rmsf_plot_${system}_${variation}_rep${rep}*.pkl"), emit: pickle

    script:
    def chain_arg = params.chain_intervals ? "--chain-intervals '${params.chain_intervals}'" : ""
    def group_selections_list = []
    if (params.group_selections) {
        if (params.group_selections instanceof Collection) {
            group_selections_list = params.group_selections as List
        } else {
            def raw = params.group_selections.toString().trim()
            if (raw.startsWith('[') && raw.endsWith(']')) {
                try {
                    group_selections_list = new JsonSlurper().parseText(raw) as List
                } catch (Exception ignored) {
                    // Support legacy non-JSON list syntax like:
                    // [chainid A, chainid B, chainid C]
                    def inner = raw.substring(1, raw.length() - 1).trim()
                    if (inner) {
                        group_selections_list = inner
                            .split(/\s*,\s*/)
                            .collect { it.toString().replaceAll(/^['\"]|['\"]$/, '') }
                            .findAll { it }
                    }
                }
            } else {
                group_selections_list = [raw]
            }
        }
    }
    def group_arg = group_selections_list ? "--group-selections ${group_selections_list.collect { "'${it.toString().replace("'", "'\\''")}'" }.join(' ')}" : ""
    def wrap_arg = "--wrap-selection none"
    def strict_arg = ""
    """
    mkdir -p ${system}/${variation}
    ln -sf \$(readlink -f ${top_file}) ${system}/${variation}/${system}_system_${variation}.${params.top_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep_1.${params.traj_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep1.${params.traj_format}

    cp -f ${outdirAbs}/rmsf/rmsf_plot_${system}_${variation}_rep${rep}*.pkl . 2>/dev/null || true

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
        ${chain_arg} ${group_arg} ${wrap_arg} ${strict_arg}

    for f in rmsf_plot_${system}_${variation}_rep1*.pkl; do
        newname=\$(echo "\$f" | sed "s/rep1/rep${rep}/")
        [ "\$f" != "\$newname" ] && mv "\$f" "\$newname" || true
    done
    """
}

process RUN_2D_RMSD {
    tag "${system}_${variation}_rep${rep}"
    label 'process_high'
    publishDir "${outdirAbs}/2d_rmsd", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(system), val(variation), val(rep), path(top_file), path(traj_file)

    output:
    path("2d_rmsd_plot_${system}_${variation}_rep${rep}.pkl"), emit: pickle

    script:
    def wrap_arg = "--wrap-selection none"
    def strict_arg = ""
    """
    mkdir -p ${system}/${variation}
    ln -sf \$(readlink -f ${top_file}) ${system}/${variation}/${system}_system_${variation}.${params.top_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep_1.${params.traj_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep1.${params.traj_format}

    cp -f ${outdirAbs}/2d_rmsd/2d_rmsd_plot_${system}_${variation}_rep${rep}.pkl . 2>/dev/null || true

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
        ${wrap_arg} ${strict_arg}

    mv 2d_rmsd_plot_${system}_${variation}_rep1.pkl 2d_rmsd_plot_${system}_${variation}_rep${rep}.pkl || true
    """
}

process RUN_ROG {
    tag "${system}_${variation}_rep${rep}"
    label 'process_medium'
    publishDir "${outdirAbs}/rog", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(system), val(variation), val(rep), path(top_file), path(traj_file)

    output:
    path("rog_plot_${system}_${variation}_rep${rep}.pkl"), emit: pickle

    script:
    def wrap_arg = "--wrap-selection none"
    def strict_arg = ""
    """
    mkdir -p ${system}/${variation}
    ln -sf \$(readlink -f ${top_file}) ${system}/${variation}/${system}_system_${variation}.${params.top_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep_1.${params.traj_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep1.${params.traj_format}

    cp -f ${outdirAbs}/rog/rog_plot_${system}_${variation}_rep${rep}.pkl . 2>/dev/null || true

    python ${scriptDir}/run_rog_analysis.py \\
        --systems '["${system}"]' \\
        --variations '{"${system}": ["${variation}"]}' \\
        --num-replicates 1 \\
        --start-frame ${params.start_frame} \\
        --traj-format ${params.traj_format} \\
        --top-format ${params.top_format} \\
        --selection '${params.rog_selection}' \\
        --time-unit ${params.time_unit} \\
        ${wrap_arg} ${strict_arg}

    mv rog_plot_${system}_${variation}_rep1.pkl rog_plot_${system}_${variation}_rep${rep}.pkl || true
    """
}

process RUN_HBONDS {
    tag "${system}_${variation}_rep${rep}"
    label 'process_high'
    cpus { params.parallel_backend != 'serial' && params.n_workers ? params.n_workers : (params.parallel_backend != 'serial' ? 4 : 1) }
    publishDir "${outdirAbs}/hbonds", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(system), val(variation), val(rep), path(top_file), path(traj_file)

    output:
    path("hbonds_plot_${system}_${variation}_rep${rep}.pkl"), emit: pickle

    script:
    def sel_arg = ""
    def atom_sel_arg = ""
    if (params.acceptors_sel && params.hydrogens_sel) {
        atom_sel_arg = "--atom-selections '${params.acceptors_sel}' '${params.hydrogens_sel}'"
    }

    if (params.between_pairs) {
        def betweenPairsJson = (params.between_pairs instanceof Collection || params.between_pairs instanceof Map)
            ? JsonOutput.toJson(params.between_pairs)
            : params.between_pairs.toString()
        sel_arg = "--between-pairs '${betweenPairsJson}' ${atom_sel_arg}"
    } else if (atom_sel_arg) {
        sel_arg = atom_sel_arg
    } else {
        error "H-bonds analysis requires either --acceptors_sel + --hydrogens_sel or --between_pairs"
    }
    def update_arg = params.update_selections ? "" : "--no-update-selections"
    def wrap_arg = "--wrap-selection none"
    def strict_arg = ""
    def parallel_arg = params.parallel_backend != 'serial' ? "--parallel-backend ${params.parallel_backend}" : ""
    def nworkers_arg = params.n_workers ? "--n-workers ${params.n_workers}" : ""
    """
    mkdir -p ${system}/${variation}
    ln -sf \$(readlink -f ${top_file}) ${system}/${variation}/${system}_system_${variation}.${params.top_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep_1.${params.traj_format}
    ln -sf \$(readlink -f ${traj_file}) ${system}/${variation}/${system}_production_${variation}_rep1.${params.traj_format}

    cp -f ${outdirAbs}/hbonds/hbonds_plot_${system}_${variation}_rep${rep}.pkl . 2>/dev/null || true

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
        ${sel_arg} ${wrap_arg} ${strict_arg} ${parallel_arg} ${nworkers_arg}

    mv hbonds_plot_${system}_${variation}_rep1.pkl hbonds_plot_${system}_${variation}_rep${rep}.pkl || true
    """
}

// ─── Plotting Processes ──────────────────────────────────────────────────────

process PLOT_RMSD {
    tag "${pickle.baseName}"
    label 'process_low'
    publishDir "${outdirAbs}/plots/rmsd", mode: 'copy'
    cache 'lenient'

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

process PLOT_RMSD_GROUPS {
    label 'process_low'
    publishDir "${outdirAbs}/plots/rmsd", mode: 'copy'
    cache 'lenient'

    input:
    path(plot_pngs)

    output:
    path("rmsd_comparison_*.png"), emit: plot, optional: true

    script:
    def plotGroupsArg = (params.plot_groups instanceof Map || params.plot_groups instanceof Collection)
        ? JsonOutput.toJson(params.plot_groups)
        : params.plot_groups
    """
    python ${scriptDir}/plot_group_comparisons.py \\
        --analysis rmsd \\
        --plot-groups '${plotGroupsArg}' \\
        --work-dir '${outdirAbs}/rmsd' \\
        --output-dir . \\
        --num-replicates ${params.num_replicates} \\
        --replicate-mode '${params.replicate_mode}' \\
        --dpi ${params.dpi} \\
        --rmsd-show-kde ${params.rmsd_show_kde}
    """
}

process PLOT_RMSF {
    tag "${pickle.baseName}"
    label 'process_low'
    publishDir "${outdirAbs}/plots/rmsf", mode: 'copy'
    cache 'lenient'

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

process PLOT_RMSF_GROUPS {
    label 'process_low'
    publishDir "${outdirAbs}/plots/rmsf", mode: 'copy'
    cache 'lenient'

    input:
    path(plot_pngs)

    output:
    path("rmsf_comparison_*.png"), emit: plot, optional: true

    script:
    def plotGroupsArg = (params.plot_groups instanceof Map || params.plot_groups instanceof Collection)
        ? JsonOutput.toJson(params.plot_groups)
        : params.plot_groups
    """
    python ${scriptDir}/plot_group_comparisons.py \\
        --analysis rmsf \\
        --plot-groups '${plotGroupsArg}' \\
        --work-dir '${outdirAbs}/rmsf' \\
        --output-dir . \\
        --num-replicates ${params.num_replicates} \\
        --replicate-mode '${params.replicate_mode}' \\
        --target-selection '${params.target_selection}' \\
        --dpi ${params.dpi}
    """
}

process PLOT_2D_RMSD {
    tag "${pickle.baseName}"
    label 'process_low'
    publishDir "${outdirAbs}/plots/2d_rmsd", mode: 'copy'
    cache 'lenient'

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
    publishDir "${outdirAbs}/plots/rog", mode: 'copy'
    cache 'lenient'

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
    publishDir "${outdirAbs}/plots/hbonds", mode: 'copy'
    cache 'lenient'

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

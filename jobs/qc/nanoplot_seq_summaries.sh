#!/bin/bash -l

#============================================================================================
# Title: NanoPlot and NanoComp Analysis
# Description: Performs QC analysis on ONT sequencing data using NanoPlot
#              and NanoComp tools
# Author: Renato Santos
# Usage: sbatch nanoplot_seq_summaries.sh
#============================================================================================

#SBATCH --job-name=nanoplot_seq_summaries
#SBATCH --partition=nd_bioinformatics_cpu,cpu,interruptible_cpu
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-01:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/ont-benchmark/jobs/qc

# Constants
readonly BASE_DIR="/scratch/prj/ppn_als_longread/ont-benchmark"
readonly SUMMARY_DIR="${BASE_DIR}/seq_summaries"
readonly NANOPLOT_OUTPUT_DIR="${BASE_DIR}/qc/nanoplot/seq_summaries"
readonly NANOCOMP_OUTPUT_DIR="${BASE_DIR}/qc/nanocomp/seq_summaries"
readonly NANOPLOT_IMAGE="docker://quay.io/biocontainers/nanoplot:1.42.0--pyhdfd78af_0"
readonly NANOCOMP_IMAGE="docker://quay.io/biocontainers/nanocomp:1.23.1--pyhdfd78af_0"

# Environment setup
export SINGULARITY_CACHEDIR="/scratch/users/${USER}/singularity/"


# Function to run NanoPlot analysis
run_nanoplot() {
    local summary_file="$1"
    local output_dir="$2"
    local is_barcoded="$3"

    local base_cmd=(
        srun -n1 -N1 singularity exec --contain --bind "${BASE_DIR}"
        "${NANOPLOT_IMAGE}"
        NanoPlot
        --summary "${summary_file}"
        --outdir "${output_dir}"
        --N50
        --loglength
        --color royalblue
        --colormap Viridis
        --plots kde
        --dpi 300
        --tsv_stats
        --info_in_report
        --threads 2
    )

    if [[ "${is_barcoded}" == true ]]; then
        base_cmd+=("--barcoded")
    fi

    "${base_cmd[@]}" &
}


# Function to run NanoComp analysis
run_nanocomp() {
    local summary_files=("$1")
    local summary_names=("$2")

    srun -n1 -N1 singularity exec --contain --bind "${BASE_DIR}" \
        "${NANOCOMP_IMAGE}" \
        NanoComp \
        --summary ${summary_files} \
        --outdir "${NANOCOMP_OUTPUT_DIR}" \
        --plot violin \
        --names ${summary_names} \
        --threads 2 &
}


main() {
    local summary_files=()
    local summary_names=()

    # Check if summary directory exists
    if [[ ! -d "${SUMMARY_DIR}" ]]; then
        echo "Error: Summary directory ${SUMMARY_DIR} does not exist"
        exit 1
    fi

    # Process each summary file
    for summary_file in "${SUMMARY_DIR}"/*/*.txt; do
        if [[ -f "${summary_file}" ]]; then
            local subdir_name
            subdir_name=$(basename "$(dirname "${summary_file}")")
            local nanoplot_outdir="${NANOPLOT_OUTPUT_DIR}/${subdir_name}"

            summary_files+=("${summary_file}")
            summary_names+=("${subdir_name}")

            mkdir -p "${nanoplot_outdir}"

            if [[ "${subdir_name}" == *__* ]]; then
                run_nanoplot "${summary_file}" "${nanoplot_outdir}" true
            else
                run_nanoplot "${summary_file}" "${nanoplot_outdir}" false
            fi
        fi
    done

    # Run NanoComp if we have files to process
    if [[ ${#summary_files[@]} -gt 0 ]]; then
        mkdir -p "${NANOCOMP_OUTPUT_DIR}"
        run_nanocomp "${summary_files[*]}" "${summary_names[*]}"
    else
        echo "Warning: No summary files found to process"
    fi

    wait
}


main

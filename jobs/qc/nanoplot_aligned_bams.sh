#!/bin/bash -l

#============================================================================================
# Title: NanoPlot and NanoComp Analysis for Aligned BAMs
# Description: Performs QC analysis on ONT aligned BAM/CRAM files using NanoPlot
#              and NanoComp tools
# Author: Renato Santos
# Usage: sbatch nanoplot_aligned_bams.sh
#============================================================================================

#SBATCH --job-name=nanoplot_aligned_bams
#SBATCH --partition=nd_bioinformatics_cpu,cpu,interruptible_cpu
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=512M
#SBATCH --time=0-06:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/ont-benchmark/jobs/qc

# Constants
readonly BASE_DIR="/scratch/prj/ppn_als_longread/ont-benchmark"
readonly CRAM_DIR="${BASE_DIR}/vcf"
readonly NANOPLOT_OUTPUT_DIR="${BASE_DIR}/qc/nanoplot/aligned_bams"
readonly NANOCOMP_OUTPUT_DIR="${BASE_DIR}/qc/nanocomp/aligned_bams"
readonly NANOPLOT_IMAGE="docker://quay.io/biocontainers/nanoplot:1.42.0--pyhdfd78af_0"
readonly NANOCOMP_IMAGE="docker://quay.io/biocontainers/nanocomp:1.23.1--pyhdfd78af_0"

# Environment setup
export SINGULARITY_CACHEDIR="/scratch/users/${USER}/singularity/"


# Function to join array elements with a space
join_by_space() {
    local IFS=' '
    echo "$*"
}


# Function to run NanoPlot analysis
run_nanoplot() {
    local cram_file="$1"
    local output_dir="$2"
    local pickle_file="${output_dir}/NanoPlot-data.pickle"

    # Determine if we should use pickle or cram
    local input_cmd
    if [[ -f "${pickle_file}" ]]; then
        echo "Using existing pickle file for ${output_dir}"
        input_cmd="--pickle ${pickle_file}"
    else
        echo "No pickle file found, using CRAM file for ${output_dir}"
        input_cmd="--cram ${cram_file}"
    fi

    srun -n1 -N1 singularity exec --contain --bind "${BASE_DIR}" \
        "${NANOPLOT_IMAGE}" \
        NanoPlot \
            ${input_cmd} \
            --outdir "${output_dir}" \
            --N50 \
            --loglength \
            --color royalblue \
            --colormap Viridis \
            --plots kde \
            --dpi 300 \
            --tsv_stats \
            --info_in_report \
            --store \
            --threads 8 &
}


# Function to run NanoComp analysis
run_nanocomp() {
    local -a pickle_files=("${@:1:$((($#)/2))}")
    local -a pickle_names=("${@:$((($#)/2+1))}")

    local pickle_files_str
    local pickle_names_str
    pickle_files_str=$(join_by_space "${pickle_files[@]}")
    pickle_names_str=$(join_by_space "${pickle_names[@]}")

    if [[ ${#pickle_files[@]} -gt 0 ]]; then
        mkdir -p "${NANOCOMP_OUTPUT_DIR}"
        echo "Running NanoComp with the following pickle files:"
        echo "Files: ${pickle_files_str}"
        echo "Names: ${pickle_names_str}"

        srun -n1 -N1 singularity exec --contain --bind "${BASE_DIR}" \
            "${NANOCOMP_IMAGE}" \
            NanoComp \
                --pickle ${pickle_files_str} \
                --outdir "${NANOCOMP_OUTPUT_DIR}" \
                --plot violin \
                --names ${pickle_names_str} \
                --threads 8 &
    fi
}


main() {
    local -a pickle_files=()
    local -a pickle_names=()

    # Check if CRAM directory exists
    if [[ ! -d "${CRAM_DIR}" ]]; then
        echo "Error: CRAM directory ${CRAM_DIR} does not exist"
        exit 1
    fi

    # Process each SUP CRAM file
    for cram_file in "${CRAM_DIR}"/*_sup/*.haplotagged.cram; do
        if [[ -f "${cram_file}" ]]; then
            local subdir_name
            subdir_name=$(basename "$(dirname "${cram_file}")")
            subdir_name=${subdir_name%.haplotagged}
            local nanoplot_outdir="${NANOPLOT_OUTPUT_DIR}/${subdir_name}"
            local pickle_file="${nanoplot_outdir}/NanoPlot-data.pickle"

            mkdir -p "${nanoplot_outdir}"

            # Add to arrays for NanoComp
            pickle_files+=("${pickle_file}")
            pickle_names+=("${subdir_name}")

            run_nanoplot "${cram_file}" "${nanoplot_outdir}"
        fi
    done

    wait

    if [[ ${#pickle_files[@]} -gt 0 ]]; then
        run_nanocomp "${pickle_files[@]}" "${pickle_names[@]}"
    else
        echo "No SUP samples found to process"
    fi

    wait
}


main

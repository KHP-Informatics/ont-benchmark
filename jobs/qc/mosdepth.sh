#!/bin/bash -l

#============================================================================================
# Title: Mosdepth Analysis for Aligned CRAMs
# Description: Performs coverage analysis on ONT aligned CRAM files using mosdepth
# Author: Renato Santos
# Usage: sbatch mosdepth_analysis.sh
#============================================================================================

#SBATCH --job-name=mosdepth
#SBATCH --partition=nd_bioinformatics_cpu,cpu,interruptible_cpu
#SBATCH --ntasks=14
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-01:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/ont-benchmark/jobs/qc

# Constants
readonly BASE_DIR="/scratch/prj/ppn_als_longread/ont-benchmark"
readonly INPUT_DIR="${BASE_DIR}/vcf"
readonly OUTPUT_BASE_DIR="${BASE_DIR}/qc/mosdepth"
readonly REFERENCE="${BASE_DIR}/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
readonly MOSDEPTH_IMAGE="docker://quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0"

# Environment setup
export SINGULARITY_CACHEDIR="/scratch/users/${USER}/singularity/"


# Function to run mosdepth analysis
run_mosdepth() {
    local input_cram="$1"
    local output_dir="$2"
    local prefix="$3"

    echo "Processing ${input_cram} with prefix ${prefix}"

    srun -n1 -N1 singularity exec --contain --bind "${BASE_DIR}" \
        "${MOSDEPTH_IMAGE}" \
        bash -c "cd ${output_dir} && \
        mosdepth \
            --threads 4 \
            --fasta ${REFERENCE} \
            --fast-mode \
            ${prefix} \
            ${input_cram}" & }


process_cram_files() {
    for input_cram in "${INPUT_DIR}"/*/*.haplotagged.cram; do
        if [[ -f "${input_cram}" ]]; then
            cram_subdir=$(basename "$(dirname "${input_cram}")")
            output_dir="${OUTPUT_BASE_DIR}/${cram_subdir}"
            prefix="${output_dir}/$(basename "${input_cram}" .haplotagged.cram)"

            # Ensure output directory exists
            mkdir -p "${output_dir}"

            run_mosdepth "${input_cram}" "${output_dir}" "${prefix}"
        fi
    done
}


main() {
    # Check if input directory exists
    if [[ ! -d "${INPUT_DIR}" ]]; then
        echo "Error: Input directory ${INPUT_DIR} does not exist"
        exit 1
    fi

    process_cram_files
    wait
}


main

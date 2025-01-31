#!/bin/bash -l

#============================================================================================
# Title: Picard LiftOver for VCF Files
# Description: Performs liftover of VCF files from GRCh37 to GRCh38 using Picard tools
# Author: Renato Santos
# Usage: sbatch picard_liftover.sh
#============================================================================================

#SBATCH --job-name=picard_liftover
#SBATCH --partition=nd_bioinformatics_cpu,cpu
#SBATCH --ntasks=42
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=0-01:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/ont-benchmark/jobs/illumina

# Constants
readonly BASE_DIR="/scratch/prj/ppn_als_longread"
readonly REFERENCE="${BASE_DIR}/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
readonly REFERENCE_DICT="${REFERENCE%.fna}.dict"
readonly INPUT_DIR="${BASE_DIR}/benchmark_data/illumina/grch37"
readonly OUTPUT_DIR="${BASE_DIR}/benchmark_data/illumina/grch38"
readonly REJECT_DIR="${BASE_DIR}/benchmark_data/illumina/picard_liftover/rejected_variants"
readonly LOG_DIR="${BASE_DIR}/benchmark_data/illumina/picard_liftover/logs"
readonly TMP_DIR="${BASE_DIR}/benchmark_data/illumina/picard_liftover/work"
readonly CHAIN="${BASE_DIR}/references/hg19ToHg38.over.chain.gz"
readonly PICARD_IMAGE="docker://broadinstitute/picard:3.1.1"

# Environment setup
export SINGULARITY_CACHEDIR="/scratch/users/${USER}/singularity/"


# Function to create sequence dictionary
create_sequence_dictionary() {
    if [[ ! -f "${REFERENCE_DICT}" ]]; then
        echo "Dictionary file not found. Creating dictionary file..."
        singularity exec --contain --bind /scratch "${PICARD_IMAGE}" \
            java -jar /usr/picard/picard.jar CreateSequenceDictionary \
            -REFERENCE "${REFERENCE}" \
            -OUTPUT "${REFERENCE_DICT}" \
            -TMP_DIR "${TMP_DIR}"
    fi
}


# Function to process a VCF file
liftover_vcf() {
    local input_vcf="$1"
    local output_vcf="${OUTPUT_DIR}/$(basename "${input_vcf}")"
    local reject_vcf="${REJECT_DIR}/$(basename "${input_vcf}")"
    local log_file="${LOG_DIR}/$(basename "${input_vcf}").log"

    echo "Processing ${input_vcf} and logging to ${log_file}"

    srun -n1 -N1 singularity exec --contain --bind /scratch "${PICARD_IMAGE}" \
        java -jar /usr/picard/picard.jar LiftoverVcf \
        -INPUT "${input_vcf}" \
        -OUTPUT "${output_vcf}" \
        -REJECT "${reject_vcf}" \
        -REFERENCE_SEQUENCE "${REFERENCE}" \
        -CHAIN "${CHAIN}" \
        -TMP_DIR "${TMP_DIR}" \
        -WARN_ON_MISSING_CONTIG true \
        -RECOVER_SWAPPED_REF_ALT true &> "${log_file}"
}


# Function to create required directories
create_directories() {
    mkdir -p "${OUTPUT_DIR}"
    mkdir -p "${REJECT_DIR}"
    mkdir -p "${LOG_DIR}"
    mkdir -p "${TMP_DIR}"
}


main() {
    create_directories

    create_sequence_dictionary

    # Export variables needed by parallel processes
    export OUTPUT_DIR REJECT_DIR REFERENCE CHAIN TMP_DIR LOG_DIR
    export -f liftover_vcf

    # Process all VCF files in parallel
    find "${INPUT_DIR}" -name 'LP*-DNA_*.vcf.gz' | \
        xargs -I{} -P 0 srun -n1 -N1 bash -c "liftover_vcf {}"

    wait
}

main

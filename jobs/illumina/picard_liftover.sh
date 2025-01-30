#!/bin/bash -l

#SBATCH --job-name=picard_liftover
#SBATCH --partition=nd_bioinformatics_cpu,cpu
#SBATCH --ntasks=42
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=0-01:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/jobs/qc

# Set Singularity cache directory
export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/

export TMP_DIR=/scratch/prj/ppn_als_longread/benchmark_data/illumina/picard_liftover/work

# Check if dictionary file exists, if not, create it
REFERENCE=/scratch/prj/ppn_als_longread/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
REFERENCE_INDEX="${REFERENCE}.fai"
REFERENCE_DICT=${REFERENCE%.fna}.dict

if [ ! -f ${REFERENCE_DICT} ]; then
    echo "Dictionary file not found. Creating dictionary file..."

    # docker://broadinstitute/picard:3.1.1
    singularity exec --bind /scratch:/scratch /scratch/users/k2474617/containers/picard_3.1.1.sif \
        java -jar /usr/picard/picard.jar CreateSequenceDictionary \
        -REFERENCE ${REFERENCE} \
        -OUTPUT ${REFERENCE_DICT} \
        -TMP_DIR ${TMP_DIR}
fi

mkdir -p /scratch/prj/ppn_als_longread/benchmark_data/illumina/grch38/
mkdir -p /scratch/prj/ppn_als_longread/benchmark_data/illumina/picard_liftover/rejected_variants/
mkdir -p /scratch/prj/ppn_als_longread/benchmark_data/illumina/picard_liftover/logs/

INPUT_DIR="/scratch/prj/ppn_als_longread/benchmark_data/illumina/grch37"
OUTPUT_DIR="/scratch/prj/ppn_als_longread/benchmark_data/illumina/grch38"
REJECT_DIR="/scratch/prj/ppn_als_longread/benchmark_data/illumina/picard_liftover/rejected_variants"
LOG_DIR="/scratch/prj/ppn_als_longread/benchmark_data/illumina/picard_liftover/logs"
CHAIN="/scratch/prj/ppn_als_longread/benchmark_data/illumina/picard_liftover/chain/hg19ToHg38.over.chain"

process_vcf() {
    local input_vcf=$1
    local output_vcf=${OUTPUT_DIR}/$(basename ${input_vcf})
    local reject_vcf=${REJECT_DIR}/$(basename ${input_vcf})
    local log_file=${LOG_DIR}/$(basename ${input_vcf}).log

    echo "Processing ${input_vcf} and logging to ${log_file}"

    # docker://broadinstitute/picard:3.1.1
    srun -n1 -N1 singularity exec --bind /scratch:/scratch /scratch/users/k2474617/containers/picard_3.1.1.sif \
        java -jar /usr/picard/picard.jar LiftoverVcf \
        -INPUT ${input_vcf} \
        -OUTPUT ${output_vcf} \
        -REJECT ${reject_vcf} \
        -REFERENCE_SEQUENCE ${REFERENCE} \
        -CHAIN ${CHAIN} \
        -TMP_DIR ${TMP_DIR} \
        -WARN_ON_MISSING_CONTIG true \
        -RECOVER_SWAPPED_REF_ALT true &> ${log_file}
}

export -f process_vcf
export OUTPUT_DIR REJECT_DIR REFERENCE CHAIN TMP_DIR LOG_DIR

find ${INPUT_DIR} -name 'LP*-DNA_*.vcf.gz' | \
    xargs -I{} -P 0 srun -n1 -N1 bash -c "process_vcf {}"

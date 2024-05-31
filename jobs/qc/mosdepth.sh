#!/bin/bash -l

#SBATCH --job-name=mosdepth
#SBATCH --partition=nd_bioinformatics_cpu,cpu,interruptible_cpu
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-00:20:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/jobs/qc

# Set Singularity cache directory
export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/

# Define directories and reference
INPUT_DIR="/scratch/prj/ppn_als_longread/vcf/*/"
OUTPUT_BASE_DIR="/scratch/prj/ppn_als_longread/qc/mosdepth/"
REFERENCE="/scratch/prj/ppn_als_longread/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Create output directories and process files
for input_cram in ${INPUT_DIR}*.haplotagged.cram; do
    cram_subdir=$(basename $(dirname ${input_cram}))
    prefix=${OUTPUT_BASE_DIR}/${cram_subdir}/$(basename ${input_cram} .haplotagged.cram)
    output_dir=${OUTPUT_BASE_DIR}/${cram_subdir}

    mkdir -p ${output_dir}

    echo "Processing ${input_cram} with prefix ${prefix}"

    srun -n1 -N1 singularity exec --bind /scratch:/scratch docker://quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0 bash -c "\
        cd ${output_dir} && \
        mosdepth \
            --threads 4 \
            --no-per-base \
            --by 1000 \
            --fasta ${REFERENCE} \
            --fast-mode \
            $(basename ${prefix}) \
            ${input_cram}" &
done

wait

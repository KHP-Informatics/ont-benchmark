#!/bin/bash -l

#SBATCH --job-name=A048_09_sup_wf-human-variation
#SBATCH --partition=nd_bioinformatics_cpu,cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/jobs/A048_09_sup

module load nextflow/24.10.2-gcc-13.2.0

export NXF_HOME=/scratch/users/${USER}/nextflow/
export NXF_CACHE=/scratch/users/${USER}/nextflow/cache
export NXF_TEMP=/scratch/users/${USER}/nextflow/tmp
export NXF_SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/
export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

nextflow run epi2me-labs/wf-human-variation \
    -r v2.1.0 \
    -c /scratch/prj/ppn_als_longread/config/wf-human-variation.config \
    --bam /scratch/prj/ppn_als_longread/basecalled/A048_09_sup \
    --out_dir /scratch/prj/ppn_als_longread/vcf/A048_09_sup \
    --sample_name A048_09_sup \
    --sex female \
    --basecaller_cfg dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
    --cnv false
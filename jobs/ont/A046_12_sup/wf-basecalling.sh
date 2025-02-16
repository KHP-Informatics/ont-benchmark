#!/bin/bash -l

#SBATCH --job-name=A046_12_sup_wf-basecalling
#SBATCH --partition=nd_bioinformatics_cpu,cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch_tmp/prj/ppn_als_longread/wf-basecalling/jobs/ont-benchmark/A046_12_sup

module load nextflow/24.10.2-gcc-13.2.0

export NXF_HOME=/scratch/users/${USER}/nextflow/
export NXF_CACHE=/scratch/users/${USER}/nextflow/cache
export NXF_TEMP=/scratch/users/${USER}/nextflow/tmp
export NXF_SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/
export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

nextflow run epi2me-labs/wf-basecalling \
    -r v1.4.5 \
    -c /scratch_tmp/prj/ppn_als_longread/wf-basecalling/config/wf-basecalling.config \
    --input /scratch_tmp/prj/ppn_als_longread/wf-basecalling/pod5/ont-benchmark/A046_12 \
    --out_dir /scratch_tmp/prj/ppn_als_longread/wf-basecalling/basecalled/A046_12_sup \
    --sample_name A046_12_sup \
    --basecaller_cfg dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    --remora_cfg dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v3

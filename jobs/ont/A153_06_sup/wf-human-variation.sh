#!/bin/bash -l

#SBATCH --job-name=A153_06_sup_wf-human-variation
#SBATCH --partition=nd_bioinformatics_cpu,cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch_tmp/prj/ppn_als_longread/wf-human-variation/jobs/ont-benchmark/A153_06_sup

module load nextflow/24.10.2-gcc-13.2.0

export NXF_HOME=/scratch/users/${USER}/nextflow/
export NXF_CACHE=/scratch/users/${USER}/nextflow/cache
export NXF_TEMP=/scratch/users/${USER}/nextflow/tmp
export NXF_SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/
export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

# See https://github.com/nextflow-io/nextflow/issues/2695#issuecomment-1635939435
nohup lsof +D /scratch_tmp/prj/ppn_als_longread/work -r 600 &> /dev/null &
LSOF_PID=$!
trap "kill $LSOF_PID" EXIT

nextflow run epi2me-labs/wf-human-variation \
    -r v2.6.0 \
    -c /scratch_tmp/prj/ppn_als_longread/wf-human-variation/config/wf-human-variation.config \
    --bam /scratch_tmp/prj/ppn_als_longread/wf-basecalling/basecalled/A153_06_sup \
    --out_dir /scratch_tmp/prj/ppn_als_longread/wf-human-variation/vcf/A153_06_sup \
    --sample_name A153_06_sup \
    --sex XX \
    --override_basecaller_cfg dna_r10.4.1_e8.2_400bps_sup@v5.0.0

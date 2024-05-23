#!/bin/bash -l

#SBATCH --job-name=A085_00_sup_wf-basecalling
#SBATCH --partition=nd_bioinformatics_cpu,cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/jobs/A085_00_sup

module load nextflow/23.10.0-gcc-13.2.0

# Set the Singularity and Nextflow cache directories
export NXF_HOME=/scratch/users/${USER}/nextflow/
export NXF_CACHE=/scratch/users/${USER}/nextflow/cache
export NXF_TEMP=/scratch/users/${USER}/nextflow/tmp
export NXF_SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/
export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/

# Specify Nextflow max heap size
export NFX_OPTS="-Xms512M -Xmx8G"

nextflow run epi2me-labs/wf-basecalling \
    -r v1.1.7 \
    -c /scratch/prj/ppn_als_longread/config/wf-basecalling.config \
    --input /scratch/prj/ppn_als_longread/pod5/A085_00 \
    --out_dir /scratch/prj/ppn_als_longread/basecalled/A085_00_sup \
    --sample_name A085_00_sup \
    --basecaller_cfg dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
    --remora_cfg dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mCG_5hmCG@v1

#!/bin/bash -l

#SBATCH --job-name=A154_04_hac_wf-human-variation
#SBATCH --partition=nd_bioinformatics_cpu,cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/jobs/A154_04_hac

module load nextflow/22.10.1-gcc-13.2.0

# Set the Singularity and Nextflow cache directories
export NXF_HOME=/scratch/users/${USER}/nextflow/
export NXF_CACHE=/scratch/users/${USER}/nextflow/cache
export NXF_TEMP=/scratch/users/${USER}/nextflow/tmp
export NXF_SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/
export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/

# Specify Nextflow max heap size
export NFX_OPTS="-Xms512M -Xmx8G"

nextflow run epi2me-labs/wf-human-variation \
    -r v2.1.0 \
    -c /scratch/prj/ppn_als_longread/config/wf-human-variation.config \
    --bam /scratch/prj/ppn_als_longread/basecalled/A154_04_hac \
    --out_dir /scratch/prj/ppn_als_longread/vcf/A154_04_hac \
    --sample_name A154_04_hac \
    --sex male \
    --basecaller_cfg dna_r10.4.1_e8.2_400bps_hac@v4.3.0

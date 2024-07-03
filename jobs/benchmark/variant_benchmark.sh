#!/bin/bash -l

#SBATCH --job-name=ont_benchmark
#SBATCH --partition=nd_bioinformatics_cpu,cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/

module load nextflow/23.10.0-gcc-13.2.0

# Set the Singularity and Nextflow cache directories
export NXF_HOME=/scratch/users/${USER}/nextflow/
export NXF_CACHE=/scratch/users/${USER}/nextflow/cache
export NXF_TEMP=/scratch/users/${USER}/nextflow/tmp
export NXF_SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/
export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/

# Specify Nextflow max heap size
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

nextflow run main.nf

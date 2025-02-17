#!/bin/bash -l

#SBATCH --job-name=ont_benchmark
#SBATCH --partition=nd_bioinformatics_cpu,cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=0-04:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/ont-benchmark

module load nextflow/24.10.2-gcc-13.2.0

export NXF_HOME=/scratch/users/${USER}/nextflow/
export NXF_CACHE=/scratch/users/${USER}/nextflow/cache
export NXF_TEMP=/scratch/users/${USER}/nextflow/tmp
export NXF_SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/
export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

# See https://github.com/nextflow-io/nextflow/issues/2695#issuecomment-1635939435
nohup lsof +D /scratch/prj/ppn_als_longread/ont-benchmark/work -r 600 &> /dev/null &
LSOF_PID=$!
trap "kill $LSOF_PID" EXIT

nextflow run main.nf

#!/bin/bash -l

#SBATCH --job-name=ops-jupyter
#SBATCH --partition=nd_bioinformatics_cpu,cpu,interruptible_cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-02:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/jupyter

# Load Anaconda module and activate the environment
module load anaconda3/2022.10-gcc-13.2.0
source activate jupyter

jupyter nbconvert --to html --execute ont_qc.ipynb

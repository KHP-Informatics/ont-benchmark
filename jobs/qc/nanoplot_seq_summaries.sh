#!/bin/bash -l

#SBATCH --job-name=nanoplot
#SBATCH --partition=nd_bioinformatics_cpu,cpu
#SBATCH --ntasks=11
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-01:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/jobs/qc

export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/

for summary_file in /scratch/prj/ppn_als_longread/seq_summaries/*/*.txt; do
    if [ -f "$summary_file" ]; then
        subdir_name=$(basename "$(dirname "$summary_file")")
        outdir="/scratch/prj/ppn_als_longread/qc/nanoplot/seq_summaries/${subdir_name}"

        # docker://staphb/nanoplot:1.42.0
        srun -n1 singularity exec --bind /scratch:/scratch /scratch/users/k2474617/containers/nanoplot_1.42.0.sif NanoPlot \
            --threads 4 \
            --outdir "$outdir" \
            --color royalblue \
            --colormap Viridis \
            --plots kde \
            --dpi 300 \
            --summary "$summary_file" &
    fi
done

wait

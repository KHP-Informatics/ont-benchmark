#!/bin/bash -l

#SBATCH --job-name=nanoplot_seq_summaries
#SBATCH --partition=nd_bioinformatics_cpu,cpu,interruptible_cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-01:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/jobs/qc

export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/

summary_files=()
summary_names=()

for summary_file in /scratch/prj/ppn_als_longread/seq_summaries/*/*.txt; do
    if [ -f "${summary_file}" ]; then
        subdir_name=$(basename "$(dirname "${summary_file}")")
        nanoplot_outdir="/scratch/prj/ppn_als_longread/qc/nanoplot/seq_summaries/${subdir_name}"

        summary_files+=("${summary_file}")
        summary_names+=("${subdir_name}")

        if [[ "${subdir_name}" == *__* ]]; then
            srun -n1 -N1 singularity exec --bind /scratch:/scratch docker://quay.io/biocontainers/nanoplot:1.42.0--pyhdfd78af_0 \
                NanoPlot \
                    --summary "${summary_file}" \
                    --outdir "${nanoplot_outdir}" \
                    --N50 \
                    --loglength \
                    --color royalblue \
                    --colormap Viridis \
                    --plots kde \
                    --dpi 300 \
                    --tsv_stats \
                    --info_in_report \
                    --threads 2 \
                    --barcoded &
        else
            srun -n1 -N1 singularity exec --bind /scratch:/scratch docker://quay.io/biocontainers/nanoplot:1.42.0--pyhdfd78af_0 \
                NanoPlot \
                    --summary "${summary_file}" \
                    --outdir "${nanoplot_outdir}" \
                    --N50 \
                    --loglength \
                    --color royalblue \
                    --colormap Viridis \
                    --plots kde \
                    --dpi 300 \
                    --tsv_stats \
                    --info_in_report \
                    --threads 2 &
        fi
    fi
done

nanocomp_outdir="/scratch/prj/ppn_als_longread/qc/nanocomp/seq_summaries/"

# Function to join array elements with a space
join_by_space() { local IFS=' '; echo "$*"; }

if [ ${#summary_files[@]} -gt 0 ]; then
    summary_files_str=$(join_by_space "${summary_files[@]}")
    summary_names_str=$(join_by_space "${summary_names[@]}")
    srun -n1 -N1 singularity exec --bind /scratch:/scratch docker://quay.io/biocontainers/nanocomp:1.23.1--pyhdfd78af_0 \
        NanoComp \
            --summary ${summary_files_str} \
            --outdir "${nanocomp_outdir}" \
            --plot violin \
            --names ${summary_names_str} \
            --threads 2 &
fi

wait

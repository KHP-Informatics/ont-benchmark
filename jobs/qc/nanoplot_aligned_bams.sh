#!/bin/bash -l

#SBATCH --job-name=nanoplot_aligned_bams
#SBATCH --partition=nd_bioinformatics_cpu,cpu,interruptible_cpu
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=512M
#SBATCH --time=0-06:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/jobs/qc

export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/

pickle_files=()
pickle_names=()
pickle_files_hac=()
pickle_names_hac=()
pickle_files_sup=()
pickle_names_sup=()

for cram_file in /scratch/prj/ppn_als_longread/vcf/*/*.haplotagged.cram; do
    if [ -f "${cram_file}" ]; then
        subdir_name=$(basename "$(dirname "${cram_file}")")
        subdir_name=${subdir_name%.haplotagged}

        outdir="/scratch/prj/ppn_als_longread/qc/nanoplot/aligned_bams/${subdir_name}"

        pickle_file="${outdir}/NanoPlot-data.pickle"

        if [ -f "${pickle_file}" ]; then
            input_option="--pickle ${pickle_file}"
        else
            input_option="--cram ${cram_file}"
        fi

        pickle_files+=("${pickle_file}")
        pickle_names+=("${subdir_name}")

        if [[ "${subdir_name}" == *"_hac" ]]; then
            pickle_files_hac+=("${pickle_file}")
            pickle_names_hac+=("${subdir_name}")
        elif [[ "${subdir_name}" == *"_sup" ]]; then
            pickle_files_sup+=("${pickle_file}")
            pickle_names_sup+=("${subdir_name}")
        fi

        srun -n1 -N1 singularity exec --bind /scratch:/scratch docker://quay.io/biocontainers/nanoplot:1.42.0--pyhdfd78af_0 \
            NanoPlot \
                ${input_option} \
                --outdir "${outdir}" \
                --N50 \
                --loglength \
                --color royalblue \
                --colormap Viridis \
                --plots kde \
                --dpi 300 \
                --tsv_stats \
                --info_in_report \
                --store \
                --threads 8 &
    fi
done

wait

nanocomp_outdir="/scratch/prj/ppn_als_longread/qc/nanocomp/aligned_bams/"

# Function to join array elements with a space
join_by_space() { local IFS=' '; echo "$*"; }

run_nanocomp() {
    local outdir_suffix=$1
    shift
    local pickle_files_arr=("$@")
    local pickle_names_arr=("${pickle_files_arr[@]:${#pickle_files_arr[@]}/2}")
    pickle_files_arr=("${pickle_files_arr[@]:0:${#pickle_files_arr[@]}/2}")

    local pickle_files_str=$(join_by_space "${pickle_files_arr[@]}")
    local pickle_names_str=$(join_by_space "${pickle_names_arr[@]}")

    local outdir="/scratch/prj/ppn_als_longread/qc/nanocomp/aligned_bams${outdir_suffix}"

    echo "pickle_files_arr: ${pickle_files_arr[@]}"
    echo "pickle_files_str: ${pickle_files_str}"
    echo "pickle_names_arr: ${pickle_names_arr[@]}"
    echo "pickle_names_str: ${pickle_names_str}"

    if [ ${#pickle_files_arr[@]} -gt 0 ]; then
        srun -n1 -N1 singularity exec --bind /scratch:/scratch docker://quay.io/biocontainers/nanocomp:1.23.1--pyhdfd78af_0 \
            NanoComp \
                --pickle ${pickle_files_str} \
                --outdir "${outdir}" \
                --plot violin \
                --names ${pickle_names_str} \
                --threads 8 &
    fi
}

run_nanocomp "" "${pickle_files[@]}" "${pickle_names[@]}"

run_nanocomp "_hac" "${pickle_files_hac[@]}" "${pickle_names_hac[@]}"

run_nanocomp "_sup" "${pickle_files_sup[@]}" "${pickle_names_sup[@]}"

wait

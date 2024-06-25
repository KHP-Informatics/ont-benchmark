#!/bin/bash -l

#SBATCH --job-name=variant_benchmark
#SBATCH --partition=nd_bioinformatics_cpu,cpu,interruptible_cpu
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-01:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/jobs/qc

export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/

bcftools_container="/scratch/users/k2474617/containers/bcftools_1.20.sif"
rtg_tools_container="/scratch/users/k2474617/containers/rtg-tools_3.12.1--hdfd78af_0.sif"
survivor_container="/scratch/users/k2474617/containers/survivor_1.0.7--hdcf5f25_4.sif"

vcf_dir_ont="/scratch/prj/ppn_als_longread/vcf"
output_dir_ont="/scratch/prj/ppn_als_longread/benchmark_data/ont/snp_indel"
vcf_dir_illumina="/scratch/prj/ppn_als_longread/benchmark_data/illumina/grch38"
output_dir_illumina="/scratch/prj/ppn_als_longread/benchmark_data/illumina/snp_indel"
reference_fasta="/scratch/prj/ppn_als_longread/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
reference_base=$(basename ${reference_fasta} .fna)
reference_sdf="/scratch/prj/ppn_als_longread/references/${reference_base}.sdf"
sample_id_file="/scratch/prj/ppn_als_longread/sample_ids.csv"
eval_output_base_dir="/scratch/prj/ppn_als_longread/qc/variant_benchmark"
roc_output_base_dir="/scratch/prj/ppn_als_longread/qc/roc_plots"
sv_output_dir_ont="/scratch/prj/ppn_als_longread/benchmark_data/ont/sv"
sv_output_dir_illumina="/scratch/prj/ppn_als_longread/benchmark_data/illumina/sv"

create_sdf_if_missing() {
    if [ ! -d "${reference_sdf}" ]; then
        echo "Reference SDF not found. Creating SDF..."
        srun -n1 -N1 singularity exec --bind /scratch:/scratch ${rtg_tools_container} rtg format -o ${reference_sdf} ${reference_fasta}
    else
        echo "Reference SDF already exists."
    fi
}

create_directories() {
    mkdir -p ${output_dir_ont}
    mkdir -p ${output_dir_illumina}
    mkdir -p ${eval_output_base_dir}
    mkdir -p ${roc_output_base_dir}
    mkdir -p ${sv_output_dir_ont}
    mkdir -p ${sv_output_dir_illumina}
}

declare -A sample_id_map
read_sample_id_mapping() {
    while IFS=',' read -r ont_id illumina_id; do
        ont_id=$(echo $ont_id | sed 's/_hac$//;s/_sup$//' | tr -d '[:space:]')
        illumina_id=$(echo $illumina_id | tr -d '[:space:]')
        if [[ -n "${ont_id}" && -n "${illumina_id}" ]]; then
            sample_id_map[$ont_id]=$illumina_id
        fi
    done < <(tail -n +2 ${sample_id_file})
}

process_ont_vcf_files() {
    local vcf_dir=$1
    local output_dir=$2
    local filter_expression=""

    shopt -s nullglob
    for vcf_file in "${vcf_dir}"/*/*.wf_snp.vcf.gz; do
        sample_id=$(basename "$(dirname "${vcf_file}")")
        echo "Processing ONT VCF file: ${vcf_file}, Sample ID: ${sample_id}"
        ont_id=$(echo "${sample_id}" | sed 's/_hac$//;s/_sup$//' | tr -d '[:space:]')
        illumina_id="${sample_id_map[${ont_id}]}"
        echo "ONT ID: ${ont_id}, Illumina ID: ${illumina_id}"

        if [[ -z "${illumina_id}" ]]; then
            echo "No matching Illumina sample ID found for ONT sample ID ${ont_id}"
            continue
        fi

        output_snps="${output_dir}/${sample_id}_snps.vcf.gz"
        output_indels="${output_dir}/${sample_id}_indels.vcf.gz"

        echo "Processing ONT VCF for sample ID: ${sample_id}"
        srun -n1 -N1 singularity exec --bind /scratch:/scratch ${bcftools_container} bash -c "
            bcftools view -v snps,mnps ${vcf_file} ${filter_expression} -Oz -o ${output_snps}
            bcftools index -t -f ${output_snps}
            bcftools view -v indels ${vcf_file} ${filter_expression} -Oz -o ${output_indels}
            bcftools index -t -f ${output_indels}
        " &

    done
    shopt -u nullglob
}

process_illumina_vcf_files() {
    local vcf_dir=$1
    local output_dir=$2
    local filter_expression="-i 'QUAL>=30'"

    shopt -s nullglob
    for vcf_file in "${vcf_dir}"/*.vcf.gz; do
        sample_id=$(basename "${vcf_file}" .vcf.gz)

        if [[ "${vcf_file}" == *".SV.vcf.gz" || "${vcf_file}" == *".genome.vcf.gz" ]]; then
            continue
        fi

        echo "Processing Illumina VCF file: ${vcf_file}, Sample ID: ${sample_id}"

        output_snps="${output_dir}/${sample_id}_snps.vcf.gz"
        output_indels="${output_dir}/${sample_id}_indels.vcf.gz"

        echo "Filtering Illumina VCF for sample ID: ${sample_id}"
        srun -n1 -N1 singularity exec --bind /scratch:/scratch ${bcftools_container} bash -c "
            bcftools view -v snps,mnps ${vcf_file} ${filter_expression} -Oz -o ${output_snps}
            bcftools index -t -f ${output_snps}
            bcftools view -v indels ${vcf_file} ${filter_expression} -Oz -o ${output_indels}
            bcftools index -t -f ${output_indels}
        " &

    done
    shopt -u nullglob
}

benchmark_variants() {
    local input_dir=$1
    local baseline_dir=$2
    local eval_output_suffix=$3
    local roc_suffix=$4
    local roc_array_name=$5

    for input_vcf in ${input_dir}/*_${eval_output_suffix}.vcf.gz; do
        sample_id=$(basename ${input_vcf} _${eval_output_suffix}.vcf.gz)
        ont_id=$(echo ${sample_id} | sed 's/_hac$//;s/_sup$//' | tr -d '[:space:]')
        illumina_id=${sample_id_map[${ont_id}]}

        if [[ -z "${illumina_id}" ]]; then
            echo "No matching Illumina sample ID found for ONT sample ID ${ont_id}"
            continue
        fi

        illumina_suffix=$(echo ${eval_output_suffix} | sed 's/^.*_//')

        baseline_vcf="${baseline_dir}/${illumina_id}_${illumina_suffix}.vcf.gz"
        eval_output_dir="${eval_output_base_dir}/${sample_id}_${eval_output_suffix}"

        if [ ! -f "${baseline_vcf}" ]; then
            echo "The baseline file ${baseline_vcf} does not exist."
            continue
        fi

        echo "Running rtg vcfeval for sample ID: ${sample_id} (${eval_output_suffix})"
        srun -n1 -N1 singularity exec --bind /scratch:/scratch ${rtg_tools_container} rtg vcfeval \
            --baseline=${baseline_vcf} \
            --calls=${input_vcf} \
            --template=${reference_sdf} \
            --output=${eval_output_dir} \
            --output-mode=roc-only &

        if [[ ${eval_output_suffix} == hac* ]]; then
            eval "${roc_array_name}+=(\"${eval_output_dir}/${roc_suffix}\")"
        elif [[ ${eval_output_suffix} == sup* ]]; then
            eval "${roc_array_name}+=(\"${eval_output_dir}/${roc_suffix}\")"
        fi
    done
}

process_sv_files() {
    local vcf_dir=$1
    local output_dir=$2

    shopt -s nullglob
    for sv_file in "${vcf_dir}"/*/*.wf_sv.vcf.gz; do
        str_file=$(echo "${sv_file}" | sed 's/.wf_sv.vcf.gz/.wf_str.vcf.gz/')
        sample_id=$(basename "$(dirname "${sv_file}")")
        merged_sv_str_file="${output_dir}/${sample_id}_merged_sv_str.vcf.gz"

        if [ -f "${str_file}" ]; then
            echo "Merging SV and STR files for sample ID: ${sample_id}"
            srun -n1 -N1 singularity exec --bind /scratch:/scratch ${bcftools_container} bash -c "
                bcftools concat -a -Oz -o ${merged_sv_str_file} ${sv_file} ${str_file}
                bcftools index -t -f ${merged_sv_str_file}
            " &
        else
            echo "STR file not found for sample ID: ${sample_id}, using only SV file."
            cp ${sv_file} ${merged_sv_str_file}
            srun -n1 -N1 singularity exec --bind /scratch:/scratch ${bcftools_container} bash -c "
                bcftools index -t -f ${merged_sv_str_file}
            " &
        fi
    done

    for vcf_file in "${vcf_dir}"/*.SV.vcf.gz; do
        sample_id=$(basename "${vcf_file}" .SV.vcf.gz)
        echo "Processing Illumina SV file: ${vcf_file}, Sample ID: ${sample_id}"

        output_sv="${output_dir}/${sample_id}_sv.vcf.gz"

        cp ${vcf_file} ${output_sv}
        srun -n1 -N1 singularity exec --bind /scratch:/scratch ${bcftools_container} bash -c "
            bcftools index -t -f ${output_sv}
        " &
    done
    shopt -u nullglob
}

merge_sv_files_with_survivor() {
    local ont_sv_dir=$1
    local illumina_sv_dir=$2

    shopt -s nullglob

    for ont_sv_file in ${ont_sv_dir}/*_merged_sv_str.vcf.gz; do
        sample_id=$(basename ${ont_sv_file} _merged_sv_str.vcf.gz)
        ont_id=$(echo ${sample_id} | sed 's/_hac$//;s/_sup$//' | tr -d '[:space:]')
        illumina_id=${sample_id_map[${ont_id}]}

        if [[ -z "${illumina_id}" ]]; then
            echo "No matching Illumina sample ID found for ONT sample ID ${ont_id}"
            continue
        fi

        illumina_sv_file="${illumina_sv_dir}/${illumina_id}_sv.vcf.gz"

        if [ ! -f "${illumina_sv_file}" ]; then
            echo "The Illumina SV file ${illumina_sv_file} does not exist."
            continue
        fi

        sample_files_list="${eval_output_base_dir}/${sample_id}_sample_files.txt"
        merged_sv_file="${eval_output_base_dir}/${sample_id}_merged_sv.vcf"

        echo -e "${ont_sv_file}\n${illumina_sv_file}" > ${sample_files_list}

        echo "Merging ONT and Illumina SV files for sample ID: ${sample_id}"
        srun -n1 -N1 singularity exec --bind /scratch:/scratch ${survivor_container} bash -c "
            SURVIVOR merge ${sample_files_list} 1000 2 1 1 0 30 ${merged_sv_file}
        " &
    done
    shopt -u nullglob
}

generate_rocplot() {
    local input_files=("$@")
    local output_file="${input_files[-1]}"
    unset 'input_files[-1]'

    local valid_input_files=()
    for file in "${input_files[@]}"; do
        if [ -f "${file}.gz" ]; then
            valid_input_files+=("${file}.gz")
        else
            echo "Warning: ${file}.gz does not exist and will be excluded from the ROC plot."
        fi
    done

    if [ ${#valid_input_files[@]} -gt 0 ]; then
        echo "Generating ROC plot: ${output_file}"
        srun -n1 -N1 singularity exec --bind /scratch:/scratch ${rtg_tools_container} rtg rocplot \
            --png=${output_file} \
            --plain \
            --palette=classic \
            "${valid_input_files[@]}" &
    else
        echo "No valid input files for ROC plot: ${output_file}. Skipping."
    fi
}

main() {
    create_sdf_if_missing
    create_directories
    read_sample_id_mapping

    process_ont_vcf_files "${vcf_dir_ont}" "${output_dir_ont}"
    process_illumina_vcf_files "${vcf_dir_illumina}" "${output_dir_illumina}"
    process_sv_files "${vcf_dir_ont}" "${sv_output_dir_ont}"
    process_sv_files "${vcf_dir_illumina}" "${sv_output_dir_illumina}"

    wait

    declare -a hac_snps_rocs
    declare -a sup_snps_rocs
    declare -a hac_indels_rocs
    declare -a sup_indels_rocs

    benchmark_variants ${output_dir_ont} ${output_dir_illumina} "hac_snps" "snp_roc.tsv" "hac_snps_rocs"
    benchmark_variants ${output_dir_ont} ${output_dir_illumina} "sup_snps" "snp_roc.tsv" "sup_snps_rocs"
    benchmark_variants ${output_dir_ont} ${output_dir_illumina} "hac_indels" "non_snp_roc.tsv" "hac_indels_rocs"
    benchmark_variants ${output_dir_ont} ${output_dir_illumina} "sup_indels" "non_snp_roc.tsv" "sup_indels_rocs"

    wait

    merge_sv_files_with_survivor ${sv_output_dir_ont} ${sv_output_dir_illumina}

    wait

    generate_rocplot "${hac_snps_rocs[@]}" "${sup_snps_rocs[@]}" "${roc_output_base_dir}/all_hac_sup_snps.png"
    generate_rocplot "${hac_indels_rocs[@]}" "${sup_indels_rocs[@]}" "${roc_output_base_dir}/all_hac_sup_indels.png"
    generate_rocplot "${hac_snps_rocs[@]}" "${roc_output_base_dir}/hac_snps.png"
    generate_rocplot "${hac_indels_rocs[@]}" "${roc_output_base_dir}/hac_indels.png"
    generate_rocplot "${sup_snps_rocs[@]}" "${roc_output_base_dir}/sup_snps.png"
    generate_rocplot "${sup_indels_rocs[@]}" "${roc_output_base_dir}/sup_indels.png"

    wait

    echo "All tasks completed."
}

main

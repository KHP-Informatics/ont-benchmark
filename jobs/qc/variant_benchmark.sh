#!/bin/bash -l

#SBATCH --job-name=variant_benchmark
#SBATCH --partition=nd_bioinformatics_cpu,cpu,interruptible_cpu
#SBATCH --ntasks=42
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-01:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/jobs/qc

# Set Singularity cache directory
export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/

bcftools_container="/scratch/users/k2474617/containers/bcftools_1.20.sif"
rtg_tools_container="/scratch/users/k2474617/containers/rtg-tools_3.12.1--hdfd78af_0.sif"

vcf_dir_ont="/scratch/prj/ppn_als_longread/vcf"
output_dir_ont="/scratch/prj/ppn_als_longread/benchmark_data/ont"
vcf_dir_illumina="/scratch/prj/ppn_als_longread/benchmark_data/illumina/grch38"
output_dir_illumina="/scratch/prj/ppn_als_longread/benchmark_data/illumina/filtered"
reference_fasta="/scratch/prj/ppn_als_longread/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
reference_base=$(basename ${reference_fasta} .fna)
reference_sdf="/scratch/prj/ppn_als_longread/references/${reference_base}.sdf"
sample_id_file="/scratch/prj/ppn_als_longread/sample_ids.csv"
eval_output_base_dir="/scratch/prj/ppn_als_longread/qc/variant_benchmark"
roc_output_base_dir="/scratch/prj/ppn_als_longread/qc/roc_plots"

# Check and create SDF
if [ ! -d "${reference_sdf}" ]; then
    echo "Reference SDF not found. Creating SDF..."
    srun -n1 -N1 singularity exec --bind /scratch:/scratch ${rtg_tools_container} rtg format -o ${reference_sdf} ${reference_fasta}
else
    echo "Reference SDF already exists."
fi

# Create output directories
mkdir -p ${output_dir_ont}
mkdir -p ${output_dir_illumina}
mkdir -p ${eval_output_base_dir}
mkdir -p ${roc_output_base_dir}

# Read sample IDs mapping
declare -A sample_id_map
while IFS=',' read -r ont_id illumina_id; do
    ont_id=$(echo $ont_id | sed 's/_hac$//;s/_sup$//' | tr -d '[:space:]')
    illumina_id=$(echo $illumina_id | tr -d '[:space:]')
    sample_id_map[$ont_id]=$illumina_id
done < <(tail -n +2 ${sample_id_file})

# Process ONT VCF files
for vcf_file in ${vcf_dir_ont}/*/*.wf_snp.vcf.gz; do
    sample_id=$(basename $(dirname ${vcf_file}))
    ont_id=$(echo ${sample_id} | sed 's/_hac$//;s/_sup$//' | tr -d '[:space:]')
    illumina_id=${sample_id_map[${ont_id}]}

    if [[ -z "${illumina_id}" ]]; then
        echo "No matching Illumina sample ID found for ONT sample ID ${ont_id}"
        continue
    fi

    output_snps="${output_dir_ont}/${sample_id}_snps.vcf.gz"
    output_indels="${output_dir_ont}/${sample_id}_indels.vcf.gz"

    echo "Processing ONT VCF for sample ID: ${sample_id}"
    srun -n1 -N1 singularity exec --bind /scratch:/scratch ${bcftools_container} bash -c "
        bcftools view -v snps,mnps ${vcf_file} -Oz -o ${output_snps}
        bcftools index -t -f ${output_snps}
        bcftools view -v indels ${vcf_file} -Oz -o ${output_indels}
        bcftools index -t -f ${output_indels}
    " &
done

# Process Illumina VCF files
for vcf_file in ${vcf_dir_illumina}/*.vcf.gz; do
    base_name=$(basename ${vcf_file})
    if [[ ${base_name} != *.SV.vcf.gz && ${base_name} != *.genome.vcf.gz ]]; then
        sample_id="${base_name%.vcf.gz}"
        output_snps="${output_dir_illumina}/${sample_id}_snps.vcf.gz"
        output_indels="${output_dir_illumina}/${sample_id}_indels.vcf.gz"

        echo "Processing Illumina VCF for sample ID: ${sample_id}"
        srun -n1 -N1 singularity exec --bind /scratch:/scratch ${bcftools_container} bash -c "
            bcftools view -v snps,mnps ${vcf_file} -Oz -o ${output_snps}
            bcftools index -t -f ${output_snps}
            bcftools view -v indels -i 'QUAL>=30' ${vcf_file} -Oz -o ${output_indels}
            bcftools index -t -f ${output_indels}
        " &
    fi
done

wait

# Prepare lists for ROC plots
hac_snps_rocs=()
sup_snps_rocs=()
hac_indels_rocs=()
sup_indels_rocs=()

# Perform variant benchmarking with rtg vcfeval
for ont_snps in ${output_dir_ont}/*_snps.vcf.gz; do
    sample_id=$(basename ${ont_snps} _snps.vcf.gz)
    ont_id=$(echo ${sample_id} | sed 's/_hac$//;s/_sup$//' | tr -d '[:space:]')
    illumina_id=${sample_id_map[${ont_id}]}

    if [[ -z "${illumina_id}" ]]; then
        echo "No matching Illumina sample ID found for ONT sample ID ${ont_id}"
        continue
    fi

    illumina_snps="${output_dir_illumina}/${illumina_id}_snps.vcf.gz"
    eval_output_dir_snps="${eval_output_base_dir}/${sample_id}_snps"

    if [ ! -f "${illumina_snps}" ]; then
        echo "The baseline file ${illumina_snps} does not exist."
        continue
    fi

    echo "Running rtg vcfeval for sample ID: ${sample_id} (SNPs)"
    srun -n1 -N1 singularity exec --bind /scratch:/scratch ${rtg_tools_container} rtg vcfeval \
        --baseline=${illumina_snps} \
        --calls=${ont_snps} \
        --template=${reference_sdf} \
        --output=${eval_output_dir_snps} \
        --output-mode=roc-only &

    if [[ ${sample_id} == *_hac ]]; then
        hac_snps_rocs+=("${eval_output_dir_snps}/snp_roc.tsv")
    elif [[ ${sample_id} == *_sup ]]; then
        sup_snps_rocs+=("${eval_output_dir_snps}/snp_roc.tsv")
    fi
done

for ont_indels in ${output_dir_ont}/*_indels.vcf.gz; do
    sample_id=$(basename ${ont_indels} _indels.vcf.gz)
    ont_id=$(echo ${sample_id} | sed 's/_hac$//;s/_sup$//' | tr -d '[:space:]')
    illumina_id=${sample_id_map[${ont_id}]}

    if [[ -z "${illumina_id}" ]]; then
        echo "No matching Illumina sample ID found for ONT sample ID ${ont_id}"
        continue
    fi

    illumina_indels="${output_dir_illumina}/${illumina_id}_indels.vcf.gz"
    eval_output_dir_indels="${eval_output_base_dir}/${sample_id}_indels"

    if [ ! -f "${illumina_indels}" ]; then
        echo "The baseline file ${illumina_indels} does not exist."
        continue
    fi

    echo "Running rtg vcfeval for sample ID: ${sample_id} (Indels)"
    srun -n1 -N1 singularity exec --bind /scratch:/scratch ${rtg_tools_container} rtg vcfeval \
        --baseline=${illumina_indels} \
        --calls=${ont_indels} \
        --template=${reference_sdf} \
        --output=${eval_output_dir_indels} \
        --output-mode=roc-only &

    if [[ ${sample_id} == *_hac ]]; then
        hac_indels_rocs+=("${eval_output_dir_indels}/non_snp_roc.tsv")
    elif [[ ${sample_id} == *_sup ]]; then
        sup_indels_rocs+=("${eval_output_dir_indels}/non_snp_roc.tsv")
    fi
done

wait

# Generate ROC plots with rtg rocplot
generate_rocplot() {
    local input_files=("$@")
    local output_file="${input_files[-1]}"
    unset 'input_files[-1]'

    # Filter out non-existent .tsv.gz files
    local valid_input_files=()
    for file in "${input_files[@]}"; do
        if [ -f "${file}.gz" ]; then
            valid_input_files+=("${file}.gz")
        else
            echo "Warning: ${file}.gz does not exist and will be excluded from the ROC plot."
        fi
    done

    # Only run the plot if there are valid input files
    if [ ${#valid_input_files[@]} -gt 0 ]; then
        echo "Generating ROC plot: ${output_file}"
        srun -n1 -N1 singularity exec --bind /scratch:/scratch ${rtg_tools_container} rtg rocplot \
            --png=${output_file} \
            --plain \
            --palette=classic \
            ${valid_input_files[@]} &
    else
        echo "No valid input files for ROC plot: ${output_file}. Skipping."
    fi
}


# 1. All hac and sup snps
generate_rocplot "${hac_snps_rocs[@]}" "${sup_snps_rocs[@]}" "${roc_output_base_dir}/all_hac_sup_snps.png"

# 2. All hac and sup indels
generate_rocplot "${hac_indels_rocs[@]}" "${sup_indels_rocs[@]}" "${roc_output_base_dir}/all_hac_sup_indels.png"

# 3. Only hac snps
generate_rocplot "${hac_snps_rocs[@]}" "${roc_output_base_dir}/hac_snps.png"

# 4. Only hac indels
generate_rocplot "${hac_indels_rocs[@]}" "${roc_output_base_dir}/hac_indels.png"

# 5. Only sup snps
generate_rocplot "${sup_snps_rocs[@]}" "${roc_output_base_dir}/sup_snps.png"

# 6. Only sup indels
generate_rocplot "${sup_indels_rocs[@]}" "${roc_output_base_dir}/sup_indels.png"

wait

echo "All tasks completed."

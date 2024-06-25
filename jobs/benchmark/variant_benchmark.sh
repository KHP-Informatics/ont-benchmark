#!/bin/bash -l

#SBATCH --job-name=variant_benchmark
#SBATCH --partition=nd_bioinformatics_cpu,cpu,interruptible_cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-01:00:00
#SBATCH --output=/scratch/prj/ppn_als_longread/slurm_jobs/%j_%x.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renato.santos@kcl.ac.uk
#SBATCH --chdir /scratch/prj/ppn_als_longread/jobs/qc

export SINGULARITY_CACHEDIR=/scratch/users/${USER}/singularity/

bcftools_container="docker://quay.io/biocontainers/bcftools:1.20--h8b25389_0"
rtg_tools_container="docker://quay.io/biocontainers/rtg-tools:3.12.1--hdfd78af_0"
curl_jq_container="docker://badouralix/curl-jq:ubuntu"

reference_fasta="/scratch/prj/ppn_als_longread/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
reference_base=$(basename ${reference_fasta} .fna)
reference_sdf="/scratch/prj/ppn_als_longread/references/${reference_base}.sdf"

sample_id_file="/scratch/prj/ppn_als_longread/sample_ids.csv"
positions_file="/scratch/prj/ppn_als_longread/variant_benchmark/microarray_positions.csv"

microarray_dir="/scratch/prj/ppn_als_longread/benchmark_data/genotyping"
processed_microarray_dir="/scratch/prj/ppn_als_longread/benchmark_data/genotyping/grch38"
ont_dir="/scratch/prj/ppn_als_longread/vcf"
variant_benchmark_dir="/scratch/prj/ppn_als_longread/variant_benchmark/snps_ont_vs_microarray"

create_sdf_if_missing() {
    if [ ! -d "${reference_sdf}" ]; then
        echo "$(date) - Reference SDF not found. Creating SDF..."
        srun -n1 -N1 singularity exec --bind /scratch:/scratch ${rtg_tools_container} rtg format -o ${reference_sdf} ${reference_fasta}
        if [ $? -ne 0 ]; then
            echo "$(date) - Error creating SDF from reference FASTA"
            exit 1
        fi
    else
        echo "$(date) - Reference SDF already exists."
    fi
}

create_directories() {
    mkdir -p ${processed_microarray_dir}
    mkdir -p ${variant_benchmark_dir}
    mkdir -p $(dirname ${positions_file})
}

declare -A sample_id_map
read_sample_id_mapping() {
    while IFS=',' read -r ont_id LP_id; do
        ont_id=$(echo $ont_id | sed 's/_hac$//;s/_sup$//' | tr -d '[:space:]')
        LP_id=$(echo $LP_id | tr -d '[:space:]')
        if [[ -n "${ont_id}" && -n "${LP_id}" ]]; then
            sample_id_map[$ont_id]=$LP_id
        fi
    done < <(tail -n +2 ${sample_id_file})
}

fetch_variant_positions() {
    local variant_ids=("$@")

    if [ ! -f "${positions_file}" ]; then
        echo "variant_id,chromosome,position" > ${positions_file}
    fi

    for variant_id in "${variant_ids[@]}"; do
        if grep -q "^${variant_id}," ${positions_file}; then
            echo "$(date) - Variant ${variant_id} already fetched. Skipping."
            continue
        fi

        response_file=$(mktemp)
        response=$(singularity exec --bind /scratch:/scratch ${curl_jq_container} sh -c "curl -s -w \"%{http_code}\" -o ${response_file} https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/${variant_id}")
        http_code=$(echo "${response}" | tail -n 1)
        if [ "${http_code}" -ne 200 ]; then
            echo "$(date) - Error fetching variant ${variant_id} from API. HTTP code: ${http_code}"
            rm -f ${response_file}
            continue
        fi

        chromosome=$(singularity exec --bind /scratch:/scratch ${curl_jq_container} jq -r '.primary_snapshot_data.placements_with_allele[] | select(.is_ptlp == true and .placement_annot.assembly_name == "GRCh38") | .placement_annot.seq_id' ${response_file})
        position=$(singularity exec --bind /scratch:/scratch ${curl_jq_container} jq -r '.primary_snapshot_data.placements_with_allele[] | select(.is_ptlp == true and .placement_annot.assembly_name == "GRCh38") | .alleles[] | select(.allele.spdi.allele == .allele.spdi.allele).allele.spdi.position' ${response_file})

        if [[ -n "${chromosome}" && -n "${position}" ]]; then
            echo "${variant_id},${chromosome},${position}" >> ${positions_file}
        else
            echo "$(date) - Incomplete data for variant ${variant_id}. Skipping."
        fi

        rm -f ${response_file}
    done
}

preprocess_microarray_vcf() {
    local microarray_vcf=$1
    echo "$(date) - Preprocessing microarray VCF ${microarray_vcf}..."
    local variant_ids=($(zcat ${microarray_vcf} | grep -v '^#' | awk '{print $3}'))
    fetch_variant_positions "${variant_ids[@]}"

    zcat ${microarray_vcf} | awk -v OFS='\t' -F '\t' '
    BEGIN { while ((getline < "'${positions_file}'") > 0) { pos[$1] = $2 "\t" $3 } }
    /^#/ { print; next }
    {
        if ($3 in pos) {
            split(pos[$3], position, "\t")
            $1 = position[1]
            $2 = position[2]
            print
        }
    }' | bgzip > ${processed_microarray_dir}/$(basename ${microarray_vcf} .vcf.gz)_processed.vcf.gz

    if [ $? -ne 0 ]; then
        echo "$(date) - Error preprocessing microarray VCF ${microarray_vcf}"
        exit 1
    fi

    singularity exec --bind /scratch:/scratch ${bcftools_container} bcftools index ${processed_microarray_dir}/$(basename ${microarray_vcf} .vcf.gz)_processed.vcf.gz
}

compare_variants() {
    local ont_vcf=$1
    local microarray_processed_vcf=$2
    echo "$(date) - Comparing variants using RTG Tools..."
    singularity exec --bind /scratch:/scratch ${rtg_tools_container} rtg vcfeval \
        -b ${microarray_processed_vcf} \
        -c ${ont_vcf} \
        -t ${reference_sdf} \
        -o ${variant_benchmark_dir}/$(basename ${ont_vcf} .vcf.gz)_vcfeval_output

    if [ $? -ne 0 ]; then
        echo "$(date) - Error comparing variants for ${ont_vcf}"
        exit 1
    fi
}

main() {
    create_sdf_if_missing
    create_directories
    read_sample_id_mapping

    for ont_id in "${!sample_id_map[@]}"; do
        LP_id=${sample_id_map[$ont_id]}
        microarray_vcf=${microarray_dir}/FinalReport_InfiniumOmni2-5-8v1-4_${LP_id}.vcf.gz
        ont_vcf=${ont_dir}/${ont_id}_sup/${ont_id}_sup.wf_snp.vcf.gz

        if [[ -f "${microarray_vcf}" && -f "${ont_vcf}" ]]; then
            preprocess_microarray_vcf ${microarray_vcf}
            microarray_processed_vcf=${processed_microarray_dir}/$(basename ${microarray_vcf} .vcf.gz)_processed.vcf.gz
            compare_variants ${ont_vcf} ${microarray_processed_vcf}
        else
            echo "$(date) - Files for ${ont_id} (${LP_id}) are missing."
        fi
    done
}

main

# Filter out repetitive/low-complexity regions

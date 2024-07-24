#!/usr/bin/env nextflow

process RTG_VCFEVAL {
    tag "${sample_id}|${file_type}|${region_type}"

    publishDir "${params.outdir}/${variant_type}/rtg_vcfeval/${region_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(other_id),
          path(query_vcf), path(query_index),
          path(truth_vcf), path(truth_index),
          val(file_type), val(region_type), path(reference_sdf)
    val(variant_type)

    output:
    path "${rtg_dir}/allele_non_snp_roc.tsv.gz", emit: allele_non_snp_roc
    path "${rtg_dir}/allele_weighted_roc.tsv.gz", emit: allele_weighted_roc
    path "${rtg_dir}/non_snp_roc.tsv.gz", emit: non_snp_roc
    path "${rtg_dir}/summary.txt", emit: summary
    path "${rtg_dir}/${sample_id}_${variant_type}_${region_type}_weighted_roc.tsv.gz", emit: weighted_roc
    path "${rtg_dir}/allele_snp_roc.tsv.gz", emit: allele_snp_roc
    path "${rtg_dir}/phasing.txt", emit: phasing
    path "${rtg_dir}/snp_roc.tsv.gz", emit: snp_roc
    path "${rtg_dir}/vcfeval.log", emit: vcfeval_log

    script:
    rtg_dir = "${sample_id}.${variant_type}"
    """
    rtg vcfeval \
        --baseline=${truth_vcf} \
        --calls=${query_vcf} \
        --template=${reference_sdf} \
        --output=${rtg_dir} \
        --output-mode=roc-only

    cp \
        ${rtg_dir}/weighted_roc.tsv.gz \
        ${rtg_dir}/${sample_id}_${variant_type}_${region_type}_weighted_roc.tsv.gz
    """
}

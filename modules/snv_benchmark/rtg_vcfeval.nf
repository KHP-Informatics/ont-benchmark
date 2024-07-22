#!/usr/bin/env nextflow

process RTG_VCFEVAL {
    tag "${sample_id}|${file_type}"

    publishDir "${params.outdir}/snv/rtg_vcfeval/", mode: 'copy'

    input:
    tuple val(sample_id), val(other_id),
          path(query_vcf), path(query_index),
          path(truth_vcf), path(truth_index),
          val(file_type), path(reference_sdf)

    output:
    path "${rtg_dir}/allele_non_snp_roc.tsv.gz", emit: allele_non_snp_roc
    path "${rtg_dir}/allele_weighted_roc.tsv.gz", emit: allele_weighted_roc
    path "${rtg_dir}/non_snp_roc.tsv.gz", emit: non_snp_roc
    path "${rtg_dir}/summary.txt", emit: summary
    path "${rtg_dir}/weighted_roc.tsv.gz", emit: weighted_roc
    path "${rtg_dir}/allele_snp_roc.tsv.gz", emit: allele_snp_roc
    path "${rtg_dir}/phasing.txt", emit: phasing
    path "${rtg_dir}/snp_roc.tsv.gz", emit: snp_roc
    path "${rtg_dir}/vcfeval.log", emit: vcfeval_log

    script:
    rtg_dir = "${sample_id}.${file_type}"
    """
    rtg vcfeval --baseline=${truth_vcf} \
                --calls=${query_vcf} \
                --template=${reference_sdf} \
                --output=${rtg_dir} \
                --output-mode=roc-only
    """
}

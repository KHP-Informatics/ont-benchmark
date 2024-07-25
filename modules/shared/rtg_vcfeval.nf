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
    path "${rtg_dir}/${sample_id}_${variant_type}_${region_type}_weighted_roc.tsv.gz", emit: weighted_roc
    path "${rtg_dir}/summary.txt", emit: summary
    path "${rtg_dir}/vcfeval.log", emit: vcfeval_log
    path "${rtg_dir}/tp.vcf.gz", emit: tp_vcf
    path "${rtg_dir}/tp-baseline.vcf.gz", emit: tp_baseline_vcf
    path "${rtg_dir}/fp.vcf.gz", emit: fp_vcf
    path "${rtg_dir}/fn.vcf.gz", emit: fn_vcf
    path "${rtg_dir}/tp.vcf.gz.tbi", emit: tp_vcf_index
    path "${rtg_dir}/tp-baseline.vcf.gz.tbi", emit: tp_baseline_vcf_index
    path "${rtg_dir}/fp.vcf.gz.tbi", emit: fp_vcf_index
    path "${rtg_dir}/fn.vcf.gz.tbi", emit: fn_vcf_index
    path "${rtg_dir}/query.vcf.gz", emit: query_vcf
    path "${rtg_dir}/query.vcf.gz.tbi", emit: query_index
    path "${rtg_dir}/truth.vcf.gz", emit: truth_vcf
    path "${rtg_dir}/truth.vcf.gz.tbi", emit: truth_index

    script:
    rtg_dir = "${sample_id}.${variant_type}"
    """
    rtg vcfeval \
        --baseline=${truth_vcf} \
        --calls=${query_vcf} \
        --template=${reference_sdf} \
        --output=${rtg_dir} \
        --output-mode=split

    cp ${rtg_dir}/weighted_roc.tsv.gz \
        ${rtg_dir}/${sample_id}_${variant_type}_${region_type}_weighted_roc.tsv.gz
    cp ${query_vcf} ${rtg_dir}/query.vcf.gz
    cp ${query_index} ${rtg_dir}/query.vcf.gz.tbi
    cp ${truth_vcf} ${rtg_dir}/truth.vcf.gz
    cp ${truth_index} ${rtg_dir}/truth.vcf.gz.tbi
    """
}

#!/usr/bin/env nextflow

process SURVIVOR {
    tag "${ont_id}|${lp_id}"

    publishDir "${params.outdir}/sv/survivor/", mode: 'copy'

    input:
    tuple val(ont_id), val(lp_id),
        path(ont_vcf), path(illumina_vcf)

    output:
    tuple val(ont_id), val(lp_id), path("${ont_id}_${lp_id}_merged.vcf"), emit: merged_vcf

    script:
    """
    ls *.vcf > sample_files
    SURVIVOR merge sample_files 1000 1 1 1 0 30 ${ont_id}_${lp_id}_merged.vcf
    """
}

#!/usr/bin/env nextflow

process SPLIT_SNV_INDELS {
    tag "${meta.id}|${meta.type}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}_snv.vcf.gz"), emit: snv
    tuple val(meta), path("${meta.id}_indel.vcf.gz"), emit: indel

    script:
    """
    bcftools view -v snps ${vcf} -Oz -o ${meta.id}_snv.vcf.gz
    bcftools view -v mnps,indels ${vcf} -Oz -o ${meta.id}_indel.vcf.gz
    """
}

#!/usr/bin/env nextflow

process SPLIT_SNV_INDELS {
    tag "${meta.id}|${meta.type}"
    cpus 2

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.snv.vcf.gz"), emit: snv
    tuple val(meta), path("${meta.id}.indel.vcf.gz"), emit: indel

    script:
    """
    bcftools view -Ou -v snps,mnps ${vcf} | \
    bcftools norm -Ou -m-any --do-not-normalize | \
    bcftools view -Oz -o ${meta.id}.snv.vcf.gz

    bcftools view -v indels ${vcf} -Oz -o ${meta.id}.indel.vcf.gz
    """
}

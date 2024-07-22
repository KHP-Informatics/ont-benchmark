#!/usr/bin/env nextflow

process SPLIT_SNV_INDELS {
    tag "${meta.id}|${meta.type}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.snv.vcf.gz"), emit: snv
    tuple val(meta), path("${meta.id}.indel.vcf.gz"), emit: indel

    script:
    """
    bcftools norm \
        --do-not-normalize \
        -m-any ${vcf} -Ou | \
    tee >(bcftools view --types snps,mnps -Ou | \
          bcftools norm --atomize -Oz -o ${meta.id}.snv.vcf.gz) | \
    bcftools view --types indels -Oz -o ${meta.id}.indel.vcf.gz
    """
}

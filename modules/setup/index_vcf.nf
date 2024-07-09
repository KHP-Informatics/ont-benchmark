#!/usr/bin/env nextflow

process INDEX_VCF {
    tag "${meta.id}|${meta.type}|${meta.variant}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${vcf}.sorted"), path("${vcf}.sorted.tbi")

    script:
    """
    bcftools sort ${vcf} -o ${vcf}.sorted
    bcftools index --force --tbi ${vcf}.sorted
    """
}

#!/usr/bin/env nextflow

process INDEX_VCF {
    tag "${meta.id}|${meta.type}|${meta.variant}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path(vcf), path("${vcf}.tbi")

    script:
    """
    bcftools index --force --tbi ${vcf}
    """
}

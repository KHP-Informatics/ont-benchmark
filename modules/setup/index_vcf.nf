#!/usr/bin/env nextflow

process INDEX_VCF {
    tag "${meta.id}|${meta.type}|${meta.variant}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.${meta.type}.${meta.variant}.sorted.vcf.gz"), path("${meta.id}.${meta.type}.${meta.variant}.sorted.vcf.gz.tbi")

    script:
    """
    bcftools sort ${vcf} -Oz -o ${meta.id}.${meta.type}.${meta.variant}.sorted.vcf.gz
    bcftools index --force --tbi ${meta.id}.${meta.type}.${meta.variant}.sorted.vcf.gz
    """
}

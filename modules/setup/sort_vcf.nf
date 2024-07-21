#!/usr/bin/env nextflow

process SORT_VCF {
    tag "${meta.id}|${meta.type}|${meta.variant}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.${meta.type}.${meta.variant}.sorted.vcf.gz"), path("${meta.id}.${meta.type}.${meta.variant}.sorted.vcf.gz.tbi")

    script:
    def max_mem = task.memory ? "${task.memory.toMega()}M" : ''
    """
    bcftools sort \
        ${vcf} \
        --max-mem ${max_mem} \
        --write-index=tbi \
        -Oz \
        -o ${meta.id}.${meta.type}.${meta.variant}.sorted.vcf.gz
    """
}

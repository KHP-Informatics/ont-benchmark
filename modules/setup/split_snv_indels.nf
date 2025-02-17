#!/usr/bin/env nextflow

process SPLIT_SNV_INDELS {
    tag "${meta.id}|${meta.type}"

    input:
    tuple val(meta), path(vcf), path(reference_fasta)

    output:
    tuple val(meta), path("${meta.id}.snv.vcf.gz"), emit: snv
    tuple val(meta), path("${meta.id}.indel.vcf.gz"), emit: indel

    script:
    """
    bcftools norm \
        --threads ${task.cpus} \
        --fasta-ref ${reference_fasta} \
        -m-any --do-not-normalize -Ou ${vcf} | \
    tee >(bcftools view --threads ${task.cpus} --types snps,mnps -Ou | \
          bcftools norm --threads ${task.cpus} --atomize -Oz -o ${meta.id}.snv.vcf.gz) | \
    bcftools view --threads ${task.cpus} --types indels -Oz -o ${meta.id}.indel.vcf.gz
    """
}

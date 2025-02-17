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
    bcftools index --threads ${task.cpus} ${vcf}

    bcftools view \
        --threads ${task.cpus} \
        --regions '^chr[1-9,X,Y]\\|^chr[1,2][0-9]\\|^[1-9]\\|^[1,2][0-9]\\|^[X,Y]' \
        ${vcf} | \
    bcftools norm \
        --threads ${task.cpus} \
        --fasta-ref ${reference_fasta} \
        -m-any -Ou | \
    tee >(bcftools view --threads ${task.cpus} --types snps,mnps -Ou | \
          bcftools norm --threads ${task.cpus} --atomize -Oz -o ${meta.id}.snv.vcf.gz) | \
    bcftools view --threads ${task.cpus} --types indels -Oz -o ${meta.id}.indel.vcf.gz
    """
}

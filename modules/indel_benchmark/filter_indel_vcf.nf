#!/usr/bin/env nextflow

process FILTER_INDEL_VCF {
    tag "${ont_id}|${lp_id}|${size != 'aggregate' ? size + 'bp' : 'aggregate'}"

    input:
    tuple val(ont_id), val(lp_id),
        path(ont_vcf), path(ont_index),
        path(illumina_vcf), path(illumina_index),
        path(high_confidence_regions),
        val(size)

    output:
    tuple val(ont_id), val(lp_id),
        path("${ont_id}.ont.indel.${size}.hc.filtered.vcf.gz"), path("${ont_id}.ont.indel.${size}.hc.filtered.vcf.gz.tbi"),
        path("${lp_id}.illumina.indel.${size}.hc.filtered.vcf.gz"), path("${lp_id}.illumina.indel.${size}.hc.filtered.vcf.gz.tbi"),
        path("${ont_id}.ont.indel.${size}.lc.filtered.vcf.gz"), path("${ont_id}.ont.indel.${size}.lc.filtered.vcf.gz.tbi"),
        path("${lp_id}.illumina.indel.${size}.lc.filtered.vcf.gz"), path("${lp_id}.illumina.indel.${size}.lc.filtered.vcf.gz.tbi"),
        val(size)

    script:
    def size_filter = size == 'aggregate' ? '' : "| bcftools filter -i 'abs(strlen(REF)-strlen(ALT))==${size}'"
    """
    # High confidence regions with size filter
    bcftools view -Ou -T ${high_confidence_regions} ${ont_vcf} \
        ${size_filter} \
        | bcftools view -Oz -o ${ont_id}.ont.indel.${size}.hc.filtered.vcf.gz
    bcftools index -t ${ont_id}.ont.indel.${size}.hc.filtered.vcf.gz

    bcftools view -Ou -T ${high_confidence_regions} ${illumina_vcf} \
        ${size_filter} \
        | bcftools view -Oz -o ${lp_id}.illumina.indel.${size}.hc.filtered.vcf.gz
    bcftools index -t ${lp_id}.illumina.indel.${size}.hc.filtered.vcf.gz

    # Low confidence regions with size filter
    bcftools view -Ou -T ^${high_confidence_regions} ${ont_vcf} \
        ${size_filter} \
        | bcftools view -Oz -o ${ont_id}.ont.indel.${size}.lc.filtered.vcf.gz
    bcftools index -t ${ont_id}.ont.indel.${size}.lc.filtered.vcf.gz

    bcftools view -Ou -T ^${high_confidence_regions} ${illumina_vcf} \
        ${size_filter} \
        | bcftools view -Oz -o ${lp_id}.illumina.indel.${size}.lc.filtered.vcf.gz
    bcftools index -t ${lp_id}.illumina.indel.${size}.lc.filtered.vcf.gz
    """
}

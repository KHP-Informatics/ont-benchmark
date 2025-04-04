#!/usr/bin/env nextflow

process FILTER_INDEL_VCF {
    tag "${ont_id}|${lp_id}|${size != 'aggregate' ? size + 'bp' : 'aggregate'}"

    input:
    tuple val(ont_id), val(lp_id),
        path(ont_vcf), path(ont_index),
        path(illumina_vcf), path(illumina_index),
        path(low_complexity_regions),
        path(dark_genome_regions),
        val(size)

    output:
    tuple val(ont_id), val(lp_id),
        path("${ont_id}.ont.indel.${size}.hc.filtered.vcf.gz"), path("${ont_id}.ont.indel.${size}.hc.filtered.vcf.gz.tbi"),
        path("${lp_id}.illumina.indel.${size}.hc.filtered.vcf.gz"), path("${lp_id}.illumina.indel.${size}.hc.filtered.vcf.gz.tbi"),
        path("${ont_id}.ont.indel.${size}.lc.filtered.vcf.gz"), path("${ont_id}.ont.indel.${size}.lc.filtered.vcf.gz.tbi"),
        path("${lp_id}.illumina.indel.${size}.lc.filtered.vcf.gz"), path("${lp_id}.illumina.indel.${size}.lc.filtered.vcf.gz.tbi"),
        path("${ont_id}.ont.indel.${size}.dark.filtered.vcf.gz"), path("${ont_id}.ont.indel.${size}.dark.filtered.vcf.gz.tbi"),
        path("${lp_id}.illumina.indel.${size}.dark.filtered.vcf.gz"), path("${lp_id}.illumina.indel.${size}.dark.filtered.vcf.gz.tbi"),
        val(size)

    script:
    def size_filter = size == 'aggregate' ? '' : "| bcftools filter -i 'abs(strlen(REF)-strlen(ALT))==${size}'"
    """
    # High complexity regions with size filter
    bcftools view -Ou -T ^${low_complexity_regions} ${ont_vcf} \
        ${size_filter} \
        | bcftools view -Oz -o ${ont_id}.ont.indel.${size}.hc.filtered.vcf.gz
    bcftools index -t ${ont_id}.ont.indel.${size}.hc.filtered.vcf.gz

    bcftools view -Ou -T ^${low_complexity_regions} ${illumina_vcf} \
        ${size_filter} \
        | bcftools view -Oz -o ${lp_id}.illumina.indel.${size}.hc.filtered.vcf.gz
    bcftools index -t ${lp_id}.illumina.indel.${size}.hc.filtered.vcf.gz

    # Low complexity regions with size filter
    bcftools view -Ou -T ${low_complexity_regions} ${ont_vcf} \
        ${size_filter} \
        | bcftools view -Oz -o ${ont_id}.ont.indel.${size}.lc.filtered.vcf.gz
    bcftools index -t ${ont_id}.ont.indel.${size}.lc.filtered.vcf.gz

    bcftools view -Ou -T ${low_complexity_regions} ${illumina_vcf} \
        ${size_filter} \
        | bcftools view -Oz -o ${lp_id}.illumina.indel.${size}.lc.filtered.vcf.gz
    bcftools index -t ${lp_id}.illumina.indel.${size}.lc.filtered.vcf.gz

    # Dark genome regions with size filter
    bcftools view -Ou -T ${dark_genome_regions} ${ont_vcf} \
        ${size_filter} \
        | bcftools view -Oz -o ${ont_id}.ont.indel.${size}.dark.filtered.vcf.gz
    bcftools index -t ${ont_id}.ont.indel.${size}.dark.filtered.vcf.gz

    bcftools view -Ou -T ${dark_genome_regions} ${illumina_vcf} \
        ${size_filter} \
        | bcftools view -Oz -o ${lp_id}.illumina.indel.${size}.dark.filtered.vcf.gz
    bcftools index -t ${lp_id}.illumina.indel.${size}.dark.filtered.vcf.gz
    """
}

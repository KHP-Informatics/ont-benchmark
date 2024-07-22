#!/usr/bin/env nextflow

process FILTER_INDEL_VCF {
    tag "${ont_id}|${lp_id}"

    input:
    tuple val(ont_id), val(lp_id),
        path(ont_vcf), path(ont_index),
        path(illumina_vcf), path(illumina_index),
        path(high_confidence_regions)

    output:
    tuple val(ont_id), val(lp_id),
        path("${ont_id}.ont.indel.hc.filtered.vcf.gz"), path("${ont_id}.ont.indel.hc.filtered.vcf.gz.tbi"),
        path("${lp_id}.illumina.indel.hc.filtered.vcf.gz"), path("${lp_id}.illumina.indel.hc.filtered.vcf.gz.tbi"),
        path("${ont_id}.ont.indel.lc.filtered.vcf.gz"), path("${ont_id}.ont.indel.lc.filtered.vcf.gz.tbi"),
        path("${lp_id}.illumina.indel.lc.filtered.vcf.gz"), path("${lp_id}.illumina.indel.lc.filtered.vcf.gz.tbi")

    script:
    """
    echo "Initial site counts:"
    echo "ONT VCF: \$(bcftools stats ${ont_vcf} | grep 'number of records:' | cut -f4)"
    echo "Illumina VCF: \$(bcftools stats ${illumina_vcf} | grep 'number of records:' | cut -f4)"

    bcftools view -Ou -T ${high_confidence_regions} ${ont_vcf} | \
    bcftools view -Oz -o ${ont_id}.ont.indel.hc.filtered.vcf.gz
    bcftools index -t ${ont_id}.ont.indel.hc.filtered.vcf.gz
    echo "ONT VCF after HC filtering: \$(bcftools stats ${ont_id}.ont.indel.hc.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    bcftools view -Ou -T ${high_confidence_regions} ${illumina_vcf} | \
    bcftools view -Oz -o ${lp_id}.illumina.indel.hc.filtered.vcf.gz
    bcftools index -t ${lp_id}.illumina.indel.hc.filtered.vcf.gz
    echo "Illumina VCF after HC filtering: \$(bcftools stats ${lp_id}.illumina.indel.hc.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    bcftools view -Ou -T ^${high_confidence_regions} ${ont_vcf} | \
    bcftools view -Oz -o ${ont_id}.ont.indel.lc.filtered.vcf.gz
    bcftools index -t ${ont_id}.ont.indel.lc.filtered.vcf.gz
    echo "ONT VCF after LC filtering: \$(bcftools stats ${ont_id}.ont.indel.lc.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    bcftools view -Ou -T ^${high_confidence_regions} ${illumina_vcf} | \
    bcftools view -Oz -o ${lp_id}.illumina.indel.lc.filtered.vcf.gz
    bcftools index -t ${lp_id}.illumina.indel.lc.filtered.vcf.gz
    echo "Illumina VCF after LC filtering: \$(bcftools stats ${lp_id}.illumina.indel.lc.filtered.vcf.gz | grep 'number of records:' | cut -f4)"
    """
}

#!/usr/bin/env nextflow

process FILTER_SNV_VCF {
    tag "${ont_id}|${lp_id}"

    input:
    tuple val(ont_id), val(lp_id),
        path(ont_vcf), path(illumina_vcf), path(array_vcf),
        path(ont_index), path(illumina_index), path(array_index),
        path(low_complexity_regions)

    output:
    tuple val(ont_id), val(lp_id),
        path("${ont_id}.ont.snv.filtered.vcf.gz"), path("${ont_id}.ont.snv.filtered.vcf.gz.tbi"),
        path("${lp_id}.illumina.snv.filtered.vcf.gz"), path("${lp_id}.illumina.snv.filtered.vcf.gz.tbi"),
        path("${lp_id}.array.snv.filtered.vcf.gz"), path("${lp_id}.array.snv.filtered.vcf.gz.tbi")

    script:
    """
    echo "Initial site counts:"
    echo "Array VCF: \$(bcftools stats ${array_vcf} | grep 'number of records:' | cut -f4)"
    echo "ONT VCF: \$(bcftools stats ${ont_vcf} | grep 'number of records:' | cut -f4)"
    echo "Illumina VCF: \$(bcftools stats ${illumina_vcf} | grep 'number of records:' | cut -f4)"

    bcftools view -Ou -T ${low_complexity_regions} ${array_vcf} | \
    bcftools view -Oz --write-index=tbi -o ${lp_id}.array.snv.filtered.vcf.gz
    echo "Array VCF after filtering: \$(bcftools stats ${lp_id}.array.snv.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    bcftools query -f '%CHROM\t%POS\n' ${lp_id}.array.snv.filtered.vcf.gz > array_positions.tsv
    echo "Number of positions extracted: \$(wc -l < array_positions.tsv)"

    bcftools view -Ou -T ${low_complexity_regions} ${ont_vcf} | \
    bcftools view -Ou -T array_positions.tsv | \
    bcftools view -Oz --write-index=tbi -o ${ont_id}.ont.snv.filtered.vcf.gz
    echo "ONT VCF after filtering: \$(bcftools stats ${ont_id}.ont.snv.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    bcftools view -Ou -T ${low_complexity_regions} ${illumina_vcf} | \
    bcftools view -Ou -T array_positions.tsv | \
    bcftools view -Oz --write-index=tbi -o ${lp_id}.illumina.snv.filtered.vcf.gz
    echo "Illumina VCF after filtering: \$(bcftools stats ${lp_id}.illumina.snv.filtered.vcf.gz | grep 'number of records:' | cut -f4)"
    """
}

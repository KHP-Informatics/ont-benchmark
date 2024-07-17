#!/usr/bin/env nextflow

process FILTER_VCF_SITES {
    tag "${ont_id}|${lp_id}"

    input:
    tuple val(ont_id), val(lp_id),
        path(ont_vcf), path(illumina_vcf), path(array_vcf),
        path(ont_index), path(illumina_index), path(array_index)

    output:
    tuple val(ont_id), val(lp_id),
        path("${ont_id}.ont.snv.filtered.vcf.gz"), path("${ont_id}.ont.snv.filtered.vcf.gz.tbi"),
        path("${lp_id}.illumina.snv.filtered.vcf.gz"), path("${lp_id}.illumina.snv.filtered.vcf.gz.tbi"),
        path(array_vcf), path(array_index)

    script:
    """
    echo "Initial site counts:"
    echo "Array VCF: \$(bcftools stats ${array_vcf} | grep 'number of records:' | cut -f4)"
    echo "ONT VCF: \$(bcftools stats ${ont_vcf} | grep 'number of records:' | cut -f4)"
    echo "Illumina VCF: \$(bcftools stats ${illumina_vcf} | grep 'number of records:' | cut -f4)"

    bcftools query -f '%CHROM\t%POS\n' ${array_vcf} > array_positions.txt

    echo "Number of positions extracted from array VCF: \$(wc -l < array_positions.txt)"

    bcftools view -T array_positions.txt ${ont_vcf} -Oz -o ${ont_id}.ont.snv.filtered.vcf.gz

    bcftools index -t ${ont_id}.ont.snv.filtered.vcf.gz

    echo "ONT VCF after filtering: \$(bcftools stats ${ont_id}.ont.snv.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    bcftools view -T array_positions.txt ${illumina_vcf} -Oz -o ${lp_id}.illumina.snv.filtered.vcf.gz

    bcftools index -t ${lp_id}.illumina.snv.filtered.vcf.gz

    echo "Illumina VCF after filtering: \$(bcftools stats ${lp_id}.illumina.snv.filtered.vcf.gz | grep 'number of records:' | cut -f4)"
    """
}

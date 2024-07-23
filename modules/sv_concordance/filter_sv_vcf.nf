#!/usr/bin/env nextflow

process FILTER_SV_VCF {
    tag "${ont_id}|${lp_id}"

    input:
    tuple val(ont_id), val(lp_id),
        path(ont_vcf), path(ont_index),
        path(illumina_vcf), path(illumina_index)

    output:
    tuple val(ont_id), val(lp_id),
        path("${ont_id}.ont.sv.filtered.vcf"),
        path("${lp_id}.illumina.sv.filtered.vcf")

    script:
    """
    echo "Initial site counts:"
    echo "ONT VCF: \$(bcftools stats ${ont_vcf} | grep 'number of records:' | cut -f4)"
    echo "Illumina VCF: \$(bcftools stats ${illumina_vcf} | grep 'number of records:' | cut -f4)"

    bcftools view -i 'INFO/SUPPORT>=5' -Ov -o ${ont_id}.ont.sv.filtered.vcf ${ont_vcf}
    echo "ONT VCF after filtering: \$(bcftools stats ${ont_id}.ont.sv.filtered.vcf | grep 'number of records:' | cut -f4)"

    bcftools view -i '(FORMAT/PR[0:1] + FORMAT/SR[0:1]) >= 5' -Ov -o ${lp_id}.illumina.sv.filtered.vcf ${illumina_vcf}
    echo "Illumina VCF after filtering: \$(bcftools stats ${lp_id}.illumina.sv.filtered.vcf | grep 'number of records:' | cut -f4)"
    """
}

#!/usr/bin/env nextflow

process MERGE_AND_FILTER_SV_VCF {
    tag "${ont_id}|${lp_id}"

    input:
    tuple val(ont_id), val(lp_id),
        path(ont_sv_vcf), path(ont_sv_index),
        path(illumina_vcf), path(illumina_index),
        path(ont_str_vcf), path(ont_str_index)

    output:
    tuple val(ont_id), val(lp_id),
        path("${ont_id}.ont.sv_str.filtered.vcf"),
        path("${lp_id}.illumina.sv.filtered.vcf")

    script:
    """
    echo "Initial site counts:"
    echo "ONT SV VCF: \$(bcftools stats ${ont_sv_vcf} | grep 'number of records:' | cut -f4)"
    echo "ONT STR VCF: \$(bcftools stats ${ont_str_vcf} | grep 'number of records:' | cut -f4)"
    echo "Illumina VCF: \$(bcftools stats ${illumina_vcf} | grep 'number of records:' | cut -f4)"

    bcftools view --write-index=csi -i 'INFO/SUPPORT>=5' -Ou -o ${ont_id}.ont.sv.filtered.bcf ${ont_sv_vcf}
    echo "ONT SV VCF after filtering: \$(bcftools stats ${ont_id}.ont.sv.filtered.bcf | grep 'number of records:' | cut -f4)"

    bcftools view --write-index=csi -i 'FORMAT/LC >= 10' -Ou -o ${ont_id}.ont.str.filtered.bcf ${ont_str_vcf}
    echo "ONT STR VCF after filtering: \$(bcftools stats ${ont_id}.ont.str.filtered.bcf | grep 'number of records:' | cut -f4)"

    bcftools concat -a ${ont_id}.ont.sv.filtered.bcf ${ont_id}.ont.str.filtered.bcf -Ov -o ${ont_id}.ont.sv_str.filtered.vcf
    echo "ONT merged SV+STR VCF: \$(bcftools stats ${ont_id}.ont.sv_str.merged.bcf | grep 'number of records:' | cut -f4)"

    bcftools view -i '(FORMAT/PR[0:1] + FORMAT/SR[0:1]) >= 5' -Ov -o ${lp_id}.illumina.sv.filtered.vcf ${illumina_vcf}
    echo "Illumina VCF after filtering: \$(bcftools stats ${lp_id}.illumina.sv.filtered.vcf | grep 'number of records:' | cut -f4)"
    """
}

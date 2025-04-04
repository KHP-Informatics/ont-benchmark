#!/usr/bin/env nextflow

process FILTER_SNV_VCF {
    tag "${ont_id}|${lp_id}"

    input:
    tuple val(ont_id), val(lp_id),
        path(ont_vcf), path(ont_index),
        path(illumina_vcf), path(illumina_index),
        path(array_vcf), path(array_index),
        path(low_complexity_regions),
        path(dark_genome_regions)

    output:
    tuple val(ont_id), val(lp_id),
        path("${ont_id}.ont.snv.hc.filtered.vcf.gz"), path("${ont_id}.ont.snv.hc.filtered.vcf.gz.tbi"),
        path("${lp_id}.illumina.snv.hc.filtered.vcf.gz"), path("${lp_id}.illumina.snv.hc.filtered.vcf.gz.tbi"),
        path("${lp_id}.array.snv.hc.filtered.vcf.gz"), path("${lp_id}.array.snv.hc.filtered.vcf.gz.tbi"),
        path("${ont_id}.ont.snv.lc.filtered.vcf.gz"), path("${ont_id}.ont.snv.lc.filtered.vcf.gz.tbi"),
        path("${lp_id}.illumina.snv.lc.filtered.vcf.gz"), path("${lp_id}.illumina.snv.lc.filtered.vcf.gz.tbi"),
        path("${lp_id}.array.snv.lc.filtered.vcf.gz"), path("${lp_id}.array.snv.lc.filtered.vcf.gz.tbi"),
        path("${ont_id}.ont.snv.dark.filtered.vcf.gz"), path("${ont_id}.ont.snv.dark.filtered.vcf.gz.tbi"),
        path("${lp_id}.illumina.snv.dark.filtered.vcf.gz"), path("${lp_id}.illumina.snv.dark.filtered.vcf.gz.tbi"),
        path("${lp_id}.array.snv.dark.filtered.vcf.gz"), path("${lp_id}.array.snv.dark.filtered.vcf.gz.tbi")

    script:
    """
    echo "Initial site counts:"
    echo "Array VCF: \$(bcftools stats ${array_vcf} | grep 'number of records:' | cut -f4)"
    echo "ONT VCF: \$(bcftools stats ${ont_vcf} | grep 'number of records:' | cut -f4)"
    echo "Illumina VCF: \$(bcftools stats ${illumina_vcf} | grep 'number of records:' | cut -f4)"

    # High Complexity: Not in low complexity regions
    bcftools view -Ou -T ^${low_complexity_regions} ${array_vcf} | \
    bcftools view -Oz -o ${lp_id}.array.snv.hc.filtered.vcf.gz
    bcftools index -t ${lp_id}.array.snv.hc.filtered.vcf.gz
    echo "Array VCF after HC filtering: \$(bcftools stats ${lp_id}.array.snv.hc.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    bcftools query -f '%CHROM\t%POS\n' ${lp_id}.array.snv.hc.filtered.vcf.gz > array_positions_hc.tsv
    echo "Number of HC positions extracted: \$(wc -l < array_positions_hc.tsv)"

    bcftools view -Ou -T ^${low_complexity_regions} ${ont_vcf} | \
    bcftools view -Ou -T array_positions_hc.tsv | \
    bcftools view -Oz -o ${ont_id}.ont.snv.hc.filtered.vcf.gz
    bcftools index -t ${ont_id}.ont.snv.hc.filtered.vcf.gz
    echo "ONT VCF after HC filtering: \$(bcftools stats ${ont_id}.ont.snv.hc.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    bcftools view -Ou -T ^${low_complexity_regions} ${illumina_vcf} | \
    bcftools view -Ou -T array_positions_hc.tsv | \
    bcftools view -Oz -o ${lp_id}.illumina.snv.hc.filtered.vcf.gz
    bcftools index -t ${lp_id}.illumina.snv.hc.filtered.vcf.gz
    echo "Illumina VCF after HC filtering: \$(bcftools stats ${lp_id}.illumina.snv.hc.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    # Low Complexity: In low complexity regions
    bcftools view -Ou -T ${low_complexity_regions} ${array_vcf} | \
    bcftools view -Oz -o ${lp_id}.array.snv.lc.filtered.vcf.gz
    bcftools index -t ${lp_id}.array.snv.lc.filtered.vcf.gz
    echo "Array VCF after LC filtering: \$(bcftools stats ${lp_id}.array.snv.lc.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    bcftools query -f '%CHROM\t%POS\n' ${lp_id}.array.snv.lc.filtered.vcf.gz > array_positions_lc.tsv
    echo "Number of LC positions extracted: \$(wc -l < array_positions_lc.tsv)"

    bcftools view -Ou -T ${low_complexity_regions} ${ont_vcf} | \
    bcftools view -Ou -T array_positions_lc.tsv | \
    bcftools view -Oz -o ${ont_id}.ont.snv.lc.filtered.vcf.gz
    bcftools index -t ${ont_id}.ont.snv.lc.filtered.vcf.gz
    echo "ONT VCF after LC filtering: \$(bcftools stats ${ont_id}.ont.snv.lc.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    bcftools view -Ou -T ${low_complexity_regions} ${illumina_vcf} | \
    bcftools view -Ou -T array_positions_lc.tsv | \
    bcftools view -Oz -o ${lp_id}.illumina.snv.lc.filtered.vcf.gz
    bcftools index -t ${lp_id}.illumina.snv.lc.filtered.vcf.gz
    echo "Illumina VCF after LC filtering: \$(bcftools stats ${lp_id}.illumina.snv.lc.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    # Dark Genome: In dark genome regions
    bcftools view -Ou -T ${dark_genome_regions} ${array_vcf} | \
    bcftools view -Oz -o ${lp_id}.array.snv.dark.filtered.vcf.gz
    bcftools index -t ${lp_id}.array.snv.dark.filtered.vcf.gz
    echo "Array VCF after dark genome filtering: \$(bcftools stats ${lp_id}.array.snv.dark.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    bcftools query -f '%CHROM\t%POS\n' ${lp_id}.array.snv.dark.filtered.vcf.gz > array_positions_dark.tsv
    echo "Number of dark genome positions extracted: \$(wc -l < array_positions_dark.tsv)"

    bcftools view -Ou -T ${dark_genome_regions} ${ont_vcf} | \
    bcftools view -Ou -T array_positions_dark.tsv | \
    bcftools view -Oz -o ${ont_id}.ont.snv.dark.filtered.vcf.gz
    bcftools index -t ${ont_id}.ont.snv.dark.filtered.vcf.gz
    echo "ONT VCF after dark genome filtering: \$(bcftools stats ${ont_id}.ont.snv.dark.filtered.vcf.gz | grep 'number of records:' | cut -f4)"

    bcftools view -Ou -T ${dark_genome_regions} ${illumina_vcf} | \
    bcftools view -Ou -T array_positions_dark.tsv | \
    bcftools view -Oz -o ${lp_id}.illumina.snv.dark.filtered.vcf.gz
    bcftools index -t ${lp_id}.illumina.snv.dark.filtered.vcf.gz
    echo "Illumina VCF after dark genome filtering: \$(bcftools stats ${lp_id}.illumina.snv.dark.filtered.vcf.gz | grep 'number of records:' | cut -f4)"
    """
}

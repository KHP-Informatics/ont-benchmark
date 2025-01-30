#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FILTER_INDEL_VCF } from '../modules/indel_benchmark/filter_indel_vcf.nf'
include { RTG_VCFEVAL } from '../modules/shared/rtg_vcfeval.nf'
include { RTG_ROCPLOT } from '../modules/shared/rtg_rocplot.nf'

/*
========================================================================================
    INITIALISE PARAMETER CHANNELS
========================================================================================
*/

Channel
    .fromPath(params.high_confidence_regions)
    .set { high_confidence_regions_ch }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW: INDEL_BENCHMARK
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INDEL_BENCHMARK {
    take:
    indel_samples_ch
    reference_sdf_ch

    main:
    def variant_type = 'indel'

    filter_indel_vcf_ch = indel_samples_ch.combine(high_confidence_regions_ch)

    FILTER_INDEL_VCF(filter_indel_vcf_ch)

    hc_comparison = FILTER_INDEL_VCF.out.map { ont_id, lp_id, ont_vcf_hc, ont_index_hc, illumina_vcf_hc, illumina_index_hc, ont_vcf_lc, ont_index_lc, illumina_vcf_lc, illumina_index_lc ->
        tuple(ont_id, lp_id, ont_vcf_hc, ont_index_hc, illumina_vcf_hc, illumina_index_hc, 'indel', 'hc')
    }

    lc_comparison = FILTER_INDEL_VCF.out.map { ont_id, lp_id, ont_vcf_hc, ont_index_hc, illumina_vcf_hc, illumina_index_hc, ont_vcf_lc, ont_index_lc, illumina_vcf_lc, illumina_index_lc ->
        tuple(ont_id, lp_id, ont_vcf_lc, ont_index_lc, illumina_vcf_lc, illumina_index_lc, 'indel', 'lc')
    }

    comparison_ch = hc_comparison.mix(lc_comparison)

    RTG_VCFEVAL(
        comparison_ch.combine(reference_sdf_ch),
        variant_type
    )

    weighted_roc_hc = RTG_VCFEVAL.out.weighted_roc.filter { it.toString().contains("_hc_") }.collect()
    weighted_roc_lc = RTG_VCFEVAL.out.weighted_roc.filter { it.toString().contains("_lc_") }.collect()

    RTG_ROCPLOT(
        weighted_roc_hc,
        weighted_roc_lc,
        variant_type
    )
}

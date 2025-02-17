#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FILTER_SNV_VCF } from '../modules/snv_benchmark/filter_snv_vcf.nf'
include { RTG_VCFEVAL } from '../modules/shared/rtg_vcfeval.nf'

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
    WORKFLOW: SNV_BENCHMARK
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SNV_BENCHMARK {
    take:
    snv_samples_ch
    reference_sdf_ch

    main:
    def variant_type = 'snv'

    filter_snv_vcf__ch = snv_samples_ch.combine(high_confidence_regions_ch)

    FILTER_SNV_VCF(filter_snv_vcf__ch)

    comparison_ch = FILTER_SNV_VCF.out.flatMap { result ->
        def (ont_id, lp_id,
             ont_vcf_hc, ont_index_hc, illumina_vcf_hc, illumina_index_hc, array_vcf_hc, array_index_hc,
             ont_vcf_lc, ont_index_lc, illumina_vcf_lc, illumina_index_lc, array_vcf_lc, array_index_lc) = result

        return [
            tuple(ont_id, lp_id, ont_vcf_hc, ont_index_hc, array_vcf_hc, array_index_hc, 'snv', 'hc', 'aggregate'),
            tuple(lp_id, ont_id, illumina_vcf_hc, illumina_index_hc, array_vcf_hc, array_index_hc, 'snv', 'hc', 'aggregate'),
            tuple(ont_id, lp_id, ont_vcf_lc, ont_index_lc, array_vcf_lc, array_index_lc, 'snv', 'lc', 'aggregate'),
            tuple(lp_id, ont_id, illumina_vcf_lc, illumina_index_lc, array_vcf_lc, array_index_lc, 'snv', 'lc', 'aggregate')
        ]
    }
    .combine(reference_sdf_ch)

    RTG_VCFEVAL(comparison_ch, variant_type)
}

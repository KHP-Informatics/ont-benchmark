#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FILTER_INDEL_VCF } from '../modules/indel_benchmark/filter_indel_vcf.nf'
include { RTG_VCFEVAL } from '../modules/shared/rtg_vcfeval.nf'

/*
========================================================================================
    INITIALISE PARAMETER CHANNELS
========================================================================================
*/

Channel
    .fromPath(params.high_confidence_regions)
    .set { high_confidence_regions_ch }

Channel
    .of('aggregate', *1..50)
    .set { sizes_ch }

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

    filter_indel_vcf_ch = indel_samples_ch
        .combine(high_confidence_regions_ch)
        .combine(sizes_ch)
        .map { it.flatten() }

    FILTER_INDEL_VCF(filter_indel_vcf_ch)

    comparison_ch = FILTER_INDEL_VCF.out.flatMap { result ->
        def (ont_id, lp_id,
            ont_vcf_hc, ont_index_hc, illumina_vcf_hc, illumina_index_hc,
            ont_vcf_lc, ont_index_lc, illumina_vcf_lc, illumina_index_lc,
            size) = result
        return [
            tuple(ont_id, lp_id, ont_vcf_hc, ont_index_hc, illumina_vcf_hc, illumina_index_hc, 'indel', 'hc', size),
            tuple(ont_id, lp_id, ont_vcf_lc, ont_index_lc, illumina_vcf_lc, illumina_index_lc, 'indel', 'lc', size)
        ]
    }
    .combine(reference_sdf_ch)

    RTG_VCFEVAL(comparison_ch, variant_type)
}

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
    .fromPath(params.low_complexity_regions)
    .set { low_complexity_regions_ch }

Channel
    .fromPath(params.dark_genome_regions)
    .set { dark_genome_regions_ch }

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

    filter_snv_vcf__ch = snv_samples_ch
        .combine(low_complexity_regions_ch)
        .combine(dark_genome_regions_ch)

    FILTER_SNV_VCF(filter_snv_vcf__ch)

    comparison_ch = FILTER_SNV_VCF.out.flatMap { result ->
        def (ont_id, lp_id,
            ont_vcf_hc, ont_index_hc, illumina_vcf_hc, illumina_index_hc, array_vcf_hc, array_index_hc,
            ont_vcf_lc, ont_index_lc, illumina_vcf_lc, illumina_index_lc, array_vcf_lc, array_index_lc,
            ont_vcf_dark, ont_index_dark, illumina_vcf_dark, illumina_index_dark, array_vcf_dark, array_index_dark) = result

        return [
            tuple(ont_id, lp_id, ont_vcf_hc, ont_index_hc, array_vcf_hc, array_index_hc, 'snv', 'hc', 'aggregate'),
            tuple(lp_id, ont_id, illumina_vcf_hc, illumina_index_hc, array_vcf_hc, array_index_hc, 'snv', 'hc', 'aggregate'),
            tuple(ont_id, lp_id, ont_vcf_lc, ont_index_lc, array_vcf_lc, array_index_lc, 'snv', 'lc', 'aggregate'),
            tuple(lp_id, ont_id, illumina_vcf_lc, illumina_index_lc, array_vcf_lc, array_index_lc, 'snv', 'lc', 'aggregate'),
            tuple(ont_id, lp_id, ont_vcf_dark, ont_index_dark, array_vcf_dark, array_index_dark, 'snv', 'dark', 'aggregate'),
            tuple(lp_id, ont_id, illumina_vcf_dark, illumina_index_dark, array_vcf_dark, array_index_dark, 'snv', 'dark', 'aggregate')
        ]
    }
    .combine(reference_sdf_ch)

    RTG_VCFEVAL(comparison_ch, variant_type)
}

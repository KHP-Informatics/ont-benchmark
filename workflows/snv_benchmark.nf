#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FILTER_SNV_VCF } from '../modules/snv_benchmark/filter_snv_vcf.nf'
include { RTG_VCFEVAL } from '../modules/snv_benchmark/rtg_vcfeval.nf'

/*
========================================================================================
    INITIALISE PARAMETER CHANNELS
========================================================================================
*/

Channel
    .fromPath(params.low_complexity_regions)
    .set { low_complexity_regions_ch }

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
    filter_vcf_input_ch = snv_samples_ch.combine(low_complexity_regions_ch)

    FILTER_SNV_VCF(filter_vcf_input_ch)

    ont_comparison = FILTER_SNV_VCF.out.map { ont_id, lp_id, ont_vcf, ont_index, illumina_vcf, illumina_index, array_vcf, array_index ->
        tuple(ont_id, lp_id, ont_vcf, ont_index, array_vcf, array_index, 'ont')
    }

    illumina_comparison = FILTER_SNV_VCF.out.map { ont_id, lp_id, ont_vcf, ont_index, illumina_vcf, illumina_index, array_vcf, array_index ->
        tuple(lp_id, ont_id, illumina_vcf, illumina_index, array_vcf, array_index, 'illumina')
    }

    comparison_ch = ont_comparison.mix(illumina_comparison)

    RTG_VCFEVAL(
        comparison_ch.combine(reference_sdf_ch)
    )
}

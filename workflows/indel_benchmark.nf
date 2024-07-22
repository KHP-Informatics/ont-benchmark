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
    filter_indel_vcf_ch = indel_samples_ch.combine(high_confidence_regions_ch)

    FILTER_INDEL_VCF(filter_indel_vcf_ch)

}

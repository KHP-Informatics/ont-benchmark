#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FILTER_VCF_SITES } from '../modules/snv_benchmark/filter_vcf_sites.nf'
include { RTG_VCFEVAL } from '../modules/snv_benchmark/rtg_vcfeval.nf'

/*
========================================================================================
    INITIALISE PARAMETER CHANNELS
========================================================================================
*/

// N/A

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
    FILTER_VCF_SITES(snv_samples_ch)

/*
    RTG_VCFEVAL(
        FILTER_VCF_SITES.out,
        reference_sdf_ch
    )

    vcfeval_results = RTG_VCFEVAL.out.eval_dir

    emit:
    vcfeval_results
*/
}

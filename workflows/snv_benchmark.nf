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

    ont_comparison = FILTER_VCF_SITES.out.map { ont_id, lp_id, ont_vcf, ont_index, illumina_vcf, illumina_index, array_vcf, array_index ->
        tuple(ont_id, lp_id, ont_vcf, ont_index, array_vcf, array_index, 'ont')
    }

    illumina_comparison = FILTER_VCF_SITES.out.map { ont_id, lp_id, ont_vcf, ont_index, illumina_vcf, illumina_index, array_vcf, array_index ->
        tuple(lp_id, ont_id, illumina_vcf, illumina_index, array_vcf, array_index, 'illumina')
    }

    comparison_ch = ont_comparison.mix(illumina_comparison)

    RTG_VCFEVAL(
        comparison_ch.combine(reference_sdf_ch)
    )
}

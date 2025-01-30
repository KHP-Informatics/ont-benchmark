#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FILTER_SNV_VCF } from '../modules/snv_benchmark/filter_snv_vcf.nf'
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

    ont_comparison_hc = FILTER_SNV_VCF.out.map { ont_id, lp_id, ont_vcf_hc, ont_index_hc, illumina_vcf_hc, illumina_index_hc, array_vcf_hc, array_index_hc, ont_vcf_lc, ont_index_lc, illumina_vcf_lc, illumina_index_lc, array_vcf_lc, array_index_lc ->
        tuple(ont_id, lp_id, ont_vcf_hc, ont_index_hc, array_vcf_hc, array_index_hc, 'ont', 'hc')
    }

    illumina_comparison_hc = FILTER_SNV_VCF.out.map { ont_id, lp_id, ont_vcf_hc, ont_index_hc, illumina_vcf_hc, illumina_index_hc, array_vcf_hc, array_index_hc, ont_vcf_lc, ont_index_lc, illumina_vcf_lc, illumina_index_lc, array_vcf_lc, array_index_lc ->
        tuple(lp_id, ont_id, illumina_vcf_hc, illumina_index_hc, array_vcf_hc, array_index_hc, 'illumina', 'hc')
    }

    ont_comparison_lc = FILTER_SNV_VCF.out.map { ont_id, lp_id, ont_vcf_hc, ont_index_hc, illumina_vcf_hc, illumina_index_hc, array_vcf_hc, array_index_hc, ont_vcf_lc, ont_index_lc, illumina_vcf_lc, illumina_index_lc, array_vcf_lc, array_index_lc ->
        tuple(ont_id, lp_id, ont_vcf_lc, ont_index_lc, array_vcf_lc, array_index_lc, 'ont', 'lc')
    }

    illumina_comparison_lc = FILTER_SNV_VCF.out.map { ont_id, lp_id, ont_vcf_hc, ont_index_hc, illumina_vcf_hc, illumina_index_hc, array_vcf_hc, array_index_hc, ont_vcf_lc, ont_index_lc, illumina_vcf_lc, illumina_index_lc, array_vcf_lc, array_index_lc ->
        tuple(lp_id, ont_id, illumina_vcf_lc, illumina_index_lc, array_vcf_lc, array_index_lc, 'illumina', 'lc')
    }

    comparison_ch = ont_comparison_hc.mix(illumina_comparison_hc, ont_comparison_lc, illumina_comparison_lc)

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

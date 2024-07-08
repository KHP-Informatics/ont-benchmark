#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

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
    snv_samples
        .flatMap { ont_id, lp_id, variant, files ->
            files.findAll { it[0] in ['ont', 'illumina'] }
            .collect { type, vcf, index ->
                [ont_id, lp_id, type, vcf, index]
            }
        }
        .set { split_snv_samples }

    snv_samples
        .map { ont_id, lp_id, variant, files ->
            def microarray = files.find { it[0] == 'microarray' }
            [ont_id, lp_id, microarray[1], microarray[2]]
        }
        .set { microarray_files }

    split_snv_samples
        .combine(microarray_files, by: [0, 1])
        .set { snv_comparison_inputs }

    RTG_VCFEVAL(
        snv_samples_ch,
        reference_sdf_ch.first()
    )

    emit:
    vcfeval_results = RTG_VCFEVAL.out.eval_dir
}

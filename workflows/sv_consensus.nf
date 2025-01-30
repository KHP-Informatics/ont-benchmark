#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { MERGE_AND_FILTER_SV_VCF } from '../modules/sv_consensus/merge_and_filter_sv_vcf.nf'
include { SURVIVOR } from '../modules/sv_consensus/survivor.nf'

/*
========================================================================================
    INITIALISE PARAMETER CHANNELS
========================================================================================
*/

Channel
    .fromPath( params.reference_fasta, checkIfExists: true )
    .set { reference_fasta_ch }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW: SV_CONSENSUS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SV_CONSENSUS {
    take:
    sv_samples_ch
    str_samples_ch

    main:
    combined_samples_ch = sv_samples_ch
        .join(str_samples_ch, by: [0, 1])
        .map { ont_id,
                lp_id,
                ont_sv_vcf,
                ont_sv_index,
                illumina_vcf,
                illumina_index,
                ont_str_vcf,
                ont_str_index,
                unused1,
                unused2 ->
            tuple(
                ont_id,
                lp_id,
                ont_sv_vcf,
                ont_sv_index,
                illumina_vcf,
                illumina_index,
                ont_str_vcf,
                ont_str_index
            )
        }

    MERGE_AND_FILTER_SV_VCF(combined_samples_ch)

    SURVIVOR(MERGE_AND_FILTER_SV_VCF.out)
}

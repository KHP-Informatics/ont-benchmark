#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FILTER_SV_VCF } from '../modules/sv_concordance/filter_sv_vcf.nf'
include { SURVIVOR } from '../modules/sv_concordance/survivor.nf'

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
    WORKFLOW: SV_CONCORDANCE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SV_CONCORDANCE {
    take:
    sv_samples_ch

    main:
    FILTER_SV_VCF(sv_samples_ch)
    SURVIVOR(FILTER_SV_VCF.out)
}

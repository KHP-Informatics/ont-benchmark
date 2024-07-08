#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { SPLIT_SNV_INDELS } from '../modules/setup/split_snv_indels.nf'
include { COLLECT_UNIQUE_MICROARRAY_VARIANT_IDS } from '../modules/setup/collect_unique_microarray_variant_ids.nf'
include { CONVERT_TO_RSIDS } from '../modules/setup/convert_to_rsids.nf'
include { QUERY_RSID_POSITIONS } from '../modules/setup/query_rsid_positions.nf'
include { UPDATE_MICROARRAY_VCF } from '../modules/setup/update_microarray_vcf.nf'
include { INDEX_VCF as INDEX_INPUT_VCF } from '../modules/setup/index_vcf.nf'
include { GENERATE_SDF_REFERENCE } from '../modules/setup/generate_sdf_reference.nf'

/*
========================================================================================
    INITIALISE PARAMETER CHANNELS
========================================================================================
*/

// Sample IDs
Channel
    .fromPath('./sample_ids.csv', checkIfExists: true)
    .splitCsv(header: true)
    .map { row -> tuple(row.ont_id, row.lp_id) }
    .set { sample_ids_ch }

// Input VCF files
sample_ids_ch.flatMap { ont_id, lp_id ->
    def files = []
    def missing_files = []

    // Microarray
    def microarray_file = file("${params.microarray_dir}/FinalReport_InfiniumOmni2-5-8v1-4_${lp_id}.vcf.gz")
    if (microarray_file.exists()) {
        files << tuple([id: lp_id, type: 'microarray', variant: 'snv'], microarray_file)
    } else {
        missing_files << "Microarray file for ${lp_id}"
    }

    // Illumina
    def illumina_file = file("${params.illumina_dir}/${lp_id}.vcf.gz")
    if (illumina_file.exists()) {
        files << tuple([id: lp_id, type: 'illumina', variant: 'snv_indel'], illumina_file)
    } else {
        missing_files << "Illumina SNV/INDEL file for ${lp_id}"
    }

    def illumina_sv_file = file("${params.illumina_dir}/${lp_id}.SV.vcf.gz")
    if (illumina_sv_file.exists()) {
        files << tuple([id: lp_id, type: 'illumina', variant: 'sv'], illumina_sv_file)
    } else {
        missing_files << "Illumina SV file for ${lp_id}"
    }

    // ONT
    ['snp', 'sv', 'str', 'cnv'].each { variant ->
        def ont_file = file("${params.ont_dir}/${ont_id}_${params.basecall}/${ont_id}_${params.basecall}.wf_${variant}.vcf.gz")
        if (ont_file.exists()) {
            files << tuple([id: ont_id, type: 'ont', variant: variant], ont_file)
        } else {
            missing_files << "ONT ${variant.toUpperCase()} file for ${ont_id}"
        }
    }

    if (!missing_files.isEmpty()) {
        error "The following required files are missing:\n${missing_files.join('\n')}\nPlease ensure all required files exist before running the pipeline."
    }

    return files
    }.branch {
        microarray: it[0].type == 'microarray'
        illumina_snv_indel: it[0].type == 'illumina' && it[0].variant == 'snv_indel'
        illumina_sv: it[0].type == 'illumina' && it[0].variant == 'sv'
        ont_snv_indel: it[0].type == 'ont' && it[0].variant == 'snp'
        ont_sv: it[0].type == 'ont' && it[0].variant == 'sv'
        ont_str: it[0].type == 'ont' && it[0].variant == 'str'
        ont_cnv: it[0].type == 'ont' && it[0].variant == 'cnv'
    }.set { all_vcf_files }

microarray_ch = all_vcf_files.microarray
illumina_snv_indel_ch = all_vcf_files.illumina_snv_indel
illumina_sv_ch = all_vcf_files.illumina_sv
ont_snv_indel_ch = all_vcf_files.ont_snv_indel
ont_sv_ch = all_vcf_files.ont_sv
ont_str_ch = all_vcf_files.ont_str
ont_cnv_ch = all_vcf_files.ont_cnv

// FASTA genome reference
Channel
    .fromPath( params.reference_fasta, checkIfExists: true )
    .set { reference_fasta_ch }

Channel
    .fromPath( params.array_positions_file, checkIfExists: true )
    .set { array_positions_ch }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW: SETUP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SETUP {
    main:
    snv_indel_files = illumina_snv_indel_ch.mix(ont_snv_indel_ch)
    SPLIT_SNV_INDELS(snv_indel_files)

    snv_ch = SPLIT_SNV_INDELS.out.snv.map { meta, file ->
        meta.variant = 'snv'
        [meta, file]
    }

    indel_ch = SPLIT_SNV_INDELS.out.indel.map { meta, file ->
        meta.variant = 'indel'
        [meta, file]
    }

    microarray_vcfs_ch = microarray_ch.map { meta, vcf -> vcf }.collect()
    COLLECT_UNIQUE_MICROARRAY_VARIANT_IDS(microarray_vcfs_ch)

    CONVERT_TO_RSIDS(
        COLLECT_UNIQUE_MICROARRAY_VARIANT_IDS.out.unique_variant_ids,
        array_positions_ch
    )

    QUERY_RSID_POSITIONS(
        CONVERT_TO_RSIDS.out.unique_rsids,
        array_positions_ch
    )

    UPDATE_MICROARRAY_VCF(
        microarray_ch
            .combine(QUERY_RSID_POSITIONS.out.rsid_positions)
            .combine(reference_fasta_ch)
    )

    all_vcf_files_ch = illumina_sv_ch
        .mix(ont_sv_ch)
        .mix(ont_str_ch)
        .mix(ont_cnv_ch)
        .mix(SPLIT_SNV_INDELS.out.snv)
        .mix(SPLIT_SNV_INDELS.out.indel)
        .mix(UPDATE_MICROARRAY_VCF.out.pos_vcf)

    INDEX_INPUT_VCF(
        all_vcf_files_ch
    )

    def group_by_sample_and_variant = INDEX_INPUT_VCF.out
        .groupTuple()
        .flatMap { entry ->
            def (id, files) = entry
            files.groupBy { it.meta.id }.collect { ont_id, variants ->
                variants.groupBy { it.meta.variant }.collect { variant, variant_files ->
                    [ont_id, variant_files.first().meta.lp_id, variant, variant_files.collect { [it.meta.type, it.path, it.index] }]
                }
            }
        }
        .groupTuple(by: [0, 1, 2])
        .map { ont_id, lp_id, variant, file_groups ->
            [ont_id, lp_id, variant, file_groups.flatten(1)]
        }

    // Create separate channels for each variant type
    group_by_sample_and_variant
        .branch {
            snv: it[2] == 'snv'
            indel: it[2] == 'indel'
            sv: it[2] == 'sv'
            str: it[2] == 'str'
            cnv: it[2] == 'cnv'
        }
        .set { grouped_variants_ch }

    GENERATE_SDF_REFERENCE(
        reference_fasta_ch
    )

    emit:
    snv_samples_ch = grouped_variants_ch.snv
    indel_samples_ch = grouped_variants_ch.indel
    sv_samples_ch = grouped_variants_ch.sv
    str_samples_ch = grouped_variants_ch.str
    cnv_samples_ch = grouped_variants_ch.cnv
    reference_sdf_ch = GENERATE_SDF_REFERENCE.out.reference_sdf
}

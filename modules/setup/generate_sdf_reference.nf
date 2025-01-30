#!/usr/bin/env nextflow

process GENERATE_SDF_REFERENCE {
    input:
    path reference_fasta

    output:
    path "${reference_base}.sdf", type: 'dir', emit: reference_sdf

    script:
    reference_base = reference_fasta.baseName
    """
    rtg format \
        --format=fasta \
        --output ${reference_base}.sdf \
        ${reference_fasta}
    """
}
#!/usr/bin/env nextflow

process RTG_ROCPLOT {
    publishDir "${params.outdir}/snv/rtg_rocplot/", mode: 'copy'

    input:
    path(weighted_roc_files)

    output:
    path "combined_snv_roc_plot.svg", emit: roc_plot

    script:
    """
    rtg rocplot --svg=combined_snv_roc_plot.svg \
        --plain \
        --palette=classic \
        ${weighted_roc_files.join(' ')}
    """
}

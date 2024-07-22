#!/usr/bin/env nextflow

process RTG_ROCPLOT {
    publishDir "${params.outdir}/snv/rtg_rocplot/", mode: 'copy'

    input:
    path(weighted_roc_files_hc)
    path(weighted_roc_files_lc)

    output:
    path "combined_roc_plot_hc.svg", emit: roc_plot_hc
    path "combined_roc_plot_lc.svg", emit: roc_plot_lc

    script:
    """
    rtg rocplot --svg=combined_roc_plot_hc.svg \
        --plain \
        --palette=classic \
        ${weighted_roc_files_hc.join(' ')}

    rtg rocplot --svg=combined_roc_plot_lc.svg \
        --plain \
        --palette=classic \
        ${weighted_roc_files_lc.join(' ')}
    """
}

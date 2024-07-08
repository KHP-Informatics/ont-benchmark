#!/usr/bin/env nextflow

process RTG_VCFEVAL {
    publishDir "${params.outdir}/snv/${ont_id}.${lp_id}/", mode: 'copy'

    input:
    tuple val(ont_id), val(lp_id), val(query_type), path(query_vcf), path(query_index), path(truth_vcf), path(truth_index)
    path reference_sdf

    output:
    path "${query_type}_eval/*", emit: eval_dir

    script:
    """
    rtg vcfeval --baseline=${truth_vcf} \
                --calls=${query_vcf} \
                --template=${reference_sdf} \
                --output=${query_type}_eval \
                --output-mode=roc-only
    """
}

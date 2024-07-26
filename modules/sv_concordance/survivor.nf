process SURVIVOR {
    tag "${ont_id}|${lp_id}"

    publishDir "${params.outdir}/sv/survivor/${ont_id}/", mode: 'copy'

    input:
    tuple val(ont_id), val(lp_id),
        path(ont_vcf), path(illumina_vcf)

    output:
    tuple val(ont_id), val(lp_id),
        path("${ont_id}_${lp_id}_merged.vcf"),
        path(ont_vcf),
        path(illumina_vcf),
        emit: sv_vcfs

    script:
    """
    ls *.vcf > sample_files
    SURVIVOR merge sample_files 1000 1 1 1 0 30 ${ont_id}_${lp_id}_merged.vcf
    """
}

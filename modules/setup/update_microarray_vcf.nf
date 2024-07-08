#!/usr/bin/env nextflow

process UPDATE_MICROARRAY_VCF {
    tag "${meta.id}"

    input:
    tuple val(meta), path(vcf)
    path rsid_positions

    output:
    tuple val(meta), path("${meta.id}.pos.vcf.gz"), emit: pos_vcf

    script:
    """
    #!/usr/bin/env python
    import pysam

    rsid_to_position = {}
    with open("${rsid_positions}", "r") as f:
        for line in f:
            rsid, position = line.strip().split("\t")
            chrom, pos = position.split(":")
            rsid_to_position[rsid] = (chrom, int(pos))

    with pysam.VariantFile("${vcf}") as vcf_in:
        vcf_in.header.info.add(
            "reference",
            "1",
            "String",
            "Reference genome used for the VCF file."
        )
        for chrom in set(chrom for chrom, _ in rsid_to_position.values()):
            vcf_in.header.contigs.add(chrom)

        with pysam.VariantFile(
            "${meta.id}.pos.vcf.gz", "wz", header=vcf_in.header
        ) as vcf_out:
            for record in vcf_in:
                if record.id in rsid_to_position:
                    chrom, pos = rsid_to_position[record.id]
                    record.chrom = chrom
                    record.pos = pos
                vcf_out.write(record)
    """
}

#!/usr/bin/env nextflow

process UPDATE_MICROARRAY_VCF {
    tag "${meta.id}"

    input:
    tuple val(meta), path(vcf), path(rsid_positions), path(reference_fasta)

    output:
    tuple val(meta), path("${meta.id}.pos.vcf.gz"), emit: pos_vcf

    script:
    """
    #!/usr/bin/env python
    import csv
    import pysam

    reference_contigs = {}
    with pysam.FastaFile("${reference_fasta}") as fasta:
        for contig in fasta.references:
            reference_contigs[contig] = fasta.get_reference_length(contig)

    variant_to_rsid_pos = {}
    with open("${rsid_positions}", "r", newline="") as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            variant, rsid, chrom, pos = row
            variant_to_rsid_pos[variant] = (rsid, chrom, int(pos))

    with pysam.VariantFile("${vcf}") as vcf_in:
        new_header = vcf_in.header.copy()
        for contig, length in reference_contigs.items():
            if contig in new_header.contigs:
                new_header.contigs[contig].length = length
            else:
                new_header.contigs.add(contig, length=length)

        new_header.add_line(f"##reference=file://${reference_fasta}")

        with pysam.VariantFile(
            "${meta.id}_intermediate.bcf", "wb", header=new_header
        ) as bcf_out:
            pass  # We only write the header at this stage

    with pysam.VariantFile("${vcf}") as vcf_in, pysam.VariantFile(
        "${meta.id}_intermediate.bcf", "rb"
    ) as bcf_in, pysam.VariantFile(
        "${meta.id}.pos.vcf.gz", "wz", header=bcf_in.header
    ) as vcf_out:
        for record in vcf_in:
            variant_id = record.id
            if variant_id in variant_to_rsid_pos:
                rsid, chrom, pos = variant_to_rsid_pos[variant_id]
                if chrom in reference_contigs and isinstance(pos, int):
                    record.chrom = chrom
                    record.pos = pos
                    record.id = rsid
                    vcf_out.write(record)
                else:
                    print(
                        f"Warning: Contig {chrom} or position {pos} not valid. Skipping record {variant_id}"
                    )
            else:
                print(
                    f"Warning: Variant {variant_id} not found in rsid_positions file or missing essential data. Skipping."
                )
    """
}

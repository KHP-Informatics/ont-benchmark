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
    import sys
    import csv
    import pysam


    def normalize_chromosome(chrom):
        return chrom if chrom.startswith("chr") else f"chr{chrom}"


    reference_contigs = {}
    with pysam.FastaFile("${reference_fasta}") as fasta:
        for contig in fasta.references:
            reference_contigs[normalize_chromosome(contig)] = fasta.get_reference_length(
                contig
            )

    variant_to_rsid_pos = {}
    with open("${rsid_positions}", "r", newline="") as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            variant, rsid, chrom, pos = row
            variant_to_rsid_pos[variant] = (rsid, normalize_chromosome(chrom), int(pos))

    with pysam.VariantFile("${vcf}") as vcf_in:
        new_header = vcf_in.header.copy()
        for contig, length in reference_contigs.items():
            if contig in new_header.contigs:
                new_header.contigs[contig].length = length
            else:
                new_header.contigs.add(contig, length=length)
        new_header.add_line(f"##reference=file://${reference_fasta}")

        with pysam.VariantFile("${meta.id}.pos.vcf.gz", "wz", header=new_header) as vcf_out:
            for record in vcf_in:
                variant_id = record.id
                if variant_id in variant_to_rsid_pos:
                    rsid, chrom, pos = variant_to_rsid_pos[variant_id]
                    if chrom in reference_contigs:
                        new_record = vcf_out.new_record()
                        new_record.chrom = chrom
                        new_record.pos = pos
                        new_record.id = rsid
                        new_record.alleles = record.alleles
                        new_record.qual = record.qual
                        new_record.filter.add(record.filter.keys()[0])
                        for key, value in record.info.items():
                            new_record.info[key] = value
                        for sample in record.samples:
                            new_record.samples[sample].update(record.samples[sample])
                        vcf_out.write(new_record)
                    else:
                        print(
                            f"Warning: Contig {chrom} not valid. Skipping record {variant_id}",
                            file=sys.stderr,
                        )
                else:
                    print(
                        f"Warning: Variant {variant_id} not found in rsid_positions file. Skipping.",
                        file=sys.stderr,
                    )
    """
}

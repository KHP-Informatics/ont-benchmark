#!/usr/bin/env nextflow

process COLLECT_UNIQUE_MICROARRAY_VARIANT_IDS {
    input:
    path vcfs

    output:
    path "unique_variant_ids.txt", emit: unique_variant_ids

    script:
    """
    #!/usr/bin/env python
    import pysam

    unique_ids = set()

    for vcf_file in "${vcfs}".split():
        with pysam.VariantFile(str(vcf_file)) as vcf:
            for record in vcf:
                variant_id = record.id
                if variant_id is not None and variant_id != ".":
                    unique_ids.add(variant_id)

    with open("unique_variant_ids.txt", "w") as outfile:
        for variant_id in sorted(unique_ids):
            outfile.write(f"{variant_id}\\n")

    print(f"Collected {len(unique_ids)} unique variant IDs.")
    """
}

#!/usr/bin/env nextflow

process CONVERT_TO_RSIDS {
    input:
    path unique_variant_ids
    path array_positions_file

    output:
    path "unique_rsids.txt", emit: unique_rsids

    script:
    """
    #!/usr/bin/env python
    import csv

    id_to_rsid = {}
    with open("${array_positions_file}", "r") as positions_file:
        reader = csv.DictReader(positions_file, delimiter="\t")
        for row in reader:
            if row["RsID"] != ".":
                id_to_rsid[row["Name"]] = row["RsID"].split(",")

    unique_rsids = set()
    with open("${unique_variant_ids}", "r") as in_file:
        for line in in_file:
            variant_id = line.strip()
            if variant_id.startswith("rs"):
                unique_rsids.add(variant_id)
            elif variant_id in id_to_rsid:
                unique_rsids.update(id_to_rsid[variant_id])
            # If the ID is not in the mapping and doesn't start with 'rs', it's skipped

    with open("unique_rsids.txt", "w") as out_file:
        for rsid in sorted(unique_rsids):
            out_file.write(f"{rsid}\\n")

    print(f"Converted {len(unique_rsids)} unique rsIDs.")
    """
}

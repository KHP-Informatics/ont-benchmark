#!/usr/bin/env nextflow

process QUERY_RSID_POSITIONS {
    input:
    path rsids_file
    val dbsnp_build

    output:
    path "rsid_positions.txt", emit: rsid_positions

    secret 'NCBI_EMAIL'
    secret 'NCBI_API_KEY'

    script:
    """
    #!/usr/bin/env python
    from Bio import Entrez
    import time
    import os

    Entrez.email = os.environ["NCBI_EMAIL"]
    Entrez.api_key = os.environ["NCBI_API_KEY"]


    def fetch_positions(rsids, dbsnp_build):
        query = " OR ".join(f"{rsid}[RSNumber]" for rsid in rsids)
        query += f' AND "{dbsnp_build}"[Build]'

        handle = Entrez.esearch(db="snp", term=query, retmax=len(rsids))
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return {}

        handle = Entrez.efetch(
            db="snp", id=",".join(record["IdList"]), rettype="json", retmode="text"
        )
        snp_records = Entrez.read(handle)
        handle.close()

        positions = {}
        for snp in snp_records["results"]:
            rsid = f"rs{snp['refsnp_id']}"
            placements = snp.get("primary_snapshot_data", {}).get(
                "placements_with_allele", []
            )

            if placements:
                for assembly in placements:
                    if (
                        assembly.get("seq_id_traits_by_assembly", [{}])[0].get(
                            "assembly_name"
                        )
                        == "GRCh38"
                    ):
                        position = assembly["placement_annot"]["seq_id_traits"]["position"]
                        positions[rsid] = position
                        break

                if rsid not in positions:
                    position = placements[0]["placement_annot"]["seq_id_traits"]["position"]
                    positions[rsid] = position

        return positions


    with open("${rsids_file}", "r") as f:
        rsids = [line.strip() for line in f if line.strip()]

    batch_size = 100
    all_positions = {}

    for i in range(0, len(rsids), batch_size):
        batch = rsids[i : i + batch_size]
        positions = fetch_positions(batch, "${dbsnp_build}")
        all_positions.update(positions)

    with open("rsid_positions.txt", "w") as f:
        for rsid, position in all_positions.items():
            f.write(f"{rsid}\\t{position}\\n")

    print(f"Queried positions for {len(all_positions)} rsIDs.")
    """
}

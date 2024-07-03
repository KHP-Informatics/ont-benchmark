#!/usr/bin/env nextflow

process QUERY_RSID_POSITIONS {
    input:
    path rsids_file

    output:
    path "rsid_positions.txt", emit: rsid_positions

    secret 'NCBI_EMAIL'
    secret 'NCBI_API_KEY'

    script:
    """
    #!/usr/bin/env python -u
    import os
    import time
    from Bio import Entrez
    import xml.etree.ElementTree as ET
    from urllib.error import HTTPError

    Entrez.email = os.environ["NCBI_EMAIL"]
    Entrez.api_key = os.environ["NCBI_API_KEY"]


    def fetch_positions(rsids):
        print(f"Fetching positions for {len(rsids)} rsIDs...")
        query = " OR ".join(f"{rsid}[Reference SNP ID]" for rsid in rsids)

        max_retries = 3
        retry_delay = 5  # seconds

        for attempt in range(max_retries):
            try:
                handle = Entrez.esearch(db="snp", term=query, retmax=10000)
                record = Entrez.read(handle)
                handle.close()

                if not record["IdList"]:
                    print("No results found for this batch.")
                    return {}

                print(f"Found {len(record['IdList'])} matching SNPs. Fetching details...")
                handle = Entrez.efetch(
                    db="snp", id=",".join(record["IdList"]), retmode="xml"
                )
                xml_data = handle.read()
                handle.close()

                positions = {}
                root = ET.fromstring(xml_data)

                ns = {"ns": "https://www.ncbi.nlm.nih.gov/SNP/docsum"}

                for doc_sum in root.findall(".//ns:DocumentSummary", ns):
                    snp_id = doc_sum.find("ns:SNP_ID", ns)
                    chrpos = doc_sum.find("ns:CHRPOS", ns)
                    if snp_id is not None and chrpos is not None and chrpos.text:
                        rs_id = f"rs{snp_id.text}"
                        positions[rs_id] = chrpos.text
                print(f"{len(positions)} positions found in this batch")

                return positions

            except HTTPError as e:
                if e.code == 400 and attempt < max_retries - 1:
                    print(f"Received HTTP 400 error. Retrying in {retry_delay} seconds...")
                    time.sleep(retry_delay)
                else:
                    raise

        return {}


    print(f"Reading rsIDs from file: ${rsids_file}")
    with open("${rsids_file}", "r") as f:
        rsids = [line.strip() for line in f if line.strip()]
    print(f"Total rsIDs to process: {len(rsids)}")

    batch_size = 1000
    all_positions = {}

    for i in range(0, len(rsids), batch_size):
        batch = rsids[i : i + batch_size]
        print(f"Processing batch {i//batch_size + 1} of {-(-len(rsids)//batch_size)}...")
        positions = fetch_positions(batch)
        all_positions.update(positions)
        print(f"Cumulative positions retrieved: {len(all_positions)}")
        time.sleep(1)

    print("Writing results to file...")
    with open("rsid_positions.txt", "w") as f:
        for rsid, position in all_positions.items():
            f.write(f"{rsid}\\t{position}\\n")
    print(f"{len(all_positions)} positions written")

    missing_rsids = set(rsids) - set(all_positions.keys())
    if missing_rsids:
        print(f"Warning: The following rsIDs were not retrieved: {missing_rsids}")
    """
}

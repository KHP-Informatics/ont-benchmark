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
    import concurrent.futures
    import threading
    from Bio import Entrez
    import xml.etree.ElementTree as ET
    from urllib.error import HTTPError
    from http.client import IncompleteRead

    Entrez.email = os.getenv("NCBI_EMAIL")
    Entrez.api_key = os.getenv("NCBI_API_KEY")

    api_semaphore = threading.Semaphore(10)


    def api_limiter():
        while True:
            api_semaphore.acquire()
            threading.Timer(0.1, api_semaphore.release).start()


    limiter_thread = threading.Thread(target=api_limiter, daemon=True)
    limiter_thread.start()


    def fetch_positions(rsids):
        with api_semaphore:  # Control API request rate
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
                        continue

                    print(
                        f"Found {len(record['IdList'])} matching SNPs. Fetching details..."
                    )
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

                except (HTTPError, IncompleteRead) as e:
                    sleep_time = retry_delay * (2**attempt)  # Exponential backoff
                    print(
                        f"{type(e).__name__} encountered: {e} - Attempt {attempt+1}/{max_retries}, retrying in {sleep_time} seconds..."
                    )
                    time.sleep(sleep_time)

            return {}


    print(f"Reading rsIDs from file: ${rsids_file}")
    with open("${rsids_file}", "r") as f:
        rsids = [line.strip() for line in f if line.strip()]
    print(f"Total rsIDs to process: {len(rsids)}")

    batch_size = 500
    all_positions = {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=${task.cpus}) as executor:
        future_to_batch = {
            executor.submit(fetch_positions, rsids[i : i + batch_size]): i
            for i in range(0, len(rsids), batch_size)
        }

        for future in concurrent.futures.as_completed(future_to_batch):
            batch_start = future_to_batch[future]
            try:
                positions = future.result()
                all_positions.update(positions)
                print(
                    f"Cumulative positions retrieved: {len(all_positions)} after processing batch starting at {batch_start}"
                )
            except Exception as e:
                print(f"Failed to fetch positions for batch starting at {batch_start}: {e}")

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

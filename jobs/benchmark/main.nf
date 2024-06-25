#!/usr/bin/env nextflow

// Input channels
Channel
    .fromPath(params.sample_ids_file)
    .splitCsv(header:true)
    .map { row ->
        def ont_id = row.ont_id.trim()
        def LP_id = row.LP_id.trim()
        def microarray_vcf = file("${params.microarray_dir}/FinalReport_InfiniumOmni2-5-8v1-4_${LP_id}.vcf.gz")
        def ont_vcf = file("${params.ont_dir}/${ont_id}_sup/${ont_id}_sup.wf_snp.vcf.gz")
        if (microarray_vcf.exists() && ont_vcf.exists()) {
            return tuple(ont_id, LP_id, microarray_vcf, ont_vcf)
        } else {
            return null
        }
    }
    .filter { it != null }
    .set { input_files_ch }

input_files_ch
    .map { it[2] }
    .collect()
    .set { microarray_vcfs_ch }

Channel
    .fromPath(params.reference_fasta)
    .set { reference_fasta_ch}

Channel
    .fromPath(params.array_positions_file)
    .set { array_positions_ch}


// Processes
process CREATE_SDF {
    input:
    path reference_fasta

    output:
    path "${reference_base}.sdf", emit: reference_sdf

    script:
    reference_base = reference_fasta.baseName
    """
    rtg format \
        --format=fasta \
        --output ${reference_base}.sdf \
        ${reference_fasta}
    """
}

process FETCH_ARRAY_POSITIONS {
    input:
    path array_positions_file
    path microarray_vcfs

    output:
    path "array_positions.json", emit: array_positions_json
    path "fetch_array_positions.log", emit: log_file

    script:
    """
    #!/usr/bin/env python3
    import json
    import pysam
    import requests
    import logging
    from collections import defaultdict
    from pathlib import Path

    logging.basicConfig(
        filename="fetch_array_positions.log",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )


    def fetch_rsid_positions(rsids):
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        positions = {}
        for rsid in rsids:
            if rsid.startswith("rs"):
                params = {"db": "snp", "id": rsid[2:], "retmode": "json"}
                response = requests.get(base_url, params=params)
                if response.status_code == 200:
                    data = response.json()
                    if "result" in data and rsid[2:] in data["result"]:
                        snp_data = data["result"][rsid[2:]]
                        assembly = snp_data.get("assembly", "")
                        if assembly == "GRCh38":
                            chrom = snp_data.get("chrpos", "").split(":")[0]
                            pos = snp_data.get("chrpos", "").split(":")[1]
                            if chrom and pos:
                                positions[rsid] = int(pos)
                                logging.info(f"Fetched position for {rsid}: {chrom}:{pos}")
                        else:
                            logging.warning(f"Skipping {rsid}: Not in GRCh38 assembly")
                else:
                    logging.error(
                        f"Failed to fetch data for {rsid}: HTTP {response.status_code}"
                    )
        return positions


    logging.info("Starting FETCH_ARRAY_POSITIONS process")

    array_positions = {}
    with open("${array_positions_file}", "r") as f:
        next(f)  # Skip header
        for line in f:
            name, rsid = line.strip().split("\t")
            if rsid != ".":
                array_positions[name] = rsid
    logging.info(f"Read {len(array_positions)} entries from array positions file")

    all_rsids = set()
    vcf_files = Path("microarray_vcfs").glob("*.vcf.gz")
    for vcf_file in vcf_files:
        logging.info(f"Processing VCF file: {vcf_file}")
        with pysam.VariantFile(str(vcf_file)) as vcf:
            for record in vcf:
                variant_id = record.id
                rsid = array_positions.get(variant_id, variant_id)
                if rsid.startswith("rs"):
                    all_rsids.add(rsid)
    logging.info(f"Collected {len(all_rsids)} unique rsIDs from all VCF files")

    positions_cache = fetch_rsid_positions(all_rsids)
    logging.info(f"Fetched positions for {len(positions_cache)} rsIDs")

    with open("array_positions.json", "w") as f:
        json.dump(positions_cache, f)
    logging.info("Saved positions to array_positions.json")

    logging.info("FETCH_ARRAY_POSITIONS process completed")
    """
}

// Workflow
workflow {
    CREATE_SDF(reference_fasta_ch)
    FETCH_ARRAY_POSITIONS(array_positions_ch, microarray_vcfs_ch)
    FETCH_ARRAY_POSITIONS.out.array_positions_json.view()
    FETCH_ARRAY_POSITIONS.out.log_file.view()
}

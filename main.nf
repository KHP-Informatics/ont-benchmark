#!/usr/bin/env nextflow

// Input channels
Channel
    .fromPath(params.sample_ids_file)
    .splitCsv(header:true)
    .map { row ->
        def ont_id = row.ont_id.trim()
        def LP_id = row.LP_id.trim()
        def microarray_vcf = file("${params.microarray_dir}/FinalReport_InfiniumOmni2-5-8v1-4_${LP_id}.vcf.gz")
        def microarray_index = file("${params.microarray_dir}/FinalReport_InfiniumOmni2-5-8v1-4_${LP_id}.vcf.gz.tbi")
        def ont_vcf = file("${params.ont_dir}/${ont_id}_sup/${ont_id}_sup.wf_snp.vcf.gz")
        def ont_index = file("${params.ont_dir}/${ont_id}_sup/${ont_id}_sup.wf_snp.vcf.gz.tbi")
        if (microarray_vcf.exists() && microarray_index.exists() && ont_vcf.exists() && ont_index.exists()) {
            return tuple(ont_id, LP_id, microarray_vcf, microarray_index, ont_vcf, ont_index)
        } else {
            return null
        }
    }
    .filter { it != null }
    .set { input_files_ch }

input_files_ch
    .map { it[2..3] }
    .collect()
    .map { vcfs ->
        def vcf_files = vcfs.findAll { it.name.endsWith('.vcf.gz') }
        def index_files = vcfs.findAll { it.name.endsWith('.vcf.gz.tbi') }
        return tuple(vcf_files, index_files)
    }
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
    tuple path(microarray_vcfs), path(microarray_indexes)
    val dbsnp_build

    output:
    path "array_positions.json", emit: array_positions_json
    path "fetch_array_positions.log", emit: log_file

    secret "NCBI_EMAIL"
    secret "NCBI_API_KEY"

    script:
    """
    #!/usr/bin/env python3
    import json
    import pysam
    import logging
    from Bio import Entrez

    logging.basicConfig(
        filename="fetch_array_positions.log",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    Entrez.email = "\$NCBI_EMAIL"
    Entrez.api_key = "\$NCBI_API_KEY"


    def process_vcf_file(vcf_file):
        rsids = set()
        logging.info(f"Processing VCF file: {vcf_file}")
        try:
            with pysam.VariantFile(str(vcf_file)) as vcf:
                record_count = 0
                for record in vcf:
                    record_count += 1
                    variant_id = record.id
                    rsid = array_positions.get(variant_id, variant_id)
                    if rsid.startswith("rs"):
                        rsids.add(rsid)
                logging.info(f"Processed {record_count} records in {vcf_file}")
        except Exception as e:
            logging.error(f"Error processing {vcf_file}: {str(e)}")
        return rsids


    def fetch_rsid_positions(rsids, build=None):
        if build is None:
            build = get_latest_dbsnp_build()

        positions = {}
        for rsid in rsids:
            if rsid.startswith("rs"):
                try:
                    handle = Entrez.esummary(db="snp", id=rsid[2:], version=str(build))
                    record = Entrez.read(handle)
                    if record:
                        snp_data = record[0]
                        chrpos = snp_data.get("CHRPOS")
                        upd_build = int(snp_data.get("UPD_BUILD", 0))
                        if chrpos and upd_build == build:
                            chrom, pos = chrpos.split(":")
                            positions[rsid] = int(pos)
                            logging.info(f"Fetched position for {rsid}: {chrom}:{pos}")
                        else:
                            logging.warning(
                                f"Skipping {rsid}: Not in build {build} or no position information"
                            )
                    else:
                        logging.warning(f"No data found for {rsid} in build {build}")
                except Exception as e:
                    logging.error(f"Failed to fetch data for {rsid}: {str(e)}")
        return positions


    def get_latest_dbsnp_build():
        try:
            handle = Entrez.einfo(db="snp")
            record = Entrez.read(handle)
            latest_build = int(record["DbInfo"]["LastMajorReleaseNumber"])
            logging.info(f"Using latest dbSNP build: {latest_build}")
            return latest_build
        except Exception as e:
            logging.error(f"Failed to fetch latest dbSNP build: {str(e)}")
            logging.warning("Defaulting to build 156")
            return 156


    array_positions = {}
    with open("${array_positions_file}", "r") as f:
        next(f)  # Skip header
        for line in f:
            name, rsid = line.strip().split("\t")
            if rsid != ".":
                array_positions[name] = rsid
    logging.info(f"Read {len(array_positions)} entries from array positions file")

    vcf_files = "${microarray_vcfs}".split()
    logging.info(f"Found {len(vcf_files)} VCF files")

    all_rsids = set()
    for vcf_file in vcf_files:
        file_rsids = process_vcf_file(vcf_file)
        all_rsids.update(file_rsids)

    logging.info(f"Collected {len(all_rsids)} unique rsIDs from all VCF files")

    positions_cache = fetch_rsid_positions(
        all_rsids, build=params.dbsnp_build if params.dbsnp_build else None
    )
    logging.info(f"Fetched positions for {len(positions_cache)} rsIDs")

    with open("array_positions.json", "w") as f:
        json.dump(positions_cache, f)
    logging.info("Saved positions to array_positions.json")
    """
}

// Workflow
workflow {
    CREATE_SDF(
        reference_fasta_ch
)

    FETCH_ARRAY_POSITIONS(
        array_positions_ch,
        microarray_vcfs_ch,
        params.dbsnp_build
)
}

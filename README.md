# Investigating the performance of Oxford Nanopore long-read sequencing with respect to Illumina microarrays and short-read sequencing

[![Nextflow](https://img.shields.io/badge/nextflow-24.10.2-brightgreen.svg?style=for-the-badge&logo=Nextflow)](https://www.nextflow.io/) [![Jupyter Badge](https://img.shields.io/badge/Jupyter-F37626?style=for-the-badge&logo=jupyter&logoColor=white)](https://jupyter.org/)[![Python Language Badge](https://img.shields.io/badge/Python-3.12-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)[![Code style: Black](https://img.shields.io/badge/code%20style-black-000000.svg?style=for-the-badge)](https://black.readthedocs.io/en/stable/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](LICENSE.md)

This repository contains the complete workflow and analysis scripts for benchmarking Oxford Nanopore Technologies (ONT) long-read sequencing against established platforms (Illumina short-read sequencing and microarrays). The project evaluates the performance of ONT for detecting various genetic variants across different genomic contexts and examines the impact of experimental factors such as multiplexing, sequencing depth, and read length.

## Abstract

Oxford Nanopore Technologies (ONT) long-read sequencing (LRS) has emerged as a promising genomic analysis tool, yet comprehensive benchmarks with established platforms across diverse datasets remain limited. This study aimed to benchmark LRS performance against Illumina short-read sequencing (SRS) and microarrays for variant detection across different genomic contexts and to evaluate the impact of experimental factors. We sequenced 14 human genomes using the three platforms and evaluated single nucleotide variants (SNVs), insertions/deletions (indels), and structural variants (SVs) detection, stratifying by high-complexity, low-complexity, and dark genome regions while assessing effects of multiplexing, depth, and read length. LRS SNV accuracy was slightly lower than that of SRS in high-complexity regions (F-measure: 0.954 vs. 0.967) but showed comparable sensitivity in low-complexity regions. LRS showed robust performance for small (1–5 bp) indels in high-complexity regions (F-measure: 0.869), but SRS agreement decreased significantly in low-complexity regions and for larger indel sizes. Within dark regions, LRS identified more indels than SRS, but showed lower base-level accuracy. LRS identified 2.86 times more SVs than SRS, excelling at detecting large variants (>6 kb), with SV detection improving with sequencing depth. Sequencing depth strongly influenced variant calling performance, whereas multiplexing effects were minimal. Our findings provide valuable insights for optimising LRS applications in genomic research and diagnostics.

## Project Structure

```bash
.
├── config/ # Workflow configuration files
├── jobs/ # Slurm job submission scripts
│ ├── benchmark/ # Benchmarking scripts
│ ├── illumina/ # Illumina processing scripts
│ ├── jupyter/ # Jupyter notebook environment setup
│ ├── ont/ # ONT processing scripts
│ └── qc/ # Quality control scripts
├── modules/ # Nextflow modules
│ ├── indel_benchmark/ # Indel analysis modules
│ ├── setup/ # Data preparation modules
│ ├── shared/ # Shared utility modules
│ ├── snv_benchmark/ # SNV analysis modules
│ └── sv_consensus/ # Structural variant consensus modules
├── references/ # Genome and positional reference files
├── workflows/ # Nextflow sub-workflows
├── main.nf # Main Nextflow workflow
├── nextflow.config # Nextflow configuration
└── ont-benchmark.ipynb # Jupyter notebook with statistical analyses
└── sample_ids.csv # ONT and Illumina sample IDs dictionary
└── seq_stats.csv # Table containing experimental records for each flowcell
```

## Setup

### Prerequisites

#### Nextflow Pipeline

- [Nextflow](https://www.nextflow.io/) (24.10.2)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)

#### Jupyter Notebook

- Please see the [conda environment](jobs/jupyter/environment.yml) for software requirements and dependencies

### Data Requirements

The analysis pipeline expects:

1. Oxford Nanopore sequencing data (processed through basecalling)
1. Illumina short-read sequencing data (aligned and variant-called)
1. Illumina microarray genotyping data

### NCBI API Key Setup

To ensure the pipeline functions correctly and to optimize access to NCBI resources, please set the following Nextflow secrets before running the workflow:

```bash
nextflow secrets set NCBI_API_KEY <your_ncbi_api_key>
nextflow secrets set NCBI_EMAIL <your_ncbi_email>
```

These secrets are necessary for accessing NCBI resources during the analysis. By default, the NCBI Datasets API and command-line tool requests are rate-limited to 5 requests per second (rps). Using an API key increases this limit to 10 rps.

For more information on obtaining and using NCBI API keys, please refer to the [NCBI Datasets API Keys Documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/api-keys/).

You can verify that the secrets have been set correctly by listing them:

```bash
nextflow secrets list
```

For more information on managing secrets in Nextflow, refer to the [Nextflow Secrets documentation](https://www.nextflow.io/docs/stable/secrets.html).

## Usage

### Running the Complete Workflow

```bash
nextflow run renatosantos98/ont-benchmark
```

or

```bash
sbatch jobs/benchmark/variant_benchmark.sh
```

# Results

Analysis results are stored in the [ont-benchmark](ont-benchmark.ipynb) jupyter notebook, organised by variant type. Each benchmark includes:

1. Precision, recall, and F-measure metrics
1. Detailed comparison between ONT, Illumina, and microarray platforms
1. Analysis of variant detection across different genomic contexts
1. Impact assessment of sequencing parameters (depth, multiplexing, read length)

# License

This project is licensed under the MIT License. You can freely use and modify the code, without warranty. See [LICENSE](LICENSE) for the full license text. The authors reserve the rights to the article content, which is currently submitted for publication.

# Citation

If you use this benchmark in your research, please cite:
[Santos R, Lee H, Williams A, Baffour-Kyei A, Breen G, Iacoangeli A. Investigating the performance of Oxford Nanopore long-read sequencing with respect to Illumina microarrays and short-read sequencing. bioRxiv. 2024 Dec 22:2024-12.](https://doi.org/10.1101/2024.12.19.629409)

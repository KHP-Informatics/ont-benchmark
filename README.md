# Investigating the performance of Oxford Nanopore long-read sequencing with respect to Illumina microarrays and short-read sequencing

[![Nextflow](https://img.shields.io/badge/nextflow-24.10.2-brightgreen.svg?style=for-the-badge&logo=Nextflow)](https://www.nextflow.io/) [![Jupyter Badge](https://img.shields.io/badge/Jupyter-F37626?style=for-the-badge&logo=jupyter&logoColor=white)](https://jupyter.org/)[![Python Language Badge](https://img.shields.io/badge/Python-3.12-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)[![Code style: Black](https://img.shields.io/badge/code%20style-black-000000.svg?style=for-the-badge)](https://black.readthedocs.io/en/stable/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](LICENSE.md)

This repository contains the complete workflow and analysis scripts for benchmarking Oxford Nanopore Technologies (ONT) long-read sequencing against established platforms (Illumina short-read sequencing and microarrays). The project evaluates the performance of ONT for detecting various genetic variants across different genomic contexts and examines the impact of experimental factors such as multiplexing, sequencing depth, and read length.

## Abstract

Oxford Nanopore Technologies (ONT) long-read sequencing (LRS) has emerged as a promising tool for genomic analysis, but comprehensive comparisons with established platforms across diverse datasets remain limited. We present a multi-platform benchmark using 14 human genomes sequenced with ONT LRS, Illumina short-read sequencing (SRS), and Illumina microarrays. Our study evaluates LRS performance for various genetic variants across genomic contexts, while also examining the impact of experimental factors such as multiplexing, depth, and read length.

In high-complexity regions, LRS demonstrated competitive yet slightly lower accuracy than SRS for SNV detection (F-measure: 0.954 vs. 0.968), with performance gaps narrowing in low-complexity regions. For indel detection, LRS showed robust performance for small indels (1-5bp) in high-complexity regions (F-measure: 0.869), but accuracy decreased significantly in low-complexity regions and for larger indels. LRS identified 2.86 times more structural variants than SRS, with superior detection of large-scale variations.

Sequencing depth strongly influenced variant calling performance across all variant types, while multiplexing effects were minimal after controlling for depth. Our findings provide valuable insights for optimising ONT LRS applications in genomic research and clinical diagnostics.

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

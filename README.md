# Bacterial WGS Pipeline

A comprehensive Nextflow pipeline for bacterial whole genome sequencing (WGS) analysis. This pipeline performs quality control, read trimming, and genome assembly for short-read, long-read, or hybrid sequencing data.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Format](#input-format)
- [Usage](#usage)
- [Pipeline Workflow](#pipeline-workflow)
- [Output](#output)
- [Parameters](#parameters)
- [Tools](#tools)
- [Troubleshooting](#troubleshooting)

## Overview

This pipeline is designed to handle bacterial genome assembly from various sequencing platforms:

- **Short-read assembly**: Uses SKESA for Illumina/short-read data
- **Long-read assembly**: Uses HYBRACTER for Oxford Nanopore or PacBio data
- **Hybrid assembly**: Combines short and long reads using HYBRACTER for improved assembly quality

The pipeline includes comprehensive quality control with FastQC for short reads, LongQC for long reads, and adaptive trimming with FASTP.

## Features

✅ **Flexible input handling**: Supports short reads only, long reads only, or hybrid sequencing data
✅ **Pre-QC analysis**: FastQC and LongQC quality assessment before trimming
✅ **Adaptive trimming**: FASTP intelligently adjusts quality thresholds based on input quality
✅ **Post-QC validation**: Quality assessment after trimming for QC tracking
✅ **MultiQC reports**: Comprehensive summary reports for all QC steps
✅ **Multiple assembly strategies**: Automatic selection of SKESA or HYBRACTER based on input type
✅ **Containerized tools**: All tools run in Docker containers for reproducibility
✅ **Scalable**: Fully parallelized analysis with support for multiple samples

## Requirements

### Software
- **Nextflow** >= 21.0.0
- **Docker** or **Singularity** (for container support)

### System Resources
Recommended specifications:
- CPU: 8+ cores (scalable)
- Memory: 32+ GB RAM
- Storage: Sufficient space for raw reads and assemblies (highly data-dependent)

## Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd bacterial-wgs-pipeline
```

2. Ensure Nextflow is installed:
```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
```

3. Verify Docker/Singularity installation for running containerized tools.

## Quick Start

1. Prepare a samplesheet CSV file (see [Input Format](#input-format))

2. Run the pipeline:
```bash
nextflow run main.nf \
    --samplesheet samples.csv \
    --outdir results
```

## Input Format

Create a `samplesheet.csv` file with the following columns:

```csv
sample_id,fastq_1,fastq_2,long_reads
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,
sample2,,/path/to/sample2_long.fastq.gz
sample3,/path/to/sample3_R1.fastq.gz,/path/to/sample3_R2.fastq.gz,/path/to/sample3_long.fastq.gz
```

**Column descriptions:**
- `sample_id`: Unique sample identifier
- `fastq_1`: Path to R1 reads (Illumina forward reads)
- `fastq_2`: Path to R2 reads (Illumina reverse reads)
- `long_reads`: Path to long-read file (Oxford Nanopore or PacBio)

**Notes:**
- Either `fastq_1`/`fastq_2` OR `long_reads` OR both must be provided
- Leave blank (empty) for missing read types
- Paths should be absolute or relative to working directory

## Usage

### Basic usage:
```bash
nextflow run main.nf --samplesheet samples.csv --outdir results
```

### With custom output directory:
```bash
nextflow run main.nf \
    --samplesheet samples.csv \
    --outdir /custom/output/path
```

### Resume a failed run:
```bash
nextflow run main.nf --samplesheet samples.csv --outdir results -resume
```

## Pipeline Workflow

```
Input Reads (Samplesheet)
    ↓
[QC Workflow]
    ├─ FastQC (Pre-trimming) → Multi-sample report (MultiQC)
    ├─ LongQC (Long reads validation)
    ├─ FASTP (Adaptive quality-based trimming)
    ├─ FastQC (Post-trimming)
    └─ LongQC (Long reads re-validation)
    ↓
[Assembly Workflow]
    ├─ SKESA (Short-read only samples)
    ├─ HYBRACTER (Long-read only or hybrid samples)
    └─ Assembly contigs output
    ↓
Final Assemblies
```

## Output

The pipeline generates results in the following directory structure:

```
results/
├── fastqc/              # FastQC reports (pre- and post-trimming)
├── longqc/              # LongQC reports
├── fastp/               # Trimming reports and statistics
├── multiqc/             # Aggregated QC reports (pre- and post-trimming)
├── skesa/               # SKESA assemblies (short-read only samples)
└── hybracter/           # HYBRACTER assemblies (long-read or hybrid samples)
```

### Key output files:
- `*_multiqc_report.html` - Interactive MultiQC report
- `*/fastqc_report.html` - Per-sample FastQC reports
- `*.fasta` - Final assembled contigs

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `samplesheet` | string | required | Path to samplesheet CSV file |
| `outdir` | string | `./results` | Output directory path |
| `skesa_cores` | integer | `12` | CPU cores for SKESA |
| `skesa_memory` | integer | `32` | Memory (GB) for SKESA |
| `skesa_min_contig` | integer | `500` | Minimum contig length (bp) for SKESA output |

### Usage examples:
```bash
# Adjust SKESA parameters
nextflow run main.nf \
    --samplesheet samples.csv \
    --skesa_cores 16 \
    --skesa_memory 64 \
    --skesa_min_contig 1000

# Specify output directory
nextflow run main.nf \
    --samplesheet samples.csv \
    --outdir /mnt/results
```

## Tools

The pipeline uses the following bioinformatics tools (all containerized):

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.12.1 | Short-read quality control |
| LongQC | latest | Long-read quality control |
| FASTP | 1.0.1 | Read trimming and filtering |
| MultiQC | latest | Report aggregation |
| SKESA | 2.4.0 | Short-read genome assembly |
| HYBRACTER | latest | Hybrid genome assembly |

## Troubleshooting

### Common issues:

**Error: "Please provide a samplesheet using the '--samplesheet' parameter."**
- Solution: Ensure you pass the `--samplesheet` parameter with a valid file path

**Error: "Docker daemon is not running"**
- Solution: Start Docker service (`docker daemon` or through system settings)

**Out of memory errors**
- Solution: Increase `skesa_memory` parameter or available system RAM

**Pipeline hangs or times out**
- Solution: Check system resources and increase timeouts in config file if needed

### Debugging:

Enable verbose logging:
```bash
nextflow run main.nf --samplesheet samples.csv -v
```

View pipeline DAG (directed acyclic graph):
```bash
nextflow run main.nf --samplesheet samples.csv -with-dag flowchart.html
```

## Citation

If you use this pipeline, please cite the individual tools used:
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [LongQC](https://github.com/yfukasawa/LongQC)
- [FASTP](https://github.com/OpenGene/fastp)
- [SKESA](https://github.com/ncbi/SKESA)
- [HYBRACTER](https://github.com/gbouras13/hybracter)
- [Nextflow](https://www.nextflow.io/)

## License

[Add your license here]

## Support

For issues or questions, please open an issue on the repository.

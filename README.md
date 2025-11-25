# Collecting variants with snippy

How to run?

`python3 vanilla/vanilla.py -c config/vanilla.config.yaml -s config/vanilla.samples.tsv`


# Bioinformatics Pipeline

This repository contains a comprehensive bioinformatics pipeline for microbial genomics. The pipeline performs read trimming, variant calling, ORF prediction, protein annotation, GenBank file generation, and aggregation of results across multiple samples.

## Features

- **Read trimming** using Cutadapt.
- **Variant calling** using Snippy.
- **ORF prediction** with EMBOSS `getorf`.
- **Protein annotation** using DIAMOND.
- **GenBank file creation** for predicted ORFs.
- **Result concatenation** for downstream analysis.
- **Parallel processing** using Python's multiprocessing.

## Requirements

- Python 3
- BioPython
- pandas
- PyYAML
- Cutadapt
- Snippy
- EMBOSS (getorf)
- DIAMOND
- gzip

## Installation

Clone this repository and ensure the required tools are installed:

```bash
git clone <repository_url>
cd <repository>
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Ensure external tools (`cutadapt`, `snippy`, `getorf`, `diamond`) are installed and available in your `$PATH`.

## Configuration

The pipeline requires a YAML configuration file specifying paths and parameters. Example:

```yaml
project: my_project
results: Variant_results
threads: 40
minithreads: 1
fastq: data/fastq
trimmed: Intermediate_files/trimmed
references: Intermediate_files/References_for_variants
snippy: Intermediate_files/Snippy
database: Databases/nr.dmnd
variants: Variant_results/all_variants.tsv
force: false
```

## Samples Table

Provide a TSV file listing the samples and their reference genomes:

```tsv
Sample	Reference_name	Reference_file
sample1	ref1	references/ref1.fna
sample2	ref2	references/ref2.fna
```

## Usage

Run the pipeline with:

```bash
python3 vanilla/bin/vanilla.py -c config/config.yaml -s config/samples.tsv
```

The pipeline will:

1. Process reference genomes (ORF prediction, DIAMOND annotation, BED and GenBank generation).
2. Trim reads for each sample.
3. Perform variant calling using Snippy.
4. Aggregate all variant results into a single TSV.

## Output Structure

- `results/` – Final results including concatenated variants.
- `references/` – Reference ORF predictions, BED files, and GenBank files.
- `snippy_results/` – Individual Snippy outputs per sample.
- `data/trimmed/` – Trimmed FASTQ files.

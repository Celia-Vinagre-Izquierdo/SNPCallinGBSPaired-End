# SNP Calling for GBS paired-end (PE) data using TASSEL

Scripts for SNPCalling GBS paired-end (PE) data and karyotype classification.

## Table of Contents
- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Contact](#contact)

## Introduction
This project contains scripts for processing GBS paired-end data and classifying karyotypes for quails. This script combines bash and R programming

## Dependencies
The following programs are required to run the scripts in this repository:

- [axe-demux](https://github.com/username/axe-demux) 
- [trim_galore](https://github.com/FelixKrueger/TrimGalore)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [bbmerge](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [TASSEL](https://www.maizegenetics.net/tassel)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [vcftools](https://vcftools.github.io/)
- [plink](https://www.cog-genomics.org/plink/)

## Workflow

### Demultiplex raw sequence data
- **axe-demux**: For demultiplexing raw sequence data.

### Remove adapters
- **trim_galore**: For removing adapters and processing paired-end reads.
- **FastQC**: For quality control of trimmed files.

### Merge paired-end reads
- **bbmerge.sh**: For merging paired-end reads.
- **FastQC**: For quality control of merged reads.

### Select >64bp with cutadapt
- **cutadapt**: For trimming adapters and selecting reads longer than 64bp.
- **FastQC**: For quality control of trimmed samples.

### Prepare your files for TASSEL
- **Custom R function**: For adding fake barcodes to FASTQ files for TASSEL5 GBSv2 compatibility.

### Use TASSEL
- **wget**: For downloading the reference genome.
- **gunzip**: For unzipping the reference genome.
- **conda activate**: For activating the conda environment.
- **bowtie2-build**: For building the reference genome index.
- **TASSEL**: For various GBS analysis steps including `GBSSeqToTagDBPlugin`, `TagExportToFastqPlugin`, `SAMToGBSdbPlugin`, `DiscoverySNPCallerPluginV2`, `ProductionSNPCallerPluginV2`, `GetTagTaxaDistFromDBPlugin`, and `genotypeSummary`.
- **bowtie2**: For aligning sequences to the reference genome.
- **gzip**: For compressing output files.

### Filter VCF and generate PCA
- **vcftools**: For filtering VCF files.
- **plink**: For performing PCA analysis.

## Contact

email: celia.vinagre.izquierdo@gmail.com

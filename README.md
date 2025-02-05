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

### To prepare your raw data to use TASSEL 
- [axe-demux](https://github.com/username/axe-demux): For demultiplexing raw sequence data.
- [trim_galore](https://github.com/FelixKrueger/TrimGalore): For removing adapters and processing paired-end reads.
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): For quality control of sequence data.
- [bbmerge.sh](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/): For merging paired-end reads for GBS data.
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/): For trimming adapters and selecting reads longer than 64bp.

After this you need to use the fake_barcode.R script and then continue in bash. 

### Once you have your sequences ready for using TASSEL to create a VCF file
- [TASSEL](https://www.maizegenetics.net/tassel): For GBS analysis.
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): For aligning sequences to the reference genome.

### To filter your VCF file
- [vcftools](https://vcftools.github.io/): For filtering VCF files.
- [plink](https://www.cog-genomics.org/plink/): For performing PCA analysis.

## Contact

email: celia.vinagre.izquierdo@gmail.com

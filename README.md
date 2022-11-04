# DUPFinder: Detection tools of Duplicates gene 

DUPfinder est un outils qui permet de caractèriser les gènes dupliqués en se basant sur la détection des variations structurales (CNV) et particulièrement des duplications dà partir un échantillon ou d'une population d'échantillons.

## Table of contents

- [Workflow of Hecaton](#overview)
- [Installation](#install)
  - [Prerequisites](#prerequisites)
  - [Gettind installation](#Getting)
- [Usage](#usage)
  - [Running multiple samples at once](#multiple)





## Workflow of DUPFinder

* Aligning reads to a reference genome using [bwa mem][bwa mem]
* Calling CNVs using the structural variant callers [Delly][delly], [Dysgu][Dysgu] and [Delly][delly]
* Post-processing each set of CNVs to keep the duplications and remove false positives [Duphold][duphold], [Bcftools][bcftools]
* Merging all sets of duplications into one large set [SURVIVOR][Survivor]
* Detection of duplication gene using the annotation file [Bedtools][bedtools]

<p align="center">
    <img title="DUPFinder Workflow" src="image/dupfinder_workflow.png" width=60%>
</p>

## Installation

### Prerequisites

DUPFinder can only be installed on Linux systems and requires Anaconda/Miniconda (Python 3.9+) to be present on the system.

/* `python3-dev`, `git`, `gawk`, `python>=3.9`, `delly`, `smoove`, `duphold`, ` SURVIVOR`, ` bcftools`, `numpy`, ` pysam`, ` bwa`, ` samtools>=1.10`, `htslib>=1.10`, `bedtools`, `pip`, `dysgu` */

All steps of Hecaton are run using the [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10`) workflow language.

## Getting Started
## <a name="install"></a>Quick installation using conda 

```bash

#Step 1. Download the dUPFinder :

git clone https://github.com/assane/dupfinder.git

#Step 2. Go to dupfinder folder

cd dupfinder

#Step 3. Find the yaml file in the folder and run :

conda env create -f dupfinder_env.yml

bash install.sh

#Step 4. Activate the environnement dupfinder_env:

conda activate dupfinder_env


## Index
```bash
# before starting, create index file for the reference genome to reduce mapping time.
 
# build index accordingly

bwa index reference.fa 

```
## <a name="Usage"></a>Usage

```
DUPfinder: Tool for detecting duplicate gene using Illumina sequencing data.

  DUPFinder version: v1.0.0
 
    Usage: nextflow run dupfinder.nf --c file.config --genome_file reference.fa --reads "pair_id_{1,2}.fastq" --annot file.bed --out Output_DUPFinder

    Command arguments DUPFinder: The following parameters need to be specified when running DUPFinder
    
	    --genome_file       Reference genome in FASTA format
	    
	    --reads             set of paired-end reads in FASTQ format. Gzipped FASTQ files are allowed
	    
	    --annot             the file containing the gene annotation: it can be in gff or bed format and must be tabulated
	    
	    --out               Output directory to which all results will be written
	    
	    --c                 Config file specifying the number of CPU cores and memory that will be assigned to DUPFinder
	   	    
	   	    
    Optional arguments:
	    -w               Working directory to which intermediate results will be written. Default: work
            -v               version
```

### <a name="multiple"></a>Running multiple samples at once
DUPFinder can be used to run multiple samples using a single command. For exemple if there existe several sample paired-end, they can all be processed using:
```bash
nextflow run dupfinder.nf --c file.config --genome_file reference.fa --reads "*_{1,2}.fastq" --annot file.bed --out Output_DUPFinder
```


















#!/bin/bash

set -e

nextflow run dupfinder.nf \
--c nextflow.config \
--genome_file test/Arabido.fasta \
--reads "test/sample_{R1,R2}.fastq.gz" \
--annot test/Araport11_gene.gff.gz \
--out test_out \
-w work

echo "test successful: Good Analysis with dupfinder pipeline V1.0.0"

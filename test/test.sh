#!/bin/bash
set -e 

mkdir -p test_sample && cd test_sample && \

nextflow run ../dupfinder/dupfinder.nf \
--c ../dupfinder/nextflow.config \
--genome_file ../test_data/Arabido.fasta \
--reads "../test_data/sample_{R1,R2}.fastq.gz" \
--annot ../test_data/Araport11_gene.gff.gz \
--out test_out \
-w test_work

echo "test successful: Goog Analyze with dupfinder pipeline"

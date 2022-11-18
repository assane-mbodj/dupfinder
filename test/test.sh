#!/bin/bash
set -e 

mkdir -p test_sample && cd test_sample && \

nextflow run ../dupfinder/dupfinder.nf \
--c ../dupfinder/nextflow.config \
--genome_file ../test/Arabido.fasta \
--reads "../test/sample_{R1,R2}.fastq.gz" \
--annot ../test/Araport11_gene.gff.gz \
--out test_out \
-w test_work

echo "test successful: Goog Analyze with dupfinder pipeline"

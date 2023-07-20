/*
*
*/

nextflow.enable.dsl=2

params.genome_file = ""
params.reads_lr = "{prefix}*.f*q*"
params.annot = ""
params.out = ""
params.help = false



include { MAPPING_LR } from './modules_lr/minimap2.nf'
include { VARIANT } from './modules_lr/calling_lr.nf'
include { VCF_MERGE } from './modules_lr/merge_vcf_lr.nf'
include { GENE_DUPLICATED } from './modules_lr/dup_gene_lr.nf'
include { VCF_FILTER } from './modules_lr/filter_vcf_lr.nf'



Channel
	.fromFilePath(params.reads_lr, checkIfExists: true )
	.ifEmpty { error "Cannot find any fastq files matching: ${params.reads}" }
	.set{reads_lr}
	
	
//Declare your variable
genome = file(params.genome_file)
genome_index_file = file(params.genome_file + ".fai")
annotation = file(params.annot)



workflow {

	MAPPING_LR ( reads_lr, genome)
	
	}

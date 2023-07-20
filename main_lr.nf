/*
*
*/

nextflow.enable.dsl=2

params.genome_file = ""
params.reads_lr = "{prefix}*.f*q*"
params.annot = ""
params.out = ""
params.help = false



include { MAPPING_LR } from './module_lr/minimap2.nf'
include { VARIANT_LR } from './module_lr/calling_lr.nf'
include { FILTER_LR } from './module_lr/filter_lr.nf'
include { MERGE_LR } from './module_lr/mergeVcf_lr.nf'
include { GENE_DUPLICATED_LR } from './module_lr/dupgene_lr.nf'


//Declare your variable
genome = file(params.genome_file)
genome_index_file = file(params.genome_file + ".fai")
annotation = file(params.annot)


Channel
	.fromPath(params.reads_lr)
	.ifEmpty { error "Cannot find any fastq files matching: ${params.reads_lr}" }
	.map { file -> tuple( file.baseName, file )}
	.set{reads_lr}
	
	



workflow {

	MAPPING_LR( reads_lr, genome)
	
	VARIANT_LR(MAPPING_LR.out, genome)
	
	FILTER_LR (MAPPING_LR.out.bam_file_ch, VARIANT_LR.out.variant_calls, genome)
	
	MERGE_LR ( FILTER_LR.out.duplication_annot_calls )
	
	GENE_DUPLICATED_LR (MERGE_LR.out.filterVcf_file, annotation)
	
	}
	
	
	/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
       DUPFinder long reads execution summary
       ---------------------------
       Completed at	: ${workflow.complete}
       Duration    	: ${workflow.duration}
       Success     	: ${workflow.success}
       workDir     	: ${workflow.workDir}
       exit status 	: ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}

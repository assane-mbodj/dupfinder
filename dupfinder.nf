/*
 *
 *
 /
 
 
/*
*
*/

nextflow.enable.dsl=2

params.genome_file = ""
params.reads_sr = "{prefix}_*{1,2}.f*q*"
params.annot = ""
params.out = ""
params.reads_lr = "{prefix}_*.f*q*"
params.help = false
params.sr = false
params.lr = false
//params.bam = false


/*
 * Add help message
 */

def helpMessage() {
    log.info"""
    =========================================
     DUPFinder version: v1.0.0
    =========================================
    Usage:
    Usage: nextflow run dupfinder.nf --c file.config --genome_file reference.fa --reads "prefix_{1,2}.fastq" --annot file.bed --out Output_DUPFinder

    Command arguments DUPFinder:
    
	    --genome_file: Reference genome in FASTA format
	    --reads_sr: set of paired-end short reads in FASTQ format. Gzipped FASTQ files are allowed
	    --reads_lr: set of single-end long reads in FASTQ format. Gzipped FASTQ files are allowed
	    --sr: allow to run the short reads version
	    --lr: allow to run the long reads version

	    --annot: the file containing the gene annotation: it can be in gff or bed format and must be tabulated
	    --out: Output directory to which all results will be written
	    --c: Config file specifying the number of CPU cores and memory that will be assigned to DUPFinder
	   	   	    
	    Optional arguments:
	    -w: Working directory to which intermediate results will be written. Default: work

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage

if (params.help){
    helpMessage()
    exit 0
}

/*
 * Input parameters validation and setup
 */

if (! params.out ) {
    println "Missing output directory parameter"
    exit 1
}
if (! params.reads_sr ) {
    println "Missing reads_sr parameter"
    exit 1
 }   
if (! params.reads_lr ) {
    println "Missing reads_lr parameter"
    exit 1
    
}
if (! params.genome_file ) {
    println "Missing genome file parameter"
    exit 1
}



include { MAPPING } from './module_sr/bwa'
include { VARIANT } from './module_sr/calling'
include { VCF_MERGE } from './module_sr/merge_vcf.nf'
include { GENE_DUPLICATED } from './module_sr/dup_gene.nf'
include { VCF_FILTER } from './module_sr/filter_vcf.nf'
//include { DYSGU } from './module_sr/dysgu.nf'
//include { DELLY } from './module_sr/delly.nf'
//include { SMOOVE; DELLY; DYSGU; VCF_MERGE } from './module_sr/smoove.nf'


include { MAPPING_LR } from './module_lr/minimap2.nf'
include { VARIANT_LR } from './module_lr/calling_lr.nf'
include { FILTER_LR } from './module_lr/filter_lr.nf'
include { MERGE_LR } from './module_lr/mergeVcf_lr.nf'
include { GENE_DUPLICATED_LR } from './module_lr/dupgene_lr.nf'




//Declare your variable
genome = file(params.genome_file)
genome_index_file = file(params.genome_file + ".fai")
genome_bwa_amb_file = file(params.genome_file + ".amb")
genome_bwa_ann_file = file(params.genome_file + ".ann")
genome_bwa_bwt_file = file(params.genome_file + ".bwt")
genome_bwa_pac_file = file(params.genome_file + ".pac")
genome_bwa_sa_file = file(params.genome_file + ".sa")
annotation = file(params.annot)



Channel
	.fromFilePairs(params.reads_sr)
//	.ifEmpty { error "Cannot find any fastq files matching: ${params.reads_sr}" }
	.set{reads_file}

Channel
	.fromPath(params.reads_lr)
//	.ifEmpty { error "Cannot find any fastq files matching: ${params.reads_lr}" }
	.map { file -> tuple( file.baseName, file )}
	.set{reads_lr}


workflow {

	if (params.sr) {
		
		MAPPING(reads_file, genome,genome_index_file,genome_bwa_amb_file,genome_bwa_ann_file,genome_bwa_bwt_file, 	genome_bwa_pac_file, genome_bwa_sa_file)
  
  
		VARIANT(MAPPING.out.bam_file_ch, genome, genome_index_file )

		VCF_FILTER(MAPPING.out.bam_file_ch, VARIANT.out.variant_calls, genome, genome_index_file )
	
		VCF_MERGE(VCF_FILTER.out.duplication_annot_calls)

		GENE_DUPLICATED(VCF_MERGE.out.filter_file, annotation)
		
	}
	else if(params.lr) {
		
		MAPPING_LR( reads_lr, genome)
	
		VARIANT_LR(MAPPING_LR.out, genome)
	
		FILTER_LR (MAPPING_LR.out.bam_file_ch, VARIANT_LR.out.variant_calls, genome)
	
		MERGE_LR ( FILTER_LR.out.duplication_annot_calls )
	
		GENE_DUPLICATED_LR (MERGE_LR.out.filterVcf_file, annotation)
		
	//	DYSGU(MAPPING.out.bam_file_ch, genome, genome_index_file )
	//	DELLY(MAPPING.out.bam_file_ch, genome, genome_index_file )
	//	SMOOVE(MAPPING.out.bam_file_ch, genome, genome_index_file )
	
	}
	
}
	
/*
========================================================================================
   Workflow DupFinder tools for long and short reads data
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
       DUPFinder tools execution summary
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
	
	

/*
*
*/

nextflow.enable.dsl=2

params.genome_file = ""
params.reads = "{prefix}_*{1,2}.f*q*"
params.annot = ""
params.out = ""
params.help = false

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
	    --reads: set of paired-end reads in FASTQ format. Gzipped FASTQ files are allowed
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
if (! params.reads ) {
    println "Missing reads parameter"
    exit 1
}
if (! params.genome_file ) {
    println "Missing genome file parameter"
    exit 1
}



include { MAPPING } from './modules/bwa'
include { VARIANT } from './modules/calling'
include { VCF_MERGE } from './modules/merge_vcf.nf'
include { GENE_DUPLICATED } from './modules/dup_gene.nf'
include { VCF_FILTER } from './modules/filter_vcf.nf'
//include { DYSGU } from './modules/dysgu.nf'
//include { DELLY } from './modules/delly.nf'
//include { SMOOVE; DELLY; DYSGU; VCF_MERGE } from './modules/smoove.nf'


Channel
	.fromFilePairs(params.reads, checkIfExists: true )
	.ifEmpty { error "Cannot find any fastq files matching: ${params.reads}" }
	.set{reads_file}


//Declare your variable
genome = file(params.genome_file)
genome_index_file = file(params.genome_file + ".fai")
genome_bwa_amb_file = file(params.genome_file + ".amb")
genome_bwa_ann_file = file(params.genome_file + ".ann")
genome_bwa_bwt_file = file(params.genome_file + ".bwt")
genome_bwa_pac_file = file(params.genome_file + ".pac")
genome_bwa_sa_file = file(params.genome_file + ".sa")
annotation = file(params.annot)






workflow {


	MAPPING(reads_file, genome,genome_index_file,genome_bwa_amb_file,genome_bwa_ann_file,genome_bwa_bwt_file, 	genome_bwa_pac_file, genome_bwa_sa_file)
  
  
	VARIANT(MAPPING.out.bam_file_ch, genome, genome_index_file )

	VCF_FILTER(MAPPING.out.bam_file_ch, VARIANT.out.variant_calls, genome, genome_index_file )
	
	VCF_MERGE(VCF_FILTER.out.duplication_annot_calls)

	GENE_DUPLICATED(VCF_MERGE.out.filter_file, annotation)
	
//VCF_MERGE( DYSGU.out.duplication_annot_calls.mix(DELLY.out.duplication_annot_calls).mix(SMOOVE.out.duplication_annot_calls)).collect()


//	DYSGU(MAPPING.out.bam_file_ch, genome, genome_index_file )
	
//	DELLY(MAPPING.out.bam_file_ch, genome, genome_index_file )
//	SMOOVE(MAPPING.out.bam_file_ch, genome, genome_index_file )
		    

}

/*
========================================================================================
   Workflow Event Handler
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


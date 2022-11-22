/*
 *Run DUPfinder
 *
 *@authors
 *Assane Mbodj <assanembodj11@gmail.com>
 *
 */

params.genome_file = ""
params.reads = "{pair_id}_{1,2}.f*q*"
params.extra_filtering = false
params.out = ""
params.help = false
params.annot = ""

/*
 * Add help message
 */

def helpMessage() {
    log.info"""
    =========================================
     DUPFinder version: v1.0.0
    =========================================
    Usage:
    Usage: nextflow run dupfinder.nf --c file.config --genome_file reference.fa --reads "pair_id_{1,2}.fastq" --annot file.bed --out Output_DUPFinder

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

//Declare your variable

read_files = Channel.fromFilePairs(params.reads)
genome_file = file(params.genome_file)
annotation = file(params.annot)
genome_index_file = file(params.genome_file + ".fai")
genome_bwa_amb_file = file(params.genome_file + ".amb")
genome_bwa_ann_file = file(params.genome_file + ".ann")
genome_bwa_bwt_file = file(params.genome_file + ".bwt")
genome_bwa_pac_file = file(params.genome_file + ".pac")
genome_bwa_sa_file = file(params.genome_file + ".sa")

/*
   Mapping reads on reference genome
*/

process align {
    publishDir "${params.out}/aligned_reads", mode:'copy'
    label "bwa"
    tag "Aligning reads with bwa mem: ${reads[0]} ${reads[1]}"

    input:
	tuple val(pair_id), file(reads) from read_files
	file genome_file
	file genome_bwa_amb_file
	file genome_bwa_ann_file
	file genome_bwa_bwt_file
	file genome_bwa_pac_file
	file genome_bwa_sa_file
	file genome_index_file

    output:
    set val(pair_id), file("${pair_id}*") into aligned_reads_ch

    script:
    """
     
    bwa mem -t ${task.cpus} -R "@RG\\tID:${pair_id}\\tSM:${pair_id}\\tLB:${pair_id}" ${genome_file} ${reads[0]} ${reads[1]} | samtools view -u -bS | samtools sort > ${pair_id}.sort.bam
    samtools index ${pair_id}.sort.bam          
    """    
}

/* 

 * Calling Variants step using three detection tools of structural variants: delly, dysgu anfd smoove

 */

process variant_calling {	
	publishDir "${params.out}/variant_calls", mode:'copy'
	label "multithreading"
	label "delly"
	label "dysgu"
	label "smoove"
	label "bcftools"
	tag "Variants calling: ${pair_id}"
	
	input: 
	tuple val(pair_id), file(reads) from aligned_reads_ch
	file genome_file  

	output:
	tuple val(pair_id), file("*.vcf"), file(reads) into variant_calls_ch
	
	script:
	"""
	
	delly call --map-qual 10 --min-clique-size 1 -g ${genome_file} ${pair_id}.sort.bam -o ${pair_id}.delly.bcf
	
	bcftools view ${pair_id}.delly.bcf > ${pair_id}.delly.vcf
	
	dysgu run --clean --mq 10 --min-size 50 -p4 ${genome_file} temp_dir ${pair_id}.sort.bam > ${pair_id}.dysgu.vcf
	
	smoove call --name ${pair_id} --outdir smoove/ --fasta ${genome_file} -p 4 ${pair_id}.sort.bam 
	
	mv smoove/${pair_id}-smoove.vcf.gz .
	mv ${pair_id}-smoove.vcf.gz ${pair_id}.smoove.vcf.gz
	bgzip -d ${pair_id}.smoove.vcf.gz
	"""
	}
	
/*	
 * Filter and annotation of variants step in keeping only the duplication on the three vcf file
 * Then using survivor to merge vcf file set to get a vcf file contening all structurals variant
 */

process vcf_filter {
	publishDir "${params.out}/duplication_annot_calls", mode:'copy'
	label "multithreading"
	label "bcftools"
	label "duphold"
	label "SURVIVOR"
	tag "Filter vcf file: ${pair_id}"
	
	input:
	tuple val(pair_id), file(vcf), file(bam) from variant_calls_ch
	file genome_file
	
	output:
	file("${pair_id}.dysgu.DUP_fc1.5.vcf")
	file("${pair_id}.delly.DUP_fc1.5.vcf")
	file("${pair_id}.smoove.DUP_fc1.5.vcf")
	file("*duphold.vcf")	
	tuple val(pair_id), file("${pair_id}.merged.DUP_survivor.vcf") into duplication_annot_calls_ch
	
	script:
	
	"""
	bcftools view -i '(SVTYPE = "DUP")' ${pair_id}.dysgu.vcf > ${pair_id}.dysgu.DUP.vcf
	
	export DUPHOLD_FLANK=2000
	export DUPHOLD_SAMPLE_NAME=${pair_id}
	
	duphold -t 4 -v ${pair_id}.dysgu.DUP.vcf -b ${pair_id}.sort.bam -f ${genome_file} -o ${pair_id}.dysgu.DUP_duphold.vcf
	
	bcftools view -i "((FMT/DHFC[0]>=1.3 & FMT/DHFFC[0]>=1.3 & FMT/DHBFC[0]>=1.3) && (INFO/PE>=1 && INFO/SR>=1))" ${pair_id}.dysgu.DUP_duphold.vcf > ${pair_id}.dysgu.DUP_fc.vcf
	
	bcftools view -i '(SVTYPE = "DUP")' ${pair_id}.delly.vcf > ${pair_id}.delly.DUP.vcf
	
	export DUPHOLD_FLANK=2000
	export DUPHOLD_SAMPLE_NAME=${pair_id}
	
	duphold -t 4 -v ${pair_id}.delly.DUP.vcf -b ${pair_id}.sort.bam -f ${genome_file} -o ${pair_id}.delly.DUP_duphold.vcf
	bcftools view -i "(FMT/DHFC[0]>=1.3 & FMT/DHFFC[0]>=1.3 & FMT/DHBFC[0]>=1.3)" ${pair_id}.delly.DUP_duphold.vcf > ${pair_id}.delly.DUP_fc.vcf
	
	bcftools view -i '(SVTYPE = "DUP")' ${pair_id}.smoove.vcf > ${pair_id}.smoove.DUP.vcf
	
	export DUPHOLD_FLANK=2000
	export DUPHOLD_SAMPLE_NAME=${pair_id}
	
	duphold -t 4 -v ${pair_id}.smoove.DUP.vcf -b ${pair_id}.sort.bam -f ${genome_file} -o ${pair_id}.smoove.DUP_duphold.vcf
	
	bcftools view -i "((FMT/DHFC[0]>=1.3 & FMT/DHFFC[0]>=1.3 & FMT/DHBFC[0]>=1.3) && (INFO/PE>=1 && INFO/SR>=1))" ${pair_id}.smoove.DUP_duphold.vcf > ${pair_id}.smoove.DUP_fc1.vcf
	
	ls ${pair_id}.dysgu.DUP_fc.vcf ${pair_id}.delly.DUP_fc.vcf ${pair_id}.smoove.DUP_fc.vcf > ${pair_id}.txt
	
	SURVIVOR merge ${pair_id}.txt 50 2 1 1 0 50 ${pair_id}.merged.DUP_survivor.vcf
	
	"""
}		

/*
 * Duplicated gene detection using the annotation file and vcf file merge to do intersection 
 */

process duplicate_Gene {

	publishDir "${params.out}/duplicated_gene", mode:'move'
	label "bcftools"
	label "bedtools"
	tag "Detection of duplicated genes"
	
	input:
	tuple val(pair_id), file(gene) from duplication_annot_calls_ch
	file(annotation)
	
	output:
	file("${pair_id}.merged.DUP_survivor.bed")
	file("${pair_id}.intersectBed.csv")
	tuple val(pair_id), file("${pair_id}.gene_duplicated.bed") into duplicated_gene_ch
	
	script:	
	
	"""
	bcftools query -f "%CHROM\\t%POS\\t%INFO/END\\t%ALT\\t%INFO/SVLEN\\n" ${pair_id}.merged.DUP_survivor.vcf > ${pair_id}.bed
	
	awk -F "\\t" '{if(\$5>500) {print \$0}}' ${pair_id}.bed > ${pair_id}.merged.DUP_survivor.bed

	
	bedtools intersect -a ${annotation} -b ${pair_id}.merged.DUP_survivor.bed -f 0.1 -wa| bedtools sort | uniq > ${pair_id}.gene_duplicated.bed

	
	bedtools intersect -a ${annotation} -b ${pair_id}.merged.DUP_survivor.bed -f 0.1 -wa -wb|uniq > ${pair_id}.intersectBed.csv
	
	"""
}

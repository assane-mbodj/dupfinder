/*
 * Détection duplicated gene using the annotaion file
 */

process GENE_DUPLICATED_LR {

	tag "Detection of duplicated genes"

	publishDir "${params.out}/duplicated_gene", mode: 'move'
	
	input:
	tuple val (prefix),  path(gene)
	file annotation
	
	output:
	path("${prefix}.merged_DUP_survivor_lr.bed")
	path("${prefix}.intersectBed_lr.csv")
	tuple val(prefix), path("${prefix}.gene_dup_lr.gff"), emit: duplicated_gene_ch
	
	script:	
	
	"""
	bcftools query -f "%CHROM\\t%POS\\t%INFO/END\\t%ALT\\t%INFO/SVLEN\\n" ${prefix}.merged_DUP_survivor_lr.vcf > ${prefix}_lr.bed
	
	awk -F "\\t" '{if(\$5>500) {print \$0}}' ${prefix}_lr.bed > ${prefix}.merged_DUP_survivor_lr.bed

	
	bedtools intersect -a ${annotation} -b ${prefix}.merged_DUP_survivor_lr.bed -f 0.40 -wa| bedtools sort | uniq > ${prefix}.gene_dup_lr.gff

	
	bedtools intersect -a ${annotation} -b ${prefix}.merged_DUP_survivor_lr.bed -f 0.40 -wa -wb|uniq > ${prefix}.intersectBed_lr.csv
	
	
	
	"""
}

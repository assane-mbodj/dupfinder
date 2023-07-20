/*
 * Détection duplicated gene using the annotaion file
 */

process GENE_DUPLICATED {

	tag "Detection of duplicated genes"

	publishDir "${params.out}/duplicated_gene", mode: 'move'
	
	input:
	tuple val (prefix),  path(gene)
	file annotation
	
	output:
	path("${prefix}.merged.DUP_survivor.bed")
	path("${prefix}.intersectBed.csv")
	tuple val(prefix), path("${prefix}.gene_duplicated.gff"), emit: duplicated_gene_ch
	
	script:	
	
	"""
	bcftools query -f "%CHROM\\t%POS\\t%INFO/END\\t%ALT\\t%INFO/SVLEN\\n" ${prefix}.merged.DUP_survivor.vcf > ${prefix}.bed
	
	awk -F "\\t" '{if(\$5>500) {print \$0}}' ${prefix}.bed > ${prefix}.merged.DUP_survivor.bed

	
	bedtools intersect -a ${annotation} -b ${prefix}.merged.DUP_survivor.bed -f 0.50 -wa| bedtools sort | uniq > ${prefix}.gene_duplicated.gff

	
	bedtools intersect -a ${annotation} -b ${prefix}.merged.DUP_survivor.bed -f 0.50 -wa -wb|uniq > ${prefix}.intersectBed.csv
	
	"""
}

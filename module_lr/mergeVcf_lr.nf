/*
 * Merging all vcf and keeping the duplication founding at least two callers
 */

process MERGE_LR {
	publishDir "${params.out}/merge_vcf", mode: 'copy'
	
	tag "Merging vcf files"
	
	input:
	tuple val(prefix), path(vcf)
	
	
	output:
		tuple val(prefix), path("${prefix}.merged_DUP_survivor_lr.vcf"), emit: filterVcf_file
		
	script:
	
	""" 
	ls ${prefix}.sniffles_DUP_fc.vcf ${prefix}.svim_DUP_fc.vcf ${prefix}.cutesv_DUP_fc.vcf > ${prefix}.txt
	
	SURVIVOR merge ${prefix}.txt 1000 2 1 1 0 50 ${prefix}.merged_DUP_survivor_lr.vcf
	
	"""
	}

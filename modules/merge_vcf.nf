/*
 * Merging all vcf and keeping the duplication founding at least two callers
 */

process VCF_MERGE {
	publishDir "${params.out}/merge_vcf", mode: 'copy'
	
	tag "Merging vcf files"
	
	input:
	tuple val(prefix), path(vcf)
	
	
	output:
		tuple val(prefix), path("${prefix}.merged.DUP_survivor.vcf"), emit: filter_file
		
	script:
	
	""" 
	ls ${prefix}.dysgu.DUP_fc.vcf ${prefix}.delly.DUP_fc.vcf ${prefix}.smoove.DUP_fc.vcf > ${prefix}.txt
	
	SURVIVOR merge ${prefix}.txt 1000 2 1 1 0 50 ${prefix}.merged.DUP_survivor.vcf
	
	"""
	}

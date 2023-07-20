/*
 * Variant caliing step using sniffles and svim tools
 */
process VARIANT_LR {	
	publishDir "${params.out}/variant_calls", mode: 'copy'

	
	tag " Variants calling ${prefix} "
	
	input: 
	tuple val(prefix), path(bam)
	file genome  

	
	output:
	tuple val(prefix), path("${prefix}*.vcf"),		emit: variant_calls

	
	
	script:
	
	"""	
	sniffles --minsvlen 50 --mapq 10 --sample-id ${prefix} --input ${prefix}.sort.bam --reference ${genome} --vcf ${prefix}.sniffles.vcf
	
	svim alignment --min_sv_size 50 --min_mapq 10 --sample ${prefix} SVIM ${prefix}.sort.bam ${genome}

	mv SVIM/variants.vcf .
	mv variants.vcf ${prefix}.svim.vcf
	
	cuteSV --threads 8 ${prefix}.sort.bam ${genome} ${prefix}.cutesv.vcf ${params.out} --sample ${prefix} --min_mapq 10

	"""
	
}


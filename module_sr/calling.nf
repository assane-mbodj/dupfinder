/*
 * Variant caliing step using dysgu, delly and smoove tools
 */
process VARIANT {	
	publishDir "${params.out}/variant_calls", mode: 'copy'

	
	tag " Variants calling "
	
	input: 
	tuple val(prefix), path(bam)
	file genome  
	file genome_index_file
	
	output:
	tuple val(prefix), path("${prefix}*.vcf"),		emit: variant_calls

	
	
	script:
	
	"""	
	delly call --map-qual 10 --min-clique-size 1 -g ${genome} ${prefix}.sort.bam -o ${prefix}.delly.bcf
	
	bcftools view ${prefix}.delly.bcf > ${prefix}.delly.vcf
	
	dysgu run --clean --mq 10 --min-size 50 -p4 ${genome} temp_dir ${prefix}.sort.bam > ${prefix}.dysgu.vcf
	
	smoove call --name ${prefix} --outdir smoove/ --fasta ${genome} -p 4 ${prefix}.sort.bam 
	
	mv smoove/${prefix}-smoove.vcf.gz .
	mv ${prefix}-smoove.vcf.gz ${prefix}.smoove.vcf.gz
	bgzip -d ${prefix}.smoove.vcf.gz
	"""
	
}


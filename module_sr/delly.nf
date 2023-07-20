process DELLY {
	publishDir "${params.out}/variant_calls", mode:'copy'
	publishDir "${params.out}/duplication_annot_calls", mode:'copy'
	
	tag { prefix }
	
	input: 
	tuple val(prefix), path(bam)
	file genome  
	output:
	tuple val(prefix), path("${prefix}.delly.vcf"),		emit: variant_calls
	tuple val(prefix), path("${prefix}.delly.DUP_fc.vcf"),        emit: duplication_annot_calls
	
	script:
	
	"""
	delly call --map-qual 10 --min-clique-size 1 -g ${genome} ${prefix}.sort.bam -o ${prefix}.delly.bcf
	
	bcftools view ${prefix}.delly.bcf > ${prefix}.delly.vcf
	
	
	bcftools view -i '(SVTYPE = "DUP")' ${prefix}.delly.vcf > ${prefix}.delly.DUP.vcf
	
	export DUPHOLD_FLANK=2000
	export DUPHOLD_SAMPLE_NAME=${prefix}.sort
	
	duphold -t 4 -v ${prefix}.delly.DUP.vcf -b ${prefix}.sort.bam -f ${genome} -o ${prefix}.delly.DUP_duphold.vcf
	bcftools view -i "((FMT/DHFFC[0]>=1.3 & FMT/DHBFC[0]>=1.3) && (INFO/PE>=1))" ${prefix}.delly.DUP_duphold.vcf > ${prefix}.delly.DUP_fc.vcf
	
	"""
	
	}

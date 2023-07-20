/*
 * Filtering and Annotaion steps for the Duplications
 */
process VCF_FILTER {
	publishDir "${params.out}/duplication_annot_calls", mode:'copy'

	tag "DUP filtering and annotation"
	
	input:
	tuple val(prefix),  path(bam)
	tuple val(prefix), path(vcf)
	file genome 
	file genome_index_file
	

	
	output:
	
	tuple val(prefix), path("${prefix}.*.DUP_fc.vcf"),        emit: duplication_annot_calls
	
	
	script:
	
	"""	
	
	# Dysgu step
	
	bcftools view -i '(SVTYPE = "DUP")' ${prefix}.dysgu.vcf > ${prefix}.dysgu.DUP.vcf
	
	export DUPHOLD_FLANK=2000
	export DUPHOLD_SAMPLE_NAME=${prefix}.sort.dysgu_reads
	
	duphold -t 4 -v ${prefix}.dysgu.DUP.vcf -b ${prefix}.sort.bam -f ${genome} -o ${prefix}.dysgu.DUP_duphold.vcf
	
	bcftools view -i "((FMT/DHFFC[0]>=1.3 & FMT/DHBFC[0]>=1.3) && (INFO/PE>=1 & INFO/SR>=1))" ${prefix}.dysgu.DUP_duphold.vcf > ${prefix}.dysgu.DUP_fc.vcf
	
	
	# Delly step
	
	bcftools view -i '(SVTYPE = "DUP")' ${prefix}.delly.vcf > ${prefix}.delly.DUP.vcf
	
	export DUPHOLD_FLANK=2000
	export DUPHOLD_SAMPLE_NAME=${prefix}.sort
	
	duphold -t 4 -v ${prefix}.delly.DUP.vcf -b ${prefix}.sort.bam -f ${genome} -o ${prefix}.delly.DUP_duphold.vcf
	bcftools view -i "((FMT/DHFFC[0]>=1.3 & FMT/DHBFC[0]>=1.3) && (INFO/PE>=1 & INFO/SR>=1))" ${prefix}.delly.DUP_duphold.vcf > ${prefix}.delly.DUP_fc.vcf
	
	
	#Smoove step
	
	bcftools view -i '(SVTYPE = "DUP")' ${prefix}.smoove.vcf > ${prefix}.smoove.DUP.vcf
	
	export DUPHOLD_FLANK=2000
	export DUPHOLD_SAMPLE_NAME=${prefix}-sort
	
	duphold -t 4 -v ${prefix}.smoove.DUP.vcf -b ${prefix}.sort.bam -f ${genome} -o ${prefix}.smoove.DUP_duphold.vcf
	
	bcftools view -i "((FMT/DHFFC[0]>=1.3 & FMT/DHBFC[0]>=1.3) && (INFO/PE>=1 & INFO/SR>=1))" ${prefix}.smoove.DUP_duphold.vcf > ${prefix}.smoove.DUP_fc.vcf
	"""
}

/*
 * Filtering and Annotaion steps for the Duplications
 */
process FILTER_LR {
	publishDir "${params.out}/duplication_annot_calls", mode:'copy'

	tag "DUP filtering and annotation"
	
	input:
	tuple val(prefix),  path(bam)
	tuple val(prefix), path(vcf)
	file genome 
	

	
	output:
	
	tuple val(prefix), path("${prefix}.*DUP_fc.vcf"),        emit: duplication_annot_calls
	
	
	script:
	
	"""	
	#Sniffles step

	
	bcftools view -i '(SVTYPE = "DUP")' ${prefix}.sniffles.vcf > ${prefix}.sniffles.DUP.vcf
	
	export DUPHOLD_FLANK=2000
	export DUPHOLD_SAMPLE_NAME=${prefix}
	
	duphold -t 8 -v ${prefix}.sniffles.DUP.vcf -b ${prefix}.sort.bam -f ${genome} -o ${prefix}.sniffles_DUP_duphold.vcf
	
	bcftools view -i "(FMT/DHFFC[0]>=1.3 & FMT/DHBFC[0]>=1.3)" ${prefix}.sniffles_DUP_duphold.vcf > ${prefix}.sniffles_DUP_fc.vcf
	
	
	#Svim step
	
	bcftools view -i '(SVTYPE = "DUP:TANDEM" || SVTYPE = "DUP:INT")' ${prefix}.svim.vcf > ${prefix}.svim_DUP.vcf
	
	export DUPHOLD_FLANK=2000
	export DUPHOLD_SAMPLE_NAME=${prefix}
	
	duphold -t 8 -v ${prefix}.svim_DUP.vcf -b ${prefix}.sort.bam -f ${genome} -o ${prefix}.svim_DUP_duphold.vcf
	
	bcftools view -i "(FMT/DHFFC[0]>=1.3 & FMT/DHBFC[0]>=1.3)" ${prefix}.svim_DUP_duphold.vcf > ${prefix}.svim_DUP_fc.vcf
	
	
	# cutesv2 step
	
	bcftools view -i '(SVTYPE = "DUP")' ${prefix}.cutesv.vcf > ${prefix}.cutesv_DUP.vcf
	
	export DUPHOLD_FLANK=2000
	export DUPHOLD_SAMPLE_NAME=${prefix}
	
	duphold -t 8 -v ${prefix}.cutesv_DUP.vcf -b ${prefix}.sort.bam -f ${genome} -o ${prefix}.cutesv_DUP_duphold.vcf
	
	bcftools view -i "(FMT/DHFFC[0]>=1.3 & FMT/DHBFC[0]>=1.3)" ${prefix}.cutesv_DUP_duphold.vcf > ${prefix}.cutesv_DUP_fc.vcf
	
	"""
	
}

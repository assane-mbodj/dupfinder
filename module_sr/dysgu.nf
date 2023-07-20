process DYSGU {

	publishDir "${params.out}/variant_calls", mode:'copy'
	publishDir "${params.out}/duplication_annot_calls", mode:'copy'
	
	tag { prefix }
	
	input: 
	tuple val(prefix), path(bam)
	file genome  
	
	
	output:
	tuple val(prefix), path("${prefix}.dysgu.vcf"),		emit: variant_calls
	tuple val(prefix), path("${prefix}.dysgu.DUP_fc.vcf"), 	emit: duplication_annot_calls
	
	script:
	
	"""
	dysgu run --clean --mq 10 --min-size 50 -p4 ${genome} temp_dir ${prefix}.sort.bam > ${prefix}.dysgu.vcf
	
	bcftools view -i '(SVTYPE = "DUP")' ${prefix}.dysgu.vcf > ${prefix}.dysgu.DUP.vcf
	
	export DUPHOLD_FLANK=2000
	export DUPHOLD_SAMPLE_NAME=${prefix}.sort.dysgu_reads
	
	duphold -t 4 -v ${prefix}.dysgu.DUP.vcf -b ${prefix}.sort.bam -f ${genome} -o ${prefix}.dysgu.DUP_duphold.vcf
	
	bcftools view -i "((FMT/DHFFC[0]>=1.3 & FMT/DHBFC[0]>=1.3) && (INFO/PE>=1 || INFO/SR>=1))" ${prefix}.dysgu.DUP_duphold.vcf > ${prefix}.dysgu.DUP_fc.vcf
	
	"""
	
	}

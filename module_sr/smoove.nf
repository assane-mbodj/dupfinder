process SMOOVE {
	publishDir "${params.out}/variant_calls", mode:'copy'
	publishDir "${params.out}/duplication_annot_calls", mode:'copy'
	
	tag { prefix }
	
	input: 
	tuple val(prefix), path(bam)
	file genome
	file genome_index_file  
	output:
	tuple val(prefix), path("${prefix}.smoove.vcf"),		emit: variant_calls
	tuple val(prefix), path("${prefix}.smoove.DUP_fc.vcf"),	emit: duplication_annot_calls
	
	script:
	
	"""
	smoove call --name ${prefix} --outdir smoove/ --fasta ${genome} -p 4 ${prefix}.sort.bam 
	
	mv smoove/${prefix}-smoove.vcf.gz .
	mv ${prefix}-smoove.vcf.gz ${prefix}.smoove.vcf.gz
	bgzip -d ${prefix}.smoove.vcf.gz
	
	
	bcftools view -i '(SVTYPE = "DUP")' ${prefix}.smoove.vcf > ${prefix}.smoove.DUP.vcf
	
	export DUPHOLD_FLANK=2000
	export DUPHOLD_SAMPLE_NAME=${prefix}-sort
	
	duphold -t 4 -v ${prefix}.smoove.DUP.vcf -b ${prefix}.sort.bam -f ${genome} -o ${prefix}.smoove.DUP_duphold.vcf
	
	bcftools view -i "((FMT/DHFFC[0]>=1.3 & FMT/DHBFC[0]>=1.3) && (INFO/PE>=1 || INFO/SR>=1))" ${prefix}.smoove.DUP_duphold.vcf > ${prefix}.smoove.DUP_fc.vcf
	
	"""
	
	}
	
	
process DELLY {
	publishDir "${params.out}/variant_calls", mode:'copy'
	publishDir "${params.out}/duplication_annot_calls", mode:'copy'
	
	tag { prefix }
	
	input: 
	tuple val(prefix), path(bam)
	file genome  
	file genome_index_file
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
	
	
process DYSGU {

	publishDir "${params.out}/variant_calls", mode:'copy'
	publishDir "${params.out}/duplication_annot_calls", mode:'copy'
	
	tag { prefix }
	
	input: 
	tuple val(prefix), path(bam)
	file genome  
	file genome_index_file
	
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
	
	
process VCF_MERGE {
	publishDir "${params.out}/merge_vcf", mode: 'copy'
	
	tag {prefix}
	
	input:
	tuple val(prefix), path(vcf)
	
	
	output:
		tuple val(prefix), path(vcf), path("${prefix}.merged.DUP_survivor.vcf"), emit: filter_file
		
	script:
	
	""" 
	ls ${prefix}.dysgu.DUP_fc.vcf ${prefix}.delly.DUP_fc.vcf ${prefix}.smoove.DUP_fc.vcf > ${prefix}.txt
	
	SURVIVOR merge ${prefix}.txt 1000 2 1 1 0 50 ${prefix}.merged.DUP_survivor.vcf
	
	"""
	}

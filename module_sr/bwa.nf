/*
 * Mapping step using bwa
 */

process MAPPING {
	
	tag "Aligning reads with bwa mem:"
       publishDir "${params.out}/Bam_file", mode : "copy"
       input:
       tuple val(prefix), path(reads)
       file genome
       file genome_index_file
       file genome_bwa_amb_file
	file genome_bwa_ann_file
	file genome_bwa_bwt_file
	file genome_bwa_pac_file
	file genome_bwa_sa_file

       output:
       tuple val(prefix), path("*.bam*"), emit: bam_file_ch
       
       
       script:
       
       """
          bwa mem -t 8 -R "@RG\\tID:${prefix}\\tSM:${prefix}\\tLB:${prefix}" ${genome} ${reads[0]} ${reads[1]} | samtools view -u -bS | samtools sort > ${prefix}.sort.bam
          samtools index ${prefix}.sort.bam
       """
}

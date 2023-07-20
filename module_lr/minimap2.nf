/*
 * Mapping step using minimap2
 */

process MAPPING_LR {
	
	tag "Aligning reads using minimap2:"
       publishDir "${params.out}/Bam_file", mode : "copy"
       input:
       tuple val(prefix), path(reads)
       path genome
       

       output:
       tuple val(prefix), path("*.bam*"), emit: bam_file_ch
       
       
       script:
       
       """
       minimap2 --MD -t 8 -ax map-ont ${genome} ${reads} | samtools view -u -bS | samtools sort > ${prefix}.sort.bam
            
          samtools index ${prefix}.sort.bam
          
       """
}

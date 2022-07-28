#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
This is a bioniformatics pipeline written to practice nextflow


workflow:

1. retrieve bam file and split into chromosomes
2. perform haplotype calling using GATK > vcf_filtering
3. Gather vcfs from all chromosome and merge into one file

*/


process SPLIT_BAM_TO_CHROMOSOMES {
  publishDir "${params.outdir}/bams"
  input:
  	val bam
  output:
   	stdout emit: verbiage
  script:
  """

  bamtools split -in $bam -reference
  """
}

process HAPLOTYPE_CALLER {
   publishDir "${params.outdir}"

   input:
     	val bam
    output:
    	stdout emit: verbiage
    	path "*"

  script:
  	"""
  samtools index ${bam}
  ${params.gatk} HaplotypeCaller --input ${bam} --reference ${params.reference} --output ${bam.baseName}.vcf > ${bam.baseName}.txt
  vcftools --vcf ${bam.baseName}.full.vcf --minGQ 15 --minDP 20 > "${bam.baseName}.vcf"
  	"""
}


process GATHER_VCF {
    publishDir "${params.outdir}"

	input:
		val vcf_path
	output:
    	stdout emit: verbiage
    	path "*"

	script:
    """
    echo ${vcf_path}
    # grep -h '^#' ${vcf_path}/*REF_22.vcf > "merge.vcf"
    # grep -hv '^#' ${vcf_path}/*.vcf >> "merge.vcf"

    cat ${vcf_path}/*.txt > merge.txt
    """
}


def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run vcf_calling_by_chromosome.ln --bam "bam" --reference_genome "reference_genome"

        Mandatory arguments:
         --bam                          BAM to run vcf calling by chromosome on
         --reference_genome				Reference genome you'd like to use

       Optional arguments:
        --outdir                       Output directory to place produced data
        --threads					   Number of threads
        --help                         This usage statement.
        """
}

workflow {
// 	Show help message
	if (params.help) {
    	helpMessage()
    	exit 0
	}

	println "$params"
	parentdir=file("$params.bam").Parent
	SPLIT_BAM_TO_CHROMOSOMES("$params.bam")
	contig = Channel.fromPath("$parentdir/*REF*.bam")
	HAPLOTYPE_CALLER(contig)
	GATHER_VCF("$workflow.projectDir/$params.outdir")

}
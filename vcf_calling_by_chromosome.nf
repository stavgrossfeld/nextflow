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
  publishDir "${params.outdir}/bams", mode: "copy"
  input:
  	path bam
  output:
   	path "*.bam"

  script:
  """
  bamtools split -in ${bam} -reference
  """
}

process HAPLOTYPE_CALLER {
   debug true
   publishDir "${params.outdir}"

   input:
   	    val(contig_file)
    output:
    	path "filtered.vcf"

  script:

	"""

	# index bam file

    samtools index $contig_file

	# haplotype calling
    ${params.gatk} HaplotypeCaller --input $contig_file --reference ${params.reference} --output "full.vcf"

	# select SNP variants
	${params.gatk} SelectVariants \
        -R ${params.reference} \
        -V "full.vcf" \
        -select-type SNP \
        -O "snp.vcf"

	# Filter on criteria
    ${params.gatk} VariantFiltration \
		-R ${params.reference} \
        -V "snp.vcf"  \
        -O "filtered.vcf"  \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"


	"""
}


process GATHER_VCF {
   debug true
   publishDir "${params.outdir}"

   input:
   		val(vcfs)
   		val(header)
//    	    path '*snp.vcf'
    output:
    	path "*.vcf"

	script:
	"""
   	grep '##' $header > results.txt
   	cat $vcfs | grep -v '##' >> results.vcf
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

	// split bam into contigs
	split_bam_ch = SPLIT_BAM_TO_CHROMOSOMES("$params.bam").collect().flatten().take(3)
	// run haplotype caller
	split_vcf_ch = HAPLOTYPE_CALLER(split_bam_ch)
				.collect()

	// gather vcf using one of the files for the header
	header_file = split_vcf_ch.flatten().first()
	header_file.view()
	split_vcf_str = split_vcf_ch.map { it -> it.join(" ") }
	GATHER_VCF(split_vcf_str, header_file)







}
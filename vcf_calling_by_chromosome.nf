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
   	stdout emit: verbiage

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
   	    val(parentdir)
    output:
    	path "*.vcf"

  script:

	"""

	# index bam file

	samtools index "${baseDir}/results/bams/${contig_file}.bam"

	# haplotype calling
	${params.gatk} HaplotypeCaller --input "${baseDir}/results/bams/${contig_file}.bam" --reference ${params.reference} --output "full.vcf"

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
// 	SPLIT_BAM_TO_CHROMOSOMES("$params.bam")

	knownChromosomes = (1..23).collect{ it.toString() } + ["X", "Y"]

	split_bam_ch = Channel.fromPath("results/bams/*REF_*.bam", checkIfExists: true)
		.map {it.baseName}


// 		| map {it.name.lastIndexOf('.')}
// 		| view()
// 	    | filter { it.name.split("_") .last() in knownChromosomes }
//     	| view()



// 		.out
//
// 	split_bam_ch = Channel.fromPath("$parentdir/*REF_*.bam")
// 						  .map { parentdir.toString()  + '/' + it.baseName }.flatten().take(1)

//
// // 	split_bam_ch.view()
	HAPLOTYPE_CALLER(split_bam_ch, parentdir)
		.collectFile(name:'result.txt', storeDir: params.outdir)
// 		.view { it.text }

// 	HAPLOTYPE_CALLER(split_bam_ch, parentdir)
// 		.collectFile(name:'result.txt', storeDir: params.outdir)
// 		.view { it.text }

}
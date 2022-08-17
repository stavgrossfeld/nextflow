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
   		tuple val(contig_id), path(bam)
    output:
//     	stdout emit: verbiage
    	path "*.vcf"

  script:
  	"""

  	echo "$bam.baseName".REF_"$contig_id".bam > "$bam.baseName".REF_"$contig_id".vcf
  # echo samtools index ${bam}
  # echo ${params.gatk} HaplotypeCaller --input ${bam} --reference ${params.reference} --output ${bam.baseName}.vcf > ${bam.baseName}.txt
  # ec ho vcftools --vcf ${bam.baseName}.full.vcf --minGQ 15 --minDP 20 > "${bam.baseName}.vcf"
  	"""
}


process GATHER_VCF {
    publishDir "${params.outdir}"

	input:
		val vcf_path
	output:

    	path "*"

	script:
    """
    echo ${vcf_path}
    # grep -h '^#' ${vcf_path}/*REF_22.vcf > "merge.vcf"
    # grep -hv '^#' ${vcf_path}/*.vcf >> "merge.vcf"

    # cat ${vcf_path}/*.txt > merge.txt
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


process blast {
  input:
  	 tuple val(x), path(y)

//     path "bam.simpleName.REF_{contig_id}.bam"


//   output:
//     tuple val(x), path('result') into blastOuts

  script:
    """
    echo $x
    echo $y
    echo $x[0]
    echo $x[1]

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

	bam_ch = Channel.fromPath("$params.bam")

// 	SPLIT_BAM_TO_CHROMOSOMES("$params.bam")
	contig_ch = Channel.of(1..23, 'X', 'Y')

	contig_ch = contig_ch.combine(bam_ch)


	gather_ch = HAPLOTYPE_CALLER(contig_ch)


	gather_ch
		.collect()
		.collectFile(name:'result.txt', storeDir: params.outdir)
// 		.view { it.text }

// 	GATHER_VCF([gather_ch])

}
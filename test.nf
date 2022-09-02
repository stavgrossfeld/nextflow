

process BAMTOOLS_SPLIT {
    publishDir "${params.outdir}/bams", mode: "copy"

    container "quay.io/biocontainers/bamtools:2.5.2--hd03093a_0"

    input:
    path bam

    output:
    path "*.bam"
//     emit: bam

    script:
    """
    bamtools split -in $bam -reference
    """
}



workflow {
	BAMTOOLS_SPLIT($params.bam)
}
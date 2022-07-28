# Running individually on one bam
``` 
Launching `vcf_calling_by_chromosome.nf` [loving_hugle] DSL2 - revision: 94528186c3

        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run vcf_calling_by_chromosome.ln --bam "bam" --reference_genome "reference_genome"

        Mandatory arguments:
         --bam                          BAM to run vcf calling by chromosome on
         --reference_genome				Reference genome you'd like to use
         --dbName                       Prefix name of the BLAST database

       Optional arguments:
        --outdir                       Output directory to place final BLAST output
        --threads					   Number of threads
        --help                         This usage statement.

```


# Running on multiple bams in production

```
reference_genome="reference_genome" s3=* plate_id=* project_dir=~/Desktop/twinstrand_project/genome_in_a_bottle/151002_ashkenazi_trio ./run_nextflow.sh
```
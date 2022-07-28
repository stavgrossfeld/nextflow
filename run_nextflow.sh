# usr/bin/env bash
set -euxo pipefail
DIR=`realpath $(dirname "$0")`

echo $DIR
# retreive plate address and download to directory
echo """s3_address = python retreive_s3_address_for_plate ${plate_id}"""
echo mkdir -p ${project_dir}
echo aws s3 sync ${s3_address} ${project_dir}/.

# save sample list
ls ${project_dir}/*.bam > sample_list.txt

echo "\n\n\n"

# for each bam in directory run vcf_calling_by_chromosome

run_nextflow() {
  basename_bam=`echo $(basename "$1") | tr '.' '_' `

  echo $basename_bam
  # create directories
  mkdir -p ${basename_bam}
  mkdir -p ${basename_bam}/bam

  cd ${basename_bam}

  echo "copy original bam to nextflow dir of bam"
  cp ${DIR}/vcf_calling_by_chromosome.nf ${DIR}/nextflow.config .
  og_bam_path=`realpath ./bam/og_bam.bam`
  echo $og_bam_path

  echo "run vcf_calling_by_chromosome.nf on copied bam \n\n\n\n"
  nextflow run vcf_calling_by_chromosome.nf --bam $og_bam_path --reference /mnt/efs/${reference_genome}
}


# Collect sample names
samples=()
for file in `cat sample_list.txt`
do
    samples+=("${file}")

done

echo "${#samples[@]} Samples found: ${samples[*]} \n\n\n\n"

for sample in ${samples[@]}
do
  run_nextflow $sample
done


######### TO DO: make parallel execution on multiple samples at once, need large instance
#export -f run_nextflow
#parallel ::: run_nextflow ${samples[@]}

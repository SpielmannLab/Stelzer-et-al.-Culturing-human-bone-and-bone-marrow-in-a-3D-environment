#! /bin/bash
# Run this script make fastq files from illumina bcl files
# Use this instead of mkfastq_main.sh, when using cell multiplexing oligo, since it requires cellranger version > 6.1.0
# STEP1. Copy the Run folder into the Cluster
# STEP2. Provide the details in the INPUT parameters section
# STEP3. Close this file and execute it using ./mkfastq_cellplex_main.sh

# INPUT Parameters
# This is the Illumina-generated run-directory /*/*/*
path_illuminaRun="/data/humangen_singlecell/Stelzer/Illumina/240430_VH00398_77_2222LVFNX/"
# This is the full path of the simplified 10x-sample sheet /*/*/*.csv
path_sampleSheet="/data/humangen_singlecell/Stelzer/Illumina/240430_VH00398_77_2222LVFNX/SampleSheet.csv"
# This is the name of the project/sample/run
id="P045_Stelzer"
# This is the place where the output will be stored. This should be an existing folder
path_out="/data/humangen_singlecell/Stelzer/Illumina/240430_VH00398_77_2222LVFNX/fastq"
# How many mismatches are allowed in the index reads. Base this on Hamming Distance.
# Default is 1. Possible are 0,1 or 2. Can also provide comma-separated values for each of the samples, but read bcl2fastq2 documentation if this is needed.
hamming_dist=1

# DO NOT CHANGE THINGS BELOW
if [[ -d ${path_out}/${id} ]]; then
	echo "Error: The directory ${path_out}/${id} already exists"
	exit 1
else
	sbatch -c 32 --mem 256GB --job-name=mkfastq_${id} mkfastq_cellplex_core.sh $path_illuminaRun $path_sampleSheet $id $path_out $hamming_dist
	echo "Message: Job sent to mkfastq pipeline"
	echo "Look for the slurm-mkfastq_${id}-xxxxxx.out to follow the progress"
fi

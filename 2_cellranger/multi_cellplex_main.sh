#! /bin/bash
# Run this script to get cell x gene matrix from fastq files, corresponding to cell multiplexing dataset.
# STEP1. Make sure you have the fastq files from the cell multiplexing as well as GEX libraries. If not, run the  cellplex_mkfastq_main.sh script
# STEP2. Fill out the two CSV sheets. One for the parameters and the other one regarding the cell multiplexing oligo sequences.
# STEP3. Provide all the details in the INPUT parameters section
# STEP4. Close this file and execute it using ./multi_cellplex_main.sh


# INPUT Parameters
# This is full address the main csv file
path_to_csv="/data/humangen_singlecell/Stelzer/fastq/config.csv"
# This is the place where output will be stored. The directory may or may not exist. Beware, contents may be over-written /*/*/*
path_results="/data/humangen_singlecell/Stelzer/count"
# Provide an id to
sample="P45" #A unique run id and output folder name a-zA-Z0-9_

# DO NOT CHANGE THINGS HERE
# if the provided path_out does not exist, create it
mkdir -p $path_results

sbatch -c 16 --job-name=cellplex_${sample} multi_cellplex_core.sh $sample $path_to_csv $path_results
echo "Message: Job cellplex_${sample} sent to cellranger multi pipline for cellplex counting"
echo "Look for slurm-cellplex_${sample}-xxxxxx.out to follow progress"

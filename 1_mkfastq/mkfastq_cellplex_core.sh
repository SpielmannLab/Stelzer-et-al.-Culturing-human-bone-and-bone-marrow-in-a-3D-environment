#! /bin/bash
### Submit this Script with: sbatch script.sh ###
# Uses cellranger v7.0.0 to allow for cellplex (using cell multiplexing oligo)

# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition shortterm:
#SBATCH --partition=shortterm
#  Use one node:
#SBATCH --nodes=1
#  Request 32 cores (hard constraint):
#SBATCH -c 32
#  Request 256GB of memory (hard constraint):
#SBATCH --mem=256GB
#  Set slurm file output nomenclature
#SBATCH --output "slurm-%x-%j.out"

path_run="$1"
path_sheet="$2"
id="$3"
path_out="$4"
hamming_dist="$5"

# Documenting the input parameters:
echo "--------------------------------------"
echo "Message: These are the run parameters:"
echo "Illumina Run Path: $path_run"
echo "10X-style simplified sample sheet: $path_sheet"
echo "id: $id"
echo "Output will be stored at: $path_out"
echo "--------------------------------------"

# Load your necessary modules (example):
module load cellranger/7.0.0
module load bcl2fastq/v2.20

#Prepare $SCRATCH
mkdir -p  $SCRATCH/run

# Move to the scratch directory, where output will be stored within a directory with the name "$id"
cd $SCRATCH

echo "Message: Copying data to scratch for fast processing"
# Copy input data to scratch for much faster execution:
cp -r $path_run/* $SCRATCH/run/
cp $path_sheet $SCRATCH/
echo "Message: Data copied to scratch for fast processing"

# Running the mkfastq pipeline
# Reference: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq
cellranger_7.0.0.simg cellranger mkfastq \
	--id="$id" \
 	--run="$SCRATCH/run/" \
 	--csv="$path_sheet" \
	--barcode-mismatches="$hamming_dist"

echo "Message: Cell Ranger execution completed"


echo "Message: Copying output to $path_out"

# Save the results from the scratch folder
cp -a $id $path_out/
cp *.mro $path_out/

echo "Message: Copied all output to $path_out"

# Clean up scratch
rm -rf $SCRATCH/*

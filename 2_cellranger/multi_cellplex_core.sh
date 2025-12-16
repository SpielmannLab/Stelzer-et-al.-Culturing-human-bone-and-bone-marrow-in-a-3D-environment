#! /bin/bash
### Submit this Script with: sbatch multi_cellplex_core.sh ###

# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=shortterm
#  Use one node:
#SBATCH --nodes=1
#  Request 32 cores (hard constraint):
#SBATCH -c 32
#  Request 256GB of memory (hard constraint):
#SBATCH --mem=256GB
# set slurm file output nomenclature:
#SBATCH --output "slurm-%x-%j.out"

sample="$1"
path_to_csv="$2"
path_results="$3"

# sample="test_hackel_simulation"
# path_to_csv="/data/humangen_singlecell/Project_2022_001_10x_Hackel_PBMC/Lawlor_stimulationData/analysis/Lawlor_HTO_RNA_config.csv"
# path_results="/data/humangen_mouse/test_area/varun/"

# Documenting the input parameters
echo "-------------------------------------"
echo "Message: These are the run parameters"
echo "Sample: ${sample}"
echo "Path containing the parameters: ${path_to_csv}"
echo "Output folder: ${path_results}"
echo "-------------------------------------"

# Load the necessary modules (Cellranger cell multiplexing requires cell ranger >v6):
module load cellranger/7.0.0

#Prepare $SCRATCH
cd "$SCRATCH"

# Run the cellranger pipeline
echo "Message: Running Cellranger multi pipeline"
cellranger_7.0.0.simg cellranger multi --id=${sample} --csv=${path_to_csv}
echo "Message: Cell Ranger execution completed"

# Save the results (only the outs) from the scratch folder and cleanup scratch
mkdir -pv ${path_results}/${sample}
cp -av $SCRATCH/${sample}/outs $path_results/${sample}/
echo "Message: Output saved at $path_results"

# Code to check the HTO when checking data
# bioawk -c fastx '{print substr($seq,2,16)}' /data/humangen_singlecell/Project_2022_001_10x_Hackel_PBMC/Lawlor_stimulationData/1-HTO_pooled_and_merged_S1_L001_R2_001.fastq.gz | head -1000 | sort | uniq

#! /bin/bash
### Submit this Script with: sbatch script.sh ###
 
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
path_refdata="$2"
path_fastqdata="$3"
path_results="$4"
intron="$5"

# Documenting the input parameters
echo "-------------------------------------"
echo "Message: These are the run parameters"
echo "Sample: ${sample}"
echo "Transcriptome reference: ${path_refdata}"
echo "Input fastq data: ${path_fastqdata}"
echo "Output folder: ${path_results}"
echo "Include introns?: ${intron}"
echo "-------------------------------------"

#Prepare $SCRATCH
mkdir -p  $SCRATCH/inputdata/{fastq,refdata} $SCRATCH/scripts/ $SCRATCH/results/

#parse if multiple fastq paths were provided.
IFS=',' read -r -a path_fastq_array <<< "$path_fastqdata"

# Copy input data into your local job directory for much faster execution:
# copy for each flowcell provided
for path_fastq in ${path_fastq_array[@]}
do
	flowid=$(basename $path_fastq)
	
	mkdir $SCRATCH/inputdata/fastq/$flowid
	cp $path_fastq/${sample}_*.gz ${SCRATCH}/inputdata/fastq/${flowid}/
done

cp -a ${path_refdata}/* $SCRATCH/inputdata/refdata/

echo "Message: Data copied to scratch for fast processing, as below"
tree $SCRATCH

# creating a comma separated string of fastq paths
path_scratch_fastq_array=($SCRATCH/inputdata/fastq/*)
data_string="${path_scratch_fastq_array[*]}"
paths="${data_string//${IFS:0:1}/,}"

# Set path to the results, where cell ranger will output the data
cd "$SCRATCH/results/"
echo "Message: The current directory has been set to $(pwd)"

# Load the necessary environment
module load singularity
export SINGULARITY_CACHEDIR="$WORK/singularity"
singularity pull docker://varunkas/cellranger:8.0.1

singularity exec docker://varunkas/cellranger:8.0.1 cellranger count --id=$sample \
		--fastqs=$paths \
		--transcriptome="$SCRATCH/inputdata/refdata" \
		--sample=$sample \
    --create-bam=true \
    --include-introns=${intron} \
    --nosecondary

echo "Message: Cell Ranger execution completed"

# Save the results from the scratch folder
mkdir -p ${path_results}/${sample}
rsync -r ${sample}/outs ${path_results}/${sample}/

echo "Message: Output saved at $path_results"

#! /bin/bash

### Submit this Script with: sbatch <script.sh> ###

# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition shortterm/debug/longterm:
#SBATCH --partition=shortterm
#  Use so many node:
#SBATCH --nodes=1
#  Request so many cores (hard constraint):
#SBATCH -c 4
#  Request so much of memory (hard constraint):
#SBATCH --mem=100GB
#  Find your job easier with a name:
#SBATCH -J "make_plots"
# set slurm file output nomenclature
#SBATCH --output "slurm-%x-%j.out"

PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default

# copy the scrits to $WORK directory - otherwise there is a jav aerror
mkdir -p $WORK/make_plots_nextflow_launchdir
rsync -r --update * $WORK/make_plots_nextflow_launchdir/
cd $WORK/make_plots_nextflow_launchdir
chmod +x bin/*
rm slurm*

module load nextflow/v22.04.1

# if -resume option not present, then clean the nextflow
if [[ $1 = -resume ]]; then
    echo "Resuming the previous nextflow run"
else
    echo "Cleaning up all the previous nextflow runs"
    echo "If you wanted to resume the previous run, use the '-resume' option"
    while [[ $(nextflow log | wc -l) -gt 1 ]]; do
        nextflow clean -f
    done
fi

# Submit the Nextflow Script:
nextflow run make_plots.nf -params-file make_plots_params.yaml -profile omics -resume

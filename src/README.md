# Read me

## Introduction
This repository contains several workflows for single cell analysis. They all need an *xxx.rds* file containing a Seurat object as an input, in addition to many parameters. The pipelines include:

- **sc_p_sample:** Run this script to carry out the basic Single-Cell-RNA analysis, starting from a raw Seurat object, qc, filtering, doublet detection, normalization, dimensionality reduction, clustering, and differential expression analysis for a single sample. Multiple samples can be run in parallel as instructed in the *sc_p_sample_params.yaml* file.  <br />
- **sc_multi_sample** Run this script to carry out the merging/integration of the multiple sample, followed by dimensionality reduction, clustering, and differential expression analysis.
- **sc_subcluster** Run this script to carry out the subclustering of the multiple clusters. You will end up with individual Seurat object (i.e., *.rds file) for each cluster that was subclustered. The subclusters can be inetgrated, if required aswell. It also performs dimensionality reduction, (sub)-clustering, as well as differential expression on each of these cluster-specific *.rds files. 
- **sc_velo** Run this script to carry out the RNA-velocity analysis on any Seurat object. Make sure the Seurat object contains "spliced" and "unspliced" count matrices (asssays).
- **sc_de_genes_group_sample** Run this script for differntial expression analysis between conditions/treatments/replicates/batches/or any other grouping parameters. The outputs are volcano plot showing the log fold changes and violin plot of the top 2 genes in each group per cluster.
- **sc_trajectory** Run this script for trajectory and pseudotime analysis of an already analysed Seurat object. It includes the conversion of the dataset into a monocle object, followed by trajectory analysis. You can decide to import the UMAP embeddings from the Seurat object, or let monocle create a new embedding.

## Installation
The *installation* folder contains all the files needed to set up conda environments to run all the above pipelines. These are the steps to be followed:
1. Conda needs to be installed as followed (if not already installed)

        srun --partition=debug --mem=50GB -c 4 --pty bash
        cd $WORK
        bash /data/humangen_mouse/scpipeline/installation/Anaconda3-2021.05-Linux-x86_64.sh
        # Provide the installation path to "$WORK/.omics/anaconda3/"

2. Next a conda environment called scVelocity needs to be created, as follows

        # Navigate to the folder where you downloaded this repository
        srun --partition=debug --mem=50GB -c 4 --pty bash
        conda env create -f installation/scVelocity.yml
        exit
        sbatch installation/install_packages.sh

3. For trajectory analysis, a separate conda environment called scTrajectory needs to be created, as follows

        # Navigate to the folder where you downloaded this repository
        srun --partition=debug --mem=50GB -c 4 --pty bash
        conda env create -f installation/scTrajectory.yml
        exit
>Tip: Try mamba instead of conda to speed up the installations quite a bit. Check the [documentation](https://mamba.readthedocs.io/en/latest/index.html) and the installation [here](https://github.com/conda-forge/miniforge#install).

## Using the scripts
Download the contents corresponding to each of the pipeline(s) along with the scripts located in the folder *src*. Copy it to anywhere on the OMICS cluster (slurm-based cluster). Edit the analysis parameters in the corresponding *xxx_params.yaml* file. Submit the following command:

    ./xxx_main.sh

*Note: you will have to run the command ```chmod +x xxx_main.sh``` the first time, to convert this file into an executable.*

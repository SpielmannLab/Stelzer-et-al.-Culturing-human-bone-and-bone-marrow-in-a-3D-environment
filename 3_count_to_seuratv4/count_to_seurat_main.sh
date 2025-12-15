#! /bin/bash
# Run this script to create a Seurat object from count matrix created by cell ranger
# It also performs velocyto. The spliced and unspliced count matrices are also added to the Seurat object
# STEP1. Provide all the details in the file "count_to_seurat_params.yaml"
# STEP2. Close this file and execute it using <./count_to_seurat_main.sh>. 
# Note: You may have to run the command <chmod +x count_to_seurat_main.sh> before hand.

sbatch --job-name=count_to_seurat count_to_seurat_core.sh
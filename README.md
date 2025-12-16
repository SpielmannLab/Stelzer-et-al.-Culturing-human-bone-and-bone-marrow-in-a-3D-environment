Following are the scripts used to analyse the single-cell data, as presented in the manuscript

Run /1_mkfastq/mkfastq_cellplex_main.sh to generate fastq files from Illumina sequencer output using bcl2fastq2 v2.20  
Run /2_cellranger/multi_cellplex_main.sh to generate gene expression count matrices from sequencing fastq files using cellranger 7.0.0  
Run /3_count_to_seuratv4/count_to_seurat_main.sh to generate a seurat v4 object from count matrix  
Run /4_sc_p_sample/sc_p_sample_main.sh to perform a standard analysis of the data, individually. This wrapper calls a Nextflow script sc_p_sample.nf, which includes multiple steps including importing the data, QC, filtering, clustering, dimensionality reduction, and DE gene analysis. This is run in the conda environment scVelocity.yml  
Run /5_sc_multi_sample/sc_multi_sample_main.sh to integrate individual seurat objects using harmony. This wrapper calls a Nextflow script sc_multi_sample.nf, which combines the data, performs dimensionality reduction, clustering and DE gene analysis. This is run in the conda environment scVelocity.yml

At this stage seurat object meta columns Annotation, Sample Type, and Sample were added using Seurat v4 and standard R packages  

Run /6_make_plots/make_plots_sbatch.sh to perform custom analysis used in this manuscript as well as generate figures presented in the manuscript.

[1] Absolute paths are used atmake certain instances, which will need to be adapted, as needed.
[2] Anaconda was used to set up the necessary environments

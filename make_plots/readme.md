# Available Meta Analysis

The scripts included in this repository needs to be downloaded in a slurm-based cluster (UzL OMICS cluster) and run using the xxx_main.sh script. You will have to convert the xxx_main.sh script into an executable, using the command.

        chmod +x xxx_main.sh

Some of the newer scripts do not have the *_main.sh script. Instead, you just submit using 

    sbatch *_sbatch.sh

or if you want to resume your job to save lots of time, use:

    sbatch *_sbatch.sh -resume

## regulon_analysis
This is to analyse the gene regulatory networks (or regulons) activate in specific cellular populations. Currently implementing SCENIC.

## cell_abundance_analysis

This is set of pipelines that is usefull to perform cell abundance analysis.

## pairwise_de_genes

This is a nextflow based script used to calculate the pairwise DE genes between any two set (groups) of cells, as defined by a key-value(s) pair. Multiple pairwise comparisons can be submitted at a time for parallel computations. Note that it is recommended to use pseudo-bulk DE genes, but this requires multiple samples per comparison. As a backup, for projects with just N=1 experimental repeat, single-cell DE analysis option is also provided

## subset_seurat_object

This is a nextflow based script used to take an rds file containing a Seurat object, and subset it based on the values of a particular meta data column. Multiple subsetting can be submitted at the same time for parallel computation. Adjust the number of cores as required in the _subset_seurat_object_core.sh_ script.

## rename_clusters (e.g., cell type annotations)

This is a nextflow based script used to modify or create a new metadata column based on a conversion table, provided in the form of a TSV file, for example during cell annotation. The TSV file needs to have two columnes, the names of the column should correspond to the "from" and "to" names of the metaadata. For example, if the names of the columns in the TSV file are "RNA_snn_res.0.1" and "Annotations", then the script will create a new metadata column called "Annotations" and store the new names for the clusters in this column. Plots can from then on be named using this new column. In theory, you can rename any metadata column values.

## make_plots

This is a nextflow based script for making different types of plots (UMAP, cell composition, violin and dotplot).

## celltype annotation

Automatically annotate celltypes based on a database in the form of an excel sheet

## celltype prioritization

Identify the cell type that has the biggest molecular change between two groups. Kind of like DE genes, but does not give you the DE genes - only prioritizes cell types.

## ssGSEA

Perform single-sample Gene Set Enrichment Analysis on Seurat object. The output is a new Seurat object with additional metadata columns (one for every gene set tested). Requires a **separate conda environment**, [see here](ssGSEA/readme.md). A pseudo-bulk version has also been implemented, where the ssGSEA scores can be averaged across cells that have the same key-value pairs [see the yaml file](ssGSEA/readme.md)

## cellxgene

Export an analysed Seurat object (i.e., _.rds file) into a _.h5ad format (i.e., AnnData format), which can be explored using cellxgene viewer. It also contains details on how to install cellxgene to be able to view the \*.h5ad file.

## cell_cell_signalling

Tease out potential cell-cell signalling based on known ligand-receptor interactions. Currently implemented using the package LIANA and NICHES, also with an option to perform differential signalling using both of these packages. Check with Varun if you have questions.

## decoupleR

Performs biologival activity analysis using decoupleR tool. 3 following analysis will be carried out:<br />
1) Pathway activity inference using PROGENy<br />
2) Pathway activity inference using msigdb<br />
3) TF activity based on CollecTRI<br />

# Installation

The scripts **make_plots**, **modify_seurat_metadata**, **pairwise_de_genes**, **rename_clusters**, **cellxgene**, and **subset_seurat_object** require the conda environment _scVelocity_. Check out [scpipeline](https://github.com/SpielmannLab/scpipeline) repository to figure out how to get this set up. **ssGSEA** on the other hand requires the conda environment _ssGSEA_ and _ggplot_essentials_ - see [here](ssGSEA/readme.md).

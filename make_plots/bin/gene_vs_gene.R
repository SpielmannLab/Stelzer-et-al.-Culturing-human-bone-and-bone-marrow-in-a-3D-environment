#!/usr/bin/env Rscript
" Plots for a given Seurat object

Usage: gene_vs_gene.R --file_sc_obj=<file> --assay=<value> --gene_pair=<value> --group_by=<value> --pt_size=<value> --width=<value> --height=<value>

Options:
    -h --help               	Show this screen.
    --version              	00.99.01
    --file_sc_obj=<file>        	*.rds file containing the Seurat Object
    --assay=<value>           Which seurat assay to use. <RNA, SCT, or integrated> 
    --gene_pair=<value>            Pair of genes to plot against each other. Comma separated
    --group_by=<value>        Metadata column to group the plotting by.
    --pt_size=<value>   Size of the points in the violin plot
    --width=<value>     Width of the plot
    --height=<value>    Heigth of the plot

"-> doc

# --- Load all the libraries
suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)
file_sc_obj <- arguments$file_sc_obj
assay <- arguments$assay
gene_pair <- arguments$gene_pair %>%
    strsplit(split = ",") %>%
    unlist()
group_by <- arguments$group_by
pt_size <- as.numeric(arguments$pt_size)
width <- as.numeric(arguments$width)
height <- as.numeric(arguments$height)

## Read the Seurat object
sc_obj <- readRDS(file_sc_obj)

# Change the seurat assay for gene expression.
DefaultAssay(sc_obj) <- assay

# Deal with null values
if (group_by %in% c("null", "NULL") || !(group_by %in% colnames(sc_obj@meta.data))) {
    group_by <- NULL
}

# check if both genes in the gene pairs are present throw an error if none present
gene_pair <- intersect(gene_pair, c(rownames(sc_obj), colnames(sc_obj@meta.data)))
if (length(gene_pair) != 2 ) stop("At least one of the genes in the gene pairs not present in data")

plot <- FeatureScatter(sc_obj,
    feature1 = gene_pair[1],
    feature2 = gene_pair[2],
    group.by = group_by,
    pt.size = pt_size)

filename = paste0(gsub(file_sc_obj,
    pattern = ".rds",
    replacement = paste0("_gvg_", gene_pair[1], "_vs_", gene_pair[2], ".pdf")))
ggsave(plot = plot,
    filename = filename,
    width = width,
    height = height)
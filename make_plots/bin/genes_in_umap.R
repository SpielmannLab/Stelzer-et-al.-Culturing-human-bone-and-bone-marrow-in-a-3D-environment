#!/usr/bin/env Rscript
" Plots for a given Seurat object

Usage: genes_in_umap --file_sc_obj=<file> --genes=<value> --split_by=<value> --assay=<value> --genes_per_file=<value> --pt_size=<value> --cell_order=<value> --width=<value> --height=<value>

Options:
    -h --help               	Show this screen.
    --version              	00.99.01
    --file_sc_obj=<file>        	*.rds file containing the Seurat Object
    --assay=<value>           Which seurat assay to use. <RNA, SCT, or integrated> 
    --genes=<value>            Genes to plot. Comma separated
    --split_by=<value>        Metadata column to split the plots by.
    --genes_per_file=<value>  How many genes to plot per file
    --pt_size=<value>    Size of the dots. 2 is a good number
    --cell_order=<value>     To order the cells based on gene expression values. Use true for low expression
    --width=<value>      Width of the plot
    --height=<value>     Height of the plot

"-> doc

# --- Load all the libraries
suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)
file_sc_obj <- arguments$file_sc_obj
genes <- arguments$genes %>%
    strsplit(split = ",") %>%
    unlist()
split_by <- arguments$split_by
assay <- arguments$assay
genes_per_file <- as.numeric(arguments$genes_per_file)
pt_size <- as.numeric(arguments$pt_size)
cell_order <- as.logical(arguments$cell_order)
width <- as.numeric(arguments$width)
height <- as.numeric(arguments$height)

## Read the Seurat object
sc_obj <- readRDS(file_sc_obj)

# Change the seurat assay for gene expression.
DefaultAssay(sc_obj) <- assay

# Deal with null values
if (split_by %in% c("null", "NULL") || !(split_by %in% colnames(sc_obj@meta.data))) {
    split_by <- NULL
}

# subset genes that are only in the dataset and throw an error if none present
genes <- intersect(genes, c(rownames(sc_obj), colnames(sc_obj@meta.data)))
if (length(genes) == 0 ) stop("Genes/features requested are not present in the data")

# Make FeaturePlot
# convert genes to a list
genes <- split(genes, ceiling(seq_along(genes) / genes_per_file))

for (i in seq_len(length(genes))) {
    genes_subset <- genes[[i]]

    plot <- FeaturePlot(sc_obj,
        features = genes_subset,
        keep.scale = "feature",
        cols = c("#bababa", "#ca0020"),
        split.by = split_by,
        raster = TRUE,
        pt.size = pt_size,
        order = cell_order
    ) &
        theme(legend.position="none",
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(), 
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            line = element_blank(),
            plot.title = element_text(family="sans", face="plain", hjust=0.5, size=15),
            aspect.ratio = 1)

    filename = paste0(gsub(file_sc_obj,
        pattern = ".rds",
        replacement = paste0("_FeaturePlot_", i, ".pdf")))
    ggsave(plot = plot, filename = filename, width = width, height = height)
}
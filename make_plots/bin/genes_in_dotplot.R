#!/usr/bin/env Rscript
" Make dotplot for provided genes for the given Seurat object

Usage: genes_in_dotplot --file_sc_obj=<file> --genes=<value> --group_by=<value> --split_by=<value> --assay=<value> --width=<value> --height=<value> --angle_x_text=<value>

Options:
    -h --help               	Show this screen.
    --version              	00.99.01
    --file_sc_obj=<file>        	*.rds file containing the Seurat Object
    --genes=<value>            Genes to plot. Comma separated
    --group_by=<value>        Metadata column to group the plotting by. 
    --split_by=<value>        Metadata column to split the plotting by. 
    --assay=<value>           Which seurat assay to use. <RNA, SCT, or integrated> 
    --width=<value>     Width of the plot
    --height=<value>    Heigth of the plot
    --angle_x_text=<value>     Angle to rotate the x-axis labels by 0..90
"-> doc

# --- Load all the libraries
suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)
file_sc_obj <- arguments$file_sc_obj
genes <- arguments$genes %>%
    strsplit(split = ",") %>%
    unlist()
group_by <- arguments$group_by
assay <- arguments$assay
split_by <- arguments$split_by
width <- as.numeric(arguments$width)
height <- as.numeric(arguments$height)
angle_x_text <- as.numeric(arguments$angle_x_text)

## Read the Seurat object
sc_obj <- readRDS(file_sc_obj)

# Change the seurat assay for gene expression.
DefaultAssay(sc_obj) <- assay

# Deal with null values
if (group_by %in% c("null", "NULL") || !(group_by %in% colnames(sc_obj@meta.data))) {
    message("Available group_by in the dataset are: ", paste(colnames(sc_obj@meta.data), collapse=";  "))
    stop("Dotplots requires the definition of <group_by>")
}

# subset genes that are only in the dataset and throw an error if none present
genes <- intersect(genes, c(rownames(sc_obj), colnames(sc_obj@meta.data)))
if (length(genes) == 0 ) stop("Genes/features requested are not present in the data")


# Deal with if a valid split_by is provided:
if ((split_by %in% c("null", "NULL")) || !(split_by %in% colnames(sc_obj@meta.data)) || (split_by == group_by)) {
    plot <- DotPlot(sc_obj, features = genes, group.by = group_by)
    # rotate the x-axis tick text if asked.
    if (angle_x_text > 0) {
        plot <- plot +  theme(axis.text.x = element_text(angle = angle_x_text, hjust = 1, vjust = 0.5))
        }
} else if (split_by %in% colnames(sc_obj@meta.data)) {
    Idents(sc_obj) <- sc_obj[[split_by]]
    plot_list <- lapply(X = levels(sc_obj), FUN = function(ident){
        plot <- DotPlot(sc_obj, features = genes, group.by = group_by, idents = ident) +
            ggtitle(ident)
        # rotate the x-axis tick text if asked.
        if (angle_x_text > 0) {
            plot <- plot +  theme(axis.text.x = element_text(angle = angle_x_text, hjust = 1, vjust = 0.5))
        }
    })
    plot <- plot_grid(plotlist = plot_list, nrow = 1)
}

filename = paste0(gsub(file_sc_obj,
    pattern = ".rds",
    replacement = paste0("_DotPlot.pdf")))
ggsave(plot = plot,
    filename = filename,
    width = width,
    height = height)

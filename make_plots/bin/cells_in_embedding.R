#!/usr/bin/env Rscript
" Make a cell embedding plot (umap, tsne, pca)

Usage: code.R --file_sc_obj=<value> --group_by=<value> --color_values=<value> --split_by=<value> --pt_size=<value> --shuffle_pts=<value> --width=<value> --height=<value> --reduction_to_use=<value> --aspect_ratio=<value>

Options:
    -h --help               	Show this screen.
    --file_sc_obj=<value>        	*.rds file containing the Seurat Object
    --group_by=<value>    	Metadata column-name to color the cells by
    --color_values=<value>  OPTIONAL. Comma-separated color values. Eg. <grey,blue,...> or <#808080,#0000FF,...>
    --split_by=<value>    	Metadata column-name to split the plots by
    --pt_size=<value>    Size of the dots. 2 is a good number
    --shuffle_pts=<value>    Default is TRUE. Whether to shuffle the points for better visibility
    --width=<value>      Width of the plot
    --height=<value>     Height of the plot
    --reduction_to_use=<value>    Which embedding to use? umap, tsne or pca
    --aspect_ratio=<value>    Aspect ratio of the plot excluding all labels etc.

" -> doc

# --- Load all the libraries
suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))

# --- Read all the arguments passed
arguments <- docopt(doc)
list2env(x = arguments, envir = environment())
pt_size <- as.numeric(pt_size)
width <- as.numeric(width)
height <- as.numeric(height)
aspect_ratio <- as.numeric(aspect_ratio)

# If split_by is not provided, set it to NULL
if (split_by %in% c("null", "NULL")) {
    split_by <- NULL
}

# If shuffle_pts is not provided, set it to TRUE, else, convert it to logical
if (shuffle_pts %in% c("null", "NULL")) {
    shuffle_pts <- TRUE
} else {
    shuffle_pts <- as.logical(shuffle_pts)
}

if (color_values %in% c("NULL", "null")) {
    color_values <- NULL
}
if (!is.null(color_values)){
    color_values <- color_values %>%
        strsplit(split = ",") %>%
        unlist()
}

# Remove unwanted variables and print environmental variables
rm(doc, arguments, help)
# Print all the environment variables
trash <- lapply(X = ls(), FUN = function(x) {message(x, ": ", get(x))})

## Read the Seurat object
sc_obj <- readRDS(file_sc_obj)

## Make the plot
plot <- DimPlot(sc_obj,
        group.by = group_by,
        cols = color_values,
        split.by = split_by,
        shuffle = shuffle_pts,
        raster = FALSE,
        pt.size = pt_size,
        reduction = reduction_to_use) +
    theme_void() +
    theme(legend.position = "right",
        aspect.ratio = aspect_ratio)

filename = gsub(basename(file_sc_obj), pattern = ".rds", replacement = paste0("_", reduction_to_use, ".pdf"))
ggsave(plot = plot, filename = filename, width = width, height = height)

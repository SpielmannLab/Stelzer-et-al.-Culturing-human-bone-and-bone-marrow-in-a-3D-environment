#!/usr/bin/env Rscript
doc <- "Ridge plots for a given Seurat object
Usage: genes_in_ridgeplot.R --file_sc_obj=<file> --genes=<value> --group_by=<value> [--split_by=<value>] [--sort_splits_by=<value>] --assay=<value> [--color_values=<value>] --bandwidth=<value> --genes_per_file=<value> --width=<value> --height=<value>

Options:
-h --help               Show this screen.
--version               00.99.01
--file_sc_obj=<file>    *.rds file containing the Seurat Object
--genes=<value>         Genes to plot. Comma separated
--group_by=<value>      Metadata column to group the plotting by.
--split_by=<value>      Metadata column to split the plots. A column per value
--sort_splits_by=<value>  OPTIONAL. <ascending,descending,ordered-comma-separated-values>
--assay=<value>         Seurat assay to use. <RNA, SCT, or integrated>
--color_values=<value>  OPTIONAL. Comma-separated values <grey,#808080,.>
--bandwidth=<value>     Bandwidth of ridgeplot. Use 0.1 for gene expression
--genes_per_file=<value>How many genes to plot per file
--width=<value>         Width of the plot
--height=<value>        Height of the plot
"

# --- Load all the libraries
suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(future))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggridges))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args = TRUE)
print(arguments)
file_sc_obj <- arguments$file_sc_obj
genes <- arguments$genes %>%
  strsplit(split = ",") %>%
  unlist()
group_by <- arguments$group_by
split_by <- arguments$split_by
sort_splits_by <- arguments$sort_splits_by %>%
  strsplit(split = ",") %>%
  unlist()
assay <- arguments$assay
color_values <- arguments$color_values
bandwidth <- arguments$bandwidth
genes_per_file <- as.numeric(arguments$genes_per_file)
width <- as.numeric(arguments$width)
height <- as.numeric(arguments$height)

## Read the Seurat object
sc_obj <- readRDS(file_sc_obj)

# Define functions
CustomRidgePlot <- function(
    sc_obj = sc_obj, features = genes_subset, group_by = group_by,
    cols = color_values) {
  cts <- GetAssayData(sc_obj, slot = "data")
  df_cts <- data.frame(Matrix::t(cts[features, , drop = FALSE])) %>%
    tibble::rownames_to_column(var = "bc")
  df_ridges <- sc_obj@meta.data %>%
    dplyr::select(all_of(c(group_by, split_by))) %>%
    tibble::rownames_to_column(var = "bc") %>%
    full_join(df_cts, by = join_by(bc))
  plots <- lapply(features, FUN = function(gene) {
    plot <- ggplot(df_ridges, aes(x = .data[[gene]], y = .data[[group_by]], fill = .data[[group_by]])) +
      geom_density_ridges(scale = 7, rel_min_height = 0, bandwidth = bandwidth) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = expansion(mult = c(0, 0.9))) +
      labs(title = gene) +
      xlab("Expression") +
      theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey92"),
        axis.ticks = element_line(color = "grey92"))
    if (!is.null(cols)) {
      plot <- plot + scale_fill_manual(values = cols)
    }
    return(plot)
  })
  names(plots) <- features
  return(plots)
}

# Change the seurat assay for gene expression.
DefaultAssay(sc_obj) <- assay

# Deal with null values
if (is.character(group_by)) {
  if (group_by %in% c("null", "NULL") || !(group_by %in% colnames(sc_obj@meta.data))) {
    group_by <- NULL
  }
}
if (is.character(split_by)) {
  if (split_by %in% c("null", "NULL") || !(split_by %in% colnames(sc_obj@meta.data))) {
    split_by <- NULL
  }
}
if (is.character(sort_splits_by)) {
  if (sort_splits_by %in% c("NULL", "null")) {
    sort_splits_by <- NULL
  } else {
    sort_splits_by <- sort_splits_by %>%
      strsplit(split = ",") %>%
      unlist()
  }
}
if (is.character(color_values)) {
  if (color_values %in% c("NULL", "null")) {
    color_values <- NULL
  } else {
    color_values <- color_values %>%
      strsplit(split = ",") %>%
      unlist()
  }
}

# subset genes that are only in the dataset and throw an error if none present
genes <- intersect(genes, c(rownames(sc_obj), colnames(sc_obj@meta.data)))
if (length(genes) == 0) stop("Genes/features requested are not present in the data")

# convert genes to a list
genes <- split(genes, ceiling(seq_along(genes) / genes_per_file))

if (is.null(split_by)) {
  # If no splitting necessary
  for (i in seq_len(length(genes))) {
    genes_subset <- genes[[i]]

    plot <- RidgePlot(sc_obj,
      features = genes_subset, group.by = group_by, cols = color_values,
      combine = FALSE
    )
    filename <- paste0(gsub(file_sc_obj, pattern = ".rds", replacement = paste0(
      "_RidgePlot_",
      i, ".pdf"
    )))
    ggsave(
      plot = plot_grid(plotlist = plot, ncol = ceiling(sqrt(length(genes_subset)))),
      filename = filename, width = width, height = height
    )

    # Generate custom ridgeplot
    plot_custom <- CustomRidgePlot(sc_obj,
      features = genes_subset, group_by = group_by,
      cols = color_values
    )
    filename <- paste0(gsub(file_sc_obj, pattern = ".rds", replacement = paste0(
      "_CustomRidgePlot_",
      i, ".pdf"
    )))
    ggsave(
      plot = plot_grid(plotlist = plot_custom, ncol = ceiling(sqrt(length(genes_subset)))),
      filename = filename, width = width, height = height
    )
  }
} else {
  # If asked to split by some key
  sc_list <- SplitObject(sc_obj, split.by = split_by)
  # Make sure sort_splits_by is one of the possible options and deal with it
  sc_list <- if (is.null(sort_splits_by)) {
    sc_list
  } else if (all(sort_splits_by %in% sc_obj@meta.data[[split_by]])) {
    sc_list[sort_splits_by]
  } else if (sort_splits_by == "ascending") {
    sc_list[sort(names(sc_list))]
  } else if (sort_splits_by == "descending") {
    sc_list[rev(sort(names(sc_list)))]
  } else {
    stop("The sort_splits_by value does not seem right")
  }
  for (i in seq_len(length(genes))) {
    genes_subset <- genes[[i]]

    plots <- lapply(sc_list, FUN = function(sc_obj) {
      plot <- RidgePlot(sc_obj,
        features = genes_subset, group.by = group_by,
        cols = color_values, combine = FALSE
      )
      # set title for every plot in the list of plots of every gene
      plot <- lapply(plot, function(p) {
        current_title <- p$labels$title
        split_value <- unique(sc_obj@meta.data[[split_by]])
        p <- p & ggtitle(paste0(split_value, ":", current_title))
      })
      return(plot)
    })

    plots <- unlist(plots, recursive = FALSE, use.names = TRUE)
    filename <- paste0(gsub(file_sc_obj, pattern = ".rds", replacement = paste0(
      "_RidgePlot_",
      i, ".pdf"
    )))
    ggsave(plot = plot_grid(
      plotlist = plots, ncol = n_distinct(sc_obj@meta.data[[split_by]]),
      byrow = FALSE
    ), filename = filename, width = width, height = height)

    # Generate custom ridgeplot
    plots_custom <- lapply(sc_list, FUN = function(sc_obj) {
      plot <- CustomRidgePlot(sc_obj,
        features = genes_subset, group_by = group_by,
        cols = color_values
      )
      # set title for every plot in the list of plots of every gene
      plot <- lapply(plot, function(p) {
        current_title <- p$labels$title
        split_value <- unique(sc_obj@meta.data[[split_by]])
        p <- p & ggtitle(paste0(split_value, ":", current_title))
      })
      return(plot)
    })
    plots_custom <- unlist(plots_custom, recursive = FALSE, use.names = TRUE)

    filename <- paste0(gsub(file_sc_obj, pattern = ".rds", replacement = paste0(
      "_CustomRidgePlot_",
      i, ".pdf"
    )))
    ggsave(plot = plot_grid(
      plotlist = plots_custom, ncol = n_distinct(sc_obj@meta.data[[split_by]]),
      byrow = FALSE
    ), filename = filename, width = width, height = height)
  }
}

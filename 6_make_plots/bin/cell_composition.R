#!/usr/bin/env Rscript
" Plots for a given Seurat object

Usage: code.R --file_sc_obj=<value> --grouping_var=<value> --xvar=<value> --correct_grouping_var_imbalance=<value> --figure_height=<value>

Options:
  -h --help               	Show this screen.
  --file_sc_obj=<value>        	*.rds file containing the Seurat Object
  --grouping_var=<value>    	Metadata column-name to stratify the cell composition by
  --xvar=<value>        Metadata column-name to be on the x-axis in the cell composition plot
  --correct_grouping_var_imbalance=<value>  TRUE or FALSE to correct for grouping_var imbalance
  --figure_height=<value>   Height of the figure

" -> doc

# --- Load all the libraries
suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args = TRUE)
file_sc_obj <- arguments$file_sc_obj
grouping_var <- arguments$grouping_var
xvar <- arguments$xvar
correct_grouping_var_imbalance <- as.logical(arguments$correct_grouping_var_imbalance)
figure_height <- as.numeric(arguments$figure_height)

## Read the Seurat object
sc_obj <- readRDS(file_sc_obj)

# Make plots only if grouping_var and xvar are found in the sc_obj metadata. Else, just output the head of the metadata
if ((grouping_var %in% colnames(sc_obj@meta.data)) &&  (xvar %in% colnames(sc_obj@meta.data))) {
    # Total count table
    count_table <- sc_obj@meta.data %>%
        group_by(across(any_of(c(xvar)))) %>%
        summarise(count = n()) %>%
        mutate(count_k=count/1000)

    plot1 <- ggplot(count_table,
    aes(x = .data[[xvar]],
        y=count_k)) +
    geom_col(fill="black") +
    xlab(xvar) +
    ylab("Total cells per cluster x 1000") +
    coord_flip()

    # Count table grouped by grouping_var
    if (!correct_grouping_var_imbalance) {
        count_table_grouped <- sc_obj@meta.data %>%
            group_by(across(any_of(c(xvar, grouping_var)))) %>%
            summarise(count = n()) %>%
            mutate(perc = count/sum(count))
        yaxis_label <- "% cells, not corrected"
    } else if (correct_grouping_var_imbalance) {
        count_table_grouped <- sc_obj@meta.data %>%
            group_by(across(any_of(c(grouping_var, xvar)))) %>%
            summarise(count = n()) %>%
            mutate(perc = count/sum(count))
        yaxis_label <- "% cells, corrected"
    }

    # save count table
    filename = gsub(basename(file_sc_obj),
        pattern = ".rds",
        replacement = "_count_table_grouped.tsv")
    write.table(count_table_grouped, file = filename, sep = "\t", row.names = FALSE)

    # Make and save plot
    plot2 <- ggplot(count_table_grouped,
    aes(x = .data[[xvar]],
        fill = .data[[grouping_var]],
        y = perc)) +
    geom_col(position = "fill") +
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    ylab(yaxis_label) +
    coord_flip() +
    xlab(NULL)

    filename = gsub(basename(file_sc_obj),
            pattern = ".rds",
            replacement = "_cell_composition.pdf")
    ggsave(plot = plot_grid(plot1 + theme_bw(),
            plot2 + theme_bw(),
            nrow = 1),
        filename = filename,
        width = 10, height = figure_height)
} else {
    table <- sc_obj@meta.data %>%
        head
    filename <- gsub(basename(file_sc_obj),
        pattern = ".rds",
        replacement = "_available_metadata.tsv")
    write.table(table, file = filename, sep = "\t", row.names = FALSE)
    warning("No grouping_var or xvar provided")
}
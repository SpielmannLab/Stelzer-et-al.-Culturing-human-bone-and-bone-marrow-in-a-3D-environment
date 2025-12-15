#!/usr/bin/env Rscript
" Plots for a given Seurat object

Usage: plots.R --file_sc_obj=<file> --genes=<value> --group_by=<value> --assay=<value> --color_values=<value> --genes_per_file=<value> --pt_size=<value> --width=<value> --height=<value> --do_statistics=<value>

Options:
    -h --help               	Show this screen.
    --version              	00.99.01
    --file_sc_obj=<file>        	*.rds file containing the Seurat Object
    --genes=<value>            Genes to plot. Comma separated
    --group_by=<value>        Metadata column to group the plotting by.
    --assay=<value>           Which seurat assay to use. <RNA, SCT, or integrated> 
    --color_values=<value>      OPTIONAL. Comma-separated color values. Eg. <grey,blue,...> or <#808080,#0000FF,...>
    --genes_per_file=<value>  How many genes to plot per file 
    --pt_size=<value>   Size of the points in the violin plot
    --width=<value>     Width of the plot
    --height=<value>    Heigth of the plot
    --do_statistics=<value>     Whether or not to do statistics. TRUE or FALSE

"-> doc

# --- Load all the libraries
suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(future))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)
file_sc_obj <- arguments$file_sc_obj
genes <- arguments$genes %>%
    strsplit(split = ",") %>%
    unlist()
group_by <- arguments$group_by
assay <- arguments$assay
color_values <- arguments$color_values 
genes_per_file <- as.numeric(arguments$genes_per_file)
pt_size <- as.numeric(arguments$pt_size)
width <- as.numeric(arguments$width)
height <- as.numeric(arguments$height)
do_statistics <- as.logical(arguments$do_statistics)

## Read the Seurat object
sc_obj <- readRDS(file_sc_obj)

# Change the seurat assay for gene expression.
DefaultAssay(sc_obj) <- assay

# Deal with null values
if (group_by %in% c("null", "NULL") || !(group_by %in% colnames(sc_obj@meta.data))) {
    group_by <- NULL
}

if (color_values %in% c("NULL", "null")) {
    color_values <- NULL
}
if (!is.null(color_values)){
    color_values <- color_values %>%
        strsplit(split = ",") %>%
        unlist()
}

# Create a custome violin plot routine for smaller and transparent points
custom_VlnPlot <- function(sc_obj = sc_obj, group.by = "genotype", color_values = NULL, features, slot = "data", do_statistics = TRUE, pt_size = 0.1) {

    plot <- list()
    stat.test <- list()
    for (feature in features){

        #Get gene expression data
        gex_data <- FetchData(sc_obj, vars = feature, slot = slot)
        gex_ident_data <- cbind(gex_data, sc_obj[[group.by]])
        colnames(gex_ident_data) <- c("expression", "group")

        plot[[feature]] <- ggplot(gex_ident_data,
            aes(x = group, y = expression)) +
            geom_violin(aes(fill = group))

        # if color values for the violin provided, use them
        if (!is.null(color_values)){
            plot[[feature]] <- plot[[feature]] + 
                scale_fill_manual(values = color_values)
        }
         
        plot[[feature]] <- plot[[feature]] +
            geom_jitter(size = pt_size, shape = 16, fill = "black", width = 0.2, alpha = 0.1, height = 0) +
            theme_classic() +
            ggtitle(feature) +
            theme(text = element_text(family = "sans", size = 12),
                plot.title = element_text(size = 10, face = "italic", hjust = 0.5),
                axis.title = element_blank(),
                axis.text = element_text(family = "sans", size = 12),
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                legend.position = "none",
                aspect.ratio = 0.8)

        # Add statistics, if asked requested
        if (do_statistics) {
            
            # If zero inflated values, skip
            total_exp <- gex_ident_data %>%
                group_by(group) %>%
                summarize(total_exp = sum(expression)) %>%
                tibble::column_to_rownames(var = "group")

            # find non-zero group pair combinations. Error because of 0 to 0 wilcox test. Remove all such 0-0 combinations.
            # first get a list of combinations
            combinations <- combn(as.character(rownames(total_exp)), m = 2)
            # convert_to list
            combinations_list <- lapply(X = seq_len(ncol(combinations)), FUN = function(x) {
                combination  <- combinations[, x]
                if (sum(total_exp[combination, ]) == 0) {
                    return()
                } else {
                    return(combination)
                }
            })
            # remove all those which has 0-0 test.
            wilcox_test_combinations <- combinations_list %>%
                purrr::compact()

            # Do statistical test
            y_step_increase <- 0.5
            stat.test[[feature]] <- gex_ident_data %>%
                wilcox_test(formula = expression ~ group,
                    p.adjust.method = "bonferroni",
                    comparisons = wilcox_test_combinations) %>%
                add_y_position(step.increase = y_step_increase)
            
            if ("p.adj.signif" %in% colnames(stat.test[[feature]])) {
                label_by <- "p.adj.signif"
            } else {
                label_by <- "p"
            }

            # Add the statistics to the test.
            plot[[feature]] <- plot[[feature]] +
                stat_pvalue_manual(stat.test[[feature]],
                    label = label_by,
                    hide.ns = TRUE,
                    tip.length = 0.01)

            stat.test[[feature]]$feature <- feature # Add a column to the stat.test data.table
        }
    }
    return(list(stat.test = stat.test, plot=plot))
}

# subset genes that are only in the dataset and throw an error if none present
genes <- intersect(genes, c(rownames(sc_obj), colnames(sc_obj@meta.data)))
if (length(genes) == 0 ) stop("Genes/features requested are not present in the data")

# convert genes to a list
genes <- split(genes, ceiling(seq_along(genes) / genes_per_file))

for (i in seq_len(length(genes))) {

    genes_subset <- genes[[i]]

    plot_n_stat <- custom_VlnPlot(sc_obj,
        features = genes_subset,
        group.by = group_by,
        color_values = color_values,
        do_statistics = do_statistics,
        pt_size = pt_size)

    plot <- plot_n_stat$plot

    if (do_statistics) {
        stat.test <- plot_n_stat$stat.test %>%
        do.call(what=rbind) %>% # convert from a list to a long data.table
        as.data.frame() %>%
        mutate(groups = paste(unlist(groups), collapse = ",")) #collapse list into a string

        filename = paste0(gsub(file_sc_obj,
            pattern = ".rds",
            replacement = paste0("_StatTest_", i, ".tsv")))
        write.table(stat.test, file = filename, sep = "\t", row.names = FALSE)
    }

    filename = paste0(gsub(file_sc_obj,
        pattern = ".rds",
        replacement = paste0("_ViolinPlot_", i, ".pdf")))
    ggsave(plot = plot_grid(plotlist = plot, ncol = ceiling(sqrt(length(genes_subset)))),
        filename = filename,
        width = width,
        height = height)
}

quit("no")

# Old script from Saranya. Has the splitting function
# Make Violin Plot
plots <- list()

for (i in seq(1:length(features))) {
    g[[i]] <- VlnPlot(sc_obj,
        features = features[i],
        split.by = condition,
        group.by = meta_col,
        split.plot = TRUE,
        pt.size = 1)

    # Do statistical testing. 
    exp <- rlang::expr(!! sym(features[i]) ~ split)
    stat.test <- g[[i]]$data %>% 
        group_by(ident) %>% 
        t_test(formula = eval(exp)) %>% 
        adjust_pvalue(method="bonferroni") %>% 
        add_significance("p.adj")
    write.table(stat.test,paste0(samplename,"-",features[i],"_VlnPlot_stats.tsv"), sep='\t')
    print(g[[i]])

}
do.call(grid.arrange,  g)
dev.off()

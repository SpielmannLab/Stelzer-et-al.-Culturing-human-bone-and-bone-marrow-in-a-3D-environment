"Take a monocle cell data set object, and perform trajectory analysis

Usage: trajectory.R --file_cds_obj=<file> --filename_prefix=<value> --clustering_k=<value> --seperate_trajectory_by_partition=<value> --close_loop=<value> --group_bys=<value> --pt_size=<value> --width=<value> --height=<value>  --do_de_genenes_in_trajectory=<value>
Options:
    -h --help			Show this screen.
    --file_cds_obj=<file>		The rds file of a monocle object (CDS) to perform trajectory analysis
    --filename_prefix=<value>		Optional. Provide a string to prefix all filenames. If <NA>, uses unique value from annotations
    --clustering_k=<value>        K parameter for clustering
    --seperate_trajectory_by_partition=<value>        Whether trajectories should be disconnected between partitions
    --close_loop=<value>      Whether circular trajectories are allowed
    --group_bys=<value>       The metadata column name to group_by
    --pt_size=<value>     Size of the points in the trajectory plots
    --width=<value>   Size of the trajectory plots
    --height=<value>      Height of the trajectory plots
    --do_de_genenes_in_trajectory=<value>        Whether or not to find genes that show seggregated expression in the UMAP plot
" -> doc

# Load libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

# --- Parameters: Read
arguments <- docopt::docopt(doc, quoted_args = TRUE)

file_cds_obj <- arguments$file_cds_obj
filename_prefix <- arguments$filename_prefix
clustering_k <- arguments$clustering_k %>%
    as.numeric()
seperate_trajectory_by_partition <- arguments$seperate_trajectory_by_partition %>%
    as.logical()
close_loop <- arguments$close_loop %>%
    as.logical()
group_bys <- arguments$group_bys %>%
    strsplit(split = ",") %>% 
    unlist()
pt_size <- arguments$pt_size %>%
    as.numeric()
width <- arguments$width %>%
    as.numeric()
height <- arguments$height %>%
    as.numeric()
do_de_genenes_in_trajectory <- arguments$do_de_genenes_in_trajectory %>%
    as.logical()

if ((filename_prefix %in% c("null", "NULL") || is.null(filename_prefix))) {
    filename_prefix <- gsub(file_cds_obj, pattern = "*.rds", replacement = "")
}

message("file_cds_obj: ", file_cds_obj)
message("filename_prefix: ", filename_prefix)
message("clustering_k: ", clustering_k)
message("seperate_trajectory_by_partition: ", seperate_trajectory_by_partition)
message("close_loop: ", close_loop)
message("group_bys: ", group_bys)
message("pt_size: ", pt_size)
message("width: ", width)
message("height: ", height)
message("do_de_genenes_in_trajectory: ", do_de_genenes_in_trajectory)

# read the monocle object file
cds <- readRDS(file_cds_obj) #This file should contain PCA and UMAP embeddings

# Cluster and partition the cells
cds <- cluster_cells(cds, reduction_method = "UMAP", k = clustering_k, random_seed=111)

# Make trajectory
seed.use <- 111
cds <- learn_graph(cds,
    learn_graph_control = list(minimal_branch_len = 20),
    use_partition = seperate_trajectory_by_partition,
    close_loop = close_loop)

# Make a figure to help choose the root node
plots <- lapply(X = group_bys, FUN = function(group_by) {
    plot_cells(cds,
        color_cells_by = group_by,
        label_cell_groups = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        graph_label_size = 2,
        cell_size = pt_size) +
        ggtitle(group_by) +
        theme_void() +
        theme(aspect.ratio = 1)
})

plot_with_principle_points <- plot_cells(cds,
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 4,
    cell_size = pt_size,
    label_principal_points = TRUE) +
    ggtitle("Root Nodes") +
    theme_void() +
    theme(aspect.ratio = 1)

plots[[length(plots) + 1]] <- plot_with_principle_points #append to the plot list. append function does not work


filename <- paste0(filename_prefix, "_trajectory", ".png")
ncol <- length(plots) %>%
    sqrt() %>%
    ceiling()
ggsave(plot <- plot_grid(plotlist = plots, ncol = ncol),
    filename = filename,
    width = width,
    height = height)

# save the analysed monocle3 cds object
filename <- paste0(filename_prefix, "_traj", ".rds")
saveRDS(cds, file = filename)

# Get genes that show localized distribution along the trajectories, using graph_test
# filter based on the qvalue, whether it is in highly variable genes and high morans_i coefficient.
# Do this only if no. cells is less than 10,000 cells to avoid error
# checkout this discussion thread: https://github.com/cole-trapnell-lab/monocle3/issues/512
if (do_de_genenes_in_trajectory && (dim(cds)[2] < 10000)) {
    de_genes_in_trajectory <- graph_test(cds, neighbor_graph = "principal_graph")
    filename <- paste0(filename_prefix, "_markers_in_trajectory", ".tsv")
    # save the file for plotting using figures_trajectory.R
    write.table(de_genes_in_trajectory, file = filename, sep = "\t", row.names = FALSE)
}

quit("no")

# get DE genes of cluters

markers <- top_markers(cds, group_cells_by = "cluster")
filename <- paste0(filename_prefix, "_markers", ".tsv")
# markers$cell_group <- as.numeric(markers$cell_group)
markers_selected <- markers %>%
    filter(fraction_expressing > 0.2) %>%
    arrange(cell_group, marker_test_q_value) %>%
    data.frame()
write.table(markers_selected, file = filename, sep = "\t", row.names = FALSE)

# make a dotplot
goi <- markers_selected %>%
    group_by(cell_group) %>%
    arrange(marker_test_q_value) %>%
    slice_head(n = 2) %>%
    pull(gene_id) %>%
    unique()
filename <- paste0(filename_prefix, "_markers", ".png")
plot <- plot_genes_by_group(cds, goi, group_cells_by = "cluster", ordering_type = "none")
ggsave(plot, filename = filename, width = 7, height = 7)
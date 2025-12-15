"Make pseudotime figures. Note trajectory.R must be run prior

Usage: pseudotime.R --file_cds_obj=<file> --filename_prefix=<value> --root_node=<value> --root_metadata_key=<value> --root_metadata_val=<value> --group_bys=<value> --pt_size=<value> --width=<value> --height=<value> --genes=<value> --gex_genes_per_file=<value> --gex_pt_size=<value> --gex_width=<value> --gex_height=<value>
Options:
    -h --help			            Show this screen.
    --file_cds_obj=<file>		    The rds file of a monocle object (CDS) to perform trajectory analysis
    --filename_prefix=<value>		Optional. Provide a string to prefix all filenames. If <NA>, uses unique value from annotations
    --root_node=<value>           If known, the root node to set as pseudotime=0
    --root_metadata_key=<value>    Metadata column for setting time zero in pseudotime
    --root_metadata_val=<value>    This cells with this metadata value will be set as time zero in pseudotime
    --group_bys=<value>           The metadata column name to group_by
    --pt_size=<value>             Size of the points in the trajectory plots
    --width=<value>               Size of the trajectory plots
    --height=<value>              Height of the trajectory plots
    --genes=<value>                 Genes to plot the expression values as a function of pseudotime
    --gex_genes_per_file=<value>    No of genes to plot per file
    --gex_pt_size=<value>           Size of dots in the gex vs pseudotime plot
    --gex_width=<value>             Size of the gex vs pseudotime plot
    --gex_height=<value>            Size of the gex vs pseudotime plot
" -> doc

# Load libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(viridis))

# --- Parameters: Read
arguments <- docopt::docopt(doc, quoted_args=TRUE)

file_cds_obj <- arguments$file_cds_obj
filename_prefix <- arguments$filename_prefix
root_node <- arguments$root_node
root_metadata_key <- arguments$root_metadata_key
root_metadata_val <- arguments$root_metadata_val
group_bys <- arguments$group_bys %>%
    strsplit(split = ",") %>% 
    unlist()
pt_size <- arguments$pt_size %>%
    as.numeric()
width <- arguments$width %>%
    as.numeric()
height <- arguments$height %>%
    as.numeric()
genes <- arguments$genes %>%
    strsplit(split = ",") %>% 
    unlist()
gex_genes_per_file <- as.numeric(arguments$gex_genes_per_file)
gex_pt_size <- as.numeric(arguments$gex_pt_size)
gex_width <- as.numeric(arguments$gex_width)
gex_height <- as.numeric(arguments$gex_height)


if ((filename_prefix %in% c("null", "NULL") || is.null(filename_prefix))) {
    filename_prefix <- gsub(file_cds_obj, pattern = "*.rds", replacement = "")
}

message("file_cds_obj: ", file_cds_obj)
message("filename_prefix: ", filename_prefix)
message("root_node: ", root_node)
message("root_metadata_key: ", root_metadata_key)
message("root_metadata_val: ", root_metadata_val)
message("group_bys: ", paste(group_bys, collapse = ", "))
message("pt_size: ", pt_size)
message("width: ", width)
message("height: ", height)
message("genes: ", paste(genes, collapse = ", "))
message("gex_genes_per_file: ", gex_genes_per_file)
message("gex_pt_size: ", gex_pt_size)
message("gex_width: ", gex_width)
message("gex_height: ", gex_height)

# read the monocle object file
cds <- readRDS(file_cds_obj) # rds file containing the cell data set

# Define function to get the root node for a particular group of cells definied by root_metadata_key/val pair
get_earliest_principal_node <- function(cds, root_metadata_key, root_metadata_val){
    cell_ids <- which(colData(cds)[, root_metadata_key] == root_metadata_val)
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
    return(root_pr_nodes)
}

# Fix root_node
if ((root_node %in% c("null", "NULL")) || is.null(root_node) || is.na(root_node)) {
    key_not_avail <- (root_metadata_key %in% c("null", "NULL")) || is.null(root_metadata_key) || is.na(root_metadata_key)
    val_not_avail <- (root_metadata_val %in% c("null", "NULL")) || is.null(root_metadata_val) || is.na(root_metadata_val)
    if (any(key_not_avail, val_not_avail)){
        stop("Neither root_node or root_metadata_key/val pair provided")
    } else {
        root_node <- get_earliest_principal_node(cds, root_metadata_key, root_metadata_val)
    }
}

message("root_node to be used: ", root_node)

# order cells by peudotime
cds <- order_cells(cds, root_pr_nodes = root_node)

# Make plots
plots1 <- lapply(X = group_bys, FUN = function(group_by) {
    plot_cells(cds,
        color_cells_by = group_by,
        label_cell_groups = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        graph_label_size = 4,
        cell_size = pt_size) +
        ggtitle(group_by) +
        theme_void() +
        theme(aspect.ratio = 1)
})

plots2 <- lapply(X = c("cluster", "pseudotime"), FUN = function(group_by) {
    plot_cells(cds,
        color_cells_by = group_by,
        label_cell_groups = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        graph_label_size = 4,
        cell_size = pt_size) +
        ggtitle(group_by) +
        theme_void() +
        theme(aspect.ratio = 1)
})

#rasterize the plots before saving it to pdf
# suppressPackageStartupMessages(library(ggrastr))
# p <- lapply(p, function(plot){
#   plot <- rasterize(plot, layers="Point", dpi=300)
#   return(plot)
# })

plots <- c(plots2, plots2) # combine all the plots

filename <- paste0(filename_prefix, "_pseudotime",  ".png")
ncol <- length(plots) %>%
    sqrt() %>%
    ceiling()
ggsave(plot <- plot_grid(plotlist = plots, ncol = ncol),
    filename = filename,
    width = width,
    height = height)

# Plot the expression of requested genes along the pseudotime
genes <- arguments$genes %>%
    strsplit(split = ",") %>%
    unlist()
genes <- intersect(genes, rownames(cds))

if (length(genes) == 0) quit("no")

genes <- split(genes, ceiling(seq_along(genes) / gex_genes_per_file))

for (i in seq_len(length(genes))) {
    # Plot the variation in the expression of these top genes with pseudotime.
    cds_subset <- cds[genes[[i]], ]
    plot <- plot_genes_in_pseudotime(cds_subset,
        color_cells_by = group_bys[1],
        cell_size = gex_pt_size,
        min_expr = 0.5,
        nrow = ceiling(sqrt(gex_genes_per_file)),
        ncol = ceiling(sqrt(gex_genes_per_file))) +
        theme(legend.position = "bottom")
    filename <- paste0(filename_prefix, "_gex_in_pseudotime_", i, ".png")
    ggsave(plot=plot, filename=filename, width = gex_width, height = gex_height)
}
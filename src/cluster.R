#
# Cluster a merged seurat3 object using several resolution points
#

"  Cluster a merged seurat3 object using several resolution points

Usage: cluster.R --scobject=<file> --res=<value> --ncores=<value> [--task=<value>]

Options:
  -h --help            Show this screen.
  --scobject=<file>    Optional. If !is.null, do not use jobname to guess scfilename.
  --res=<value>        String separated by comma with resolution for clutering
  --ncores=<value>     Number of processors to use
  --task=<value>       Harmony or Seurat integration

Author:
  Cesar Prada
e-mail:
  prada@molgen.mpg.de

" -> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, quoted_args = TRUE)
print(arguments)

# -------------------------------
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(dplyr))
suppressMessages(library(parallel))
# ------------------------------

# ---------- Functions ---------

getcluster <- function(sc, res) {
	sc <- FindClusters(sc,
			    resolution = res)

	res_cl <- colnames(sc@meta.data)[grep("_snn_res.", colnames(sc@meta.data))]
	print(res_cl)

	g <- DimPlot(sc,
		     reduction = "umap",
		     pt.size = 0.2,
		     group.by = res_cl) +
		ggtitle(res)

	clstr <- select(sc@meta.data, matches(paste0("_snn_res.", res)))
	print(head(clstr, 2))
	rownames(clstr) <- rownames(sc@meta.data)

	return(list('g' = g,
		    'cluster' = clstr))
}
# ------------------------------
sc_file <- arguments$scobject
res <- as.numeric(unlist(strsplit(arguments$res, ",")))

if (is.null(arguments$task)) {
	reduction <- 'pca'
} else if(arguments$task == 'integrate-harmony') {
	reduction <- 'harmony'
} else {
	reduction <- 'pca'
}

# -removed on 25th August 2023, because multicore for large datasets ~120k cells fails
# Memory limits
# options(future.globals.maxSize=91943040000000)
# # ---------------------------------------------------
# # Setting multicore env for some Seurat applications
# # ---------------------------------------------------
# ncores <- as.numeric(arguments$ncores)
# plan("multicore", workers = ncores)
# plan()

# --- Read cell_data_set object with merged samples
samplename <- gsub(pattern="_\\w+.rds",replacement="",sc_file)
sc <- readRDS(sc_file )
assay <- DefaultAssay(sc)
print(sc)
npcs <- ncol(sc@reductions[[reduction]])

# Remove legacy clustering metadata
new_meta <- sc@meta.data[!grepl("_snn_res", colnames(sc@meta.data))]
rownames(new_meta) <- rownames(sc@meta.data)
sc@meta.data <- new_meta[colnames(sc), ]

# --- Rounds of clustering
sc <- FindNeighbors(sc,
    assay = assay,
    reduction=reduction,
    dims = seq(npcs - 20),
    verbose = TRUE)

# -removed on 25th August 2023, because mclapply for large datasets ~120k cells fails
# cls <- mclapply(res, function(r) getcluster(sc, r), mc.cores = ncores)
cls <- lapply(res, function(r) getcluster(sc, r))

# --- Write QC plots
title <- ggdraw() +
    draw_label(
    samplename,     
    fontface = 'bold',
    hjust = 0.5
    ) +
    theme(
    plot.margin = margin(0, 0, 0, 0)
    )
sc_file <- paste0(samplename, '_clustering.rds')
rfile <- paste0(samplename, '_clustering.pdf')


cl_gplots <- lapply(cls, function(cl) cl[['g']])

pdf(rfile, height = 11.69, width = 8.27)
plot(
    plot_grid(title, plot_grid(plotlist = cl_gplots, ncol = 2),
        ncol = 1, rel_heights = c(0.02, 1))
    )
dev.off()

meta_cls <- Reduce(cbind, lapply(cls, function(cl) cl[['cluster']]))

new_meta <- cbind(sc@meta.data, meta_cls[rownames(sc@meta.data),,drop=F])
rownames(new_meta) <- rownames(sc@meta.data)

sc@meta.data <- new_meta[colnames(sc),]
saveRDS(sc, file = sc_file)
message(paste("The cell_data_set object (Seurat3)",
    sc_file, "has been created."))

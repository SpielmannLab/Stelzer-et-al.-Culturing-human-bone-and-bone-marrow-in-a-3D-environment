#
# Single  sample analysis following the Seurat approach
#

" Analyse single sample following the Seurat approach

Usage: dimensionality_reduction.R --scobject=<file> --npcs=<numeric_value> --ncores=<value> [--task=<value>]

Options:
  -h --help Show this screen.
  --scobject=<file>      Seurat RDS object.
  --npcs=<numeric_value> Number of principal components to use.
  --task=<value>         Type of Integration
  --ncores=<value>       Number of threads to use.
"-> doc

library(docopt)
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
suppressMessages(library(future))
# ------------------------------

#---------------------------
#---pca and umap
#---------------------------
dim_reduction <- function(sc,assay,npcs, reduction=reduction){
    sc <- RunPCA(sc,
                 features = VariableFeatures(object = sc),
                 assay = assay,
                 npcs = npcs+20,
                 approx=FALSE,
                 seed.use = 444)
    sc <- RunUMAP(sc, dims = seq(npcs), reduction=reduction)
    return(sc)
}

#---------------------------
#------pca visualization
#----------------------------
pca_sc <- function(sc, assay,npcs = npcs, reduction=reduction){
    DefaultAssay(sc) <- assay

    g_pcafeat <- DimHeatmap(sc,
                            dims = ceiling(seq(npcs)/4),
                            cells = 500,
                            balanced = TRUE,
                            assays = assay,
                            fast = FALSE,
                            reduction = reduction )


    g_elbow <- ElbowPlot(sc,
                         ndims = npcs+20,
                         reduction = reduction)


    return(list('elbow'=g_elbow,
                'pcafeat'=g_pcafeat))
}
#-----------------------
#--------
#----------------------
vis_meta <- function(sc) {

    metavar <- grep("RAW", colnames(sc@meta.data), ignore.case=TRUE, value=TRUE, invert=TRUE)
    metavar <- metavar[!metavar %in% c('cell_log10_depth_counts','metadata_counts','metadata_counts','rank','UMAP_scrublet_1','UMAP_scrublet_2','doublet','barcode','log10_count_p_cell','mt_counts','rb_counts','S.Score','G2M.Score')]

    metadata <- lapply(metavar, function(colgroup) {
                           if(class(sc@meta.data[[colgroup]]) %in% c("integer", "numeric")) {
                               if(length(which(is.na(sc@meta.data[[colgroup]]))) == nrow(sc@meta.data)) {
                                   sc@meta.data[[colgroup]][which(is.na(sc@meta.data[[colgroup]]))] <- 0
                               }
                               sc@meta.data[[colgroup]] <- as.numeric(sc@meta.data[[colgroup]])
                               print("Num:")
                               print(colgroup)
                               print(class(sc@meta.data[[colgroup]]))

                               FeaturePlot(sc,
                                           reduction = "umap",
                                           pt.size = 0.2,
                                           features=colgroup) +
                       #ggplot2::theme(legend.position = "none") +
                       ggtitle(colgroup) +
                       scale_colour_gradient(low = "grey", high = "blue")

                           } else if(class(sc@meta.data[[colgroup]]) %in% c("character", "factor")) {

                               if(length(which(is.na(sc@meta.data[[colgroup]]))) == nrow(sc@meta.data)) {

                                   sc@meta.data[[colgroup]][which(is.na(sc@meta.data[[colgroup]]))] <- "missing"
                               }

                               sc@meta.data[[colgroup]] <- as.factor(sc@meta.data[[colgroup]])

                               print("Fac:")
                               print(colgroup)
                               print(class(sc@meta.data[[colgroup]]))

                               DimPlot(sc,
                                       reduction = "umap",
                                       pt.size = 0.2,
                                       group.by=colgroup) +
                       #ggplot2::theme(legend.position = "none") +
                       ggtitle(colgroup)
                           }
                         })

    return(metadata)
}


# -----------------
# -----Parameters
# -----------------
sc_file <- arguments$scobject
npcs <- as.numeric(arguments$npcs)

if(is.null(arguments$task)){
	reduction <- 'pca'
} else if(arguments$task=='integrate-harmony') {
	reduction <- 'harmony'
} else {
	reduction <- 'pca'
}

# setting parameters of the future package for memory?!
l <- 25000*1024^2
options(future.globals.maxSize = l)

# ---------------------------------------------------
# Setting multicore env for some Seurat applications
# ---------------------------------------------------
ncores <- as.numeric(arguments$ncores)
plan("multicore", workers = ncores)
plan()

# ------------
# Reading merge sc file
# ------------
samplename <- gsub(pattern="_\\w+.rds",replacement="",sc_file)
sc <- readRDS(sc_file)

if(class(sc) != "Seurat") {
    stop(paste("Sorry", print(class(sc)), "objects are not supported!",
	 "Try Seurat object instead!"))
}
assay <- DefaultAssay(sc)
print(npcs)
#Dimensional reduction

sc <- dim_reduction(sc, assay ,npcs = npcs, reduction=reduction)
dim_reduc <- pca_sc(sc, assay ,npcs = npcs, reduction=reduction)

gg_file <- paste0(samplename, "_pca.pdf")
pdf(gg_file,
	width = 500*2/72, height = 500*2/72)
	plot(
	 plot_grid(plotlist = c(dim_reduc), ncol = 2)
	)
	dev.off()

metadata <- vis_meta(sc)
nc <- ceiling(sqrt(length(metadata)))


gg_file <- paste0(samplename, "_cell_metadata.pdf")

pdf(gg_file,
	width = 500*nc/72, height = 500*nc/72)
	plot(
	 plot_grid(plotlist = metadata, ncol = nc)
	)
	dev.off()
sc_file <- paste0(samplename, '_dim_reduction.rds')
saveRDS(sc, file = sc_file)

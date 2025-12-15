" Calculate the differentially expressed genes for a given Seurat object

Usage: de_genes.R --sc_file=<value> --test.use=<value> --min.cells.group=<value> --min.pct=<value> --logfc.threshold=<value> --resolution=<value> [--features=<value>] [--no.of.pages=<value>] --ncores=<value>  

Options:
  -h --help               	Show this screen.
  --version              	00.99.01
  --sc_file=<file>        	*.rds file containing the Seurat Object
  --test.use=<value>	        Which test to use. (roc, wilcoxon, or negbinom)
  --min.cells.group=<value>     Min. number of cell in clusters
  --min.pct=<value> 		Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed.
  --logfc.threshold=<value> 	Min. log fold change among cell clusters of a given gene to be considered.
  --resolution=<value>    	Clustering resolution of interest	
  --features=<value>            Genes for feature plot
  --no.of.pages=<value>      number of pages to plot the feature plot
  --ncores=<value>        	Number of cores to use

"-> doc

# --- Load all the libraries
suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
suppressMessages(library(future))
suppressMessages(library(RColorBrewer))
# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters
sc_file <- arguments$sc_file
test.use <- arguments$test.use
min.pct <- as.numeric(arguments$min.pct)
min.cells.group <- as.numeric(arguments$min.cells.group)
logfc.threshold <- as.numeric(arguments$logfc.threshold)
resolution <- as.numeric(arguments$resolution)
if(is.null(arguments$no.of.pages) || arguments$no.of.pages=='NULL' || arguments$no.of.pages =='null') {
	n <- 0
} else {	
	n <- as.numeric(arguments$no.of.pages)
}
if(is.null(arguments$features) || arguments$features=='NULL' || arguments$features =='null') {
	features <- ''
	nf <- 10
} else {
	features <- unlist(strsplit(arguments$features, ","))
	nf <- ceiling(sqrt(length(features)))
}
ncores <- as.numeric(arguments$ncores)
nf <- ceiling(sqrt(length(features)))

message("sc_file: ",sc_file)
message("min.pct: ",min.pct)
message("min.cells.group: ",min.cells.group)
message("logfc.threshold: ",logfc.threshold)
message("resolution: ",resolution)
message("features to plot:", features)
message("ncores: ",ncores)
message("test.use: ", test.use)

# --- Setting multicore env for some Seurat applications
plan("multicore", workers = ncores)
plan()

# setting parameters of the future package for memory?!
l <- 25000*1024^2
options(future.globals.maxSize = l)

# --- Run
samplename <- gsub(pattern="_\\w+.rds", replacement="", sc_file) # Get sample name from filename

# Read the Seurat object from file
sc <- readRDS(sc_file)
message("Read the following data:")
print(sc)

# Set idents to the cluster identity at the required clustering resolution
meta_col <- meta_name <-  grep(pattern=paste0("snn_res.",resolution),names(sc@meta.data), value=TRUE)
Idents(sc) <- sc@meta.data[meta_col]
print(table(Idents(sc)))
nc <- dim(table(Idents(sc)))
# Visualize the clusters as a UMAP plot
filename <- paste0(samplename, "_clusters_umap.pdf") 
pdf(filename,width=30,height=30)
	DimPlot(sc, reduction = "umap", pt.size = 0.01) + 
		theme(legend.position = "bottom") + 
		ggtitle(paste0(samplename, " @ resolution:", resolution))
dev.off()

# Find marker genes
markers <- FindAllMarkers(sc,
			  logfc.threshold = logfc.threshold,
			  min.pct = min.pct,
			  test.use = test.use,
			  min.cells.group=min.cells.group)

# Save the de genes
filename <- paste0(samplename,"_markers.tsv")
write.table(markers, file=filename, sep = "\t", row.names=FALSE)

# Chose top 10 markers for the heatmap
top_10_markers <- markers %>% group_by(cluster) %>% slice_head(n=10)


# Marker expression visualization (heatmap) [AUC]
filename <- paste0(samplename,"_heatmap.pdf")
pdf(filename, width = 200*nc/70, height = 200*nc/70)
DoHeatmap(
  sc,
  features = top_10_markers$gene,
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  slot = "scale.data",
  label = TRUE,
  raster = FALSE,
  draw.lines = TRUE,
  combine = TRUE) + 

scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
dev.off()

marker <- split(top_10_markers$gene, factor(top_10_markers$cluster, levels=unique(top_10_markers$cluster)))
print(marker)
#Plot the list of genes in UMAP space

	pdf(paste0(samplename, "_featureplot_markergenes.pdf"), width= 500*nf/70, height=500*nf/70)
	for (i in seq(1:length(marker))){
	print(FeaturePlot(
	sc,
	features=marker[[i]],order=TRUE))
}

dev.off()
print(features)
#Plot the list of genes in UMAP space
if(features !=''){
if(length(features) >0){
	print('Entering')
	features <- unique(features)
	chunk <- function(features, n) split(features, seq(1:n))
	genes <- chunk(features,n)
	pdf(paste0(samplename, "_featureplot.pdf"), width= 500*nf/70, height=500*nf/70)
	for (i in seq(1:length(genes))){
	print(FeaturePlot(
	sc,
	features=genes[[i]],order=TRUE))	    
}

dev.off()
}
}

if(features !=''){
if(length(features) > 0){
	features <- unique(features)
	chunk <- function(features, n) split(features, seq(1:n))
	genes <- chunk(features,n)
	pdf(paste0(samplename, "_dotplot.pdf"), width= 500*nf/70, height=500*nf/70)
	for (i in seq(1:length(genes))){
	print(DotPlot(sc, features = features, cols="Spectral")  + RotatedAxis())
	}
dev.off()
}
}

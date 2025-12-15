" Calculate the differentially expressed genes for a given Seurat object

Usage: de_genes_group.R --sc_file=<value> --test.use=<value> --min.cells.group=<value> --min.pct=<value> --logfc.threshold=<value> --resolution=<value> --group.by.col=<value>  --ncores=<value>

Options:
  -h --help               	Show this screen.
  --version              	00.99.01
  --sc_file=<file>        	*.rds file containing the Seurat Object
  --test.use=<value>	        Which test to use. (roc, wilcoxon, or negbinom)
  --min.cells.group=<value>     Min. number of cell in clusters
  --min.pct=<value> 		Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed.
  --logfc.threshold=<value> 	Min. log fold change among cell clusters of a given gene to be considered.
  --resolution=<value>    	Clustering resolution of interest
  --group.by.col=<value>    Column to use for differential gene expression calculation
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
suppressMessages(library(purrr))
suppressMessages(library(ggrepel))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters
sc_file <- arguments$sc_file
test.use <- arguments$test.use
min.pct <- as.numeric(arguments$min.pct)
min.cells.group <- as.numeric(arguments$min.cells.group)
logfc.threshold <- as.numeric(arguments$logfc.threshold)
resolution <- as.numeric(arguments$resolution)
group.by.col <- arguments$group.by.col
ncores <- as.numeric(arguments$ncores)


message("sc_file: ",sc_file)
message("min.pct: ",min.pct)
message("min.cells.group: ",min.cells.group)
message("logfc.threshold: ",logfc.threshold)
message("resolution: ",resolution)
message("group.by.col: ",group.by.col)
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

get_markers <- function(cluster){
sc_sub <- subset(x = sc, idents = cluster)
meta_col <- grep(pattern=group.by.col,names(sc_sub@meta.data), value=TRUE)
sc_sub <- UpdateSeuratObject(sc_sub)
Idents(sc_sub) <- meta_col
if(length(table(Idents(sc_sub)))>1){
markers <- FindAllMarkers(sc_sub,
			  logfc.threshold = logfc.threshold,
			  min.pct = min.pct,
			  test.use = test.use,
			  min.cells.group=min.cells.group) %>%
                cbind(cluster_id = cluster, .)
return(markers)
}
}

markers_group <- map_dfr(as.numeric(levels(Idents(sc))), get_markers)
filename <- paste0(samplename,"_markers_",group.by.col,".tsv")
write.table(markers_group, file=filename, sep = "\t", row.names=FALSE)

top_2_markers <- markers_group %>% group_by(cluster, cluster_id) %>% slice_head(n=2)
nf = length(top_2_markers$gene)

chunk <- function(features,n) split(features, seq(1:n))
genes <- chunk(top_2_markers$gene, length(top_2_markers$gene)/4)
pdf(paste0(samplename,"_markers_",group.by.col,"_top2.pdf"))
for (i in seq(1:length(genes))){
print(VlnPlot(sc, features=genes[[i]], split.by=group.by.col,pt.size=0) + theme(legend.position ='bottom'))
}
dev.off()

markers_group$diffexpressed <- "no"
markers_group$diffexpressed[markers_group$avg_log2FC > 1.5 & markers_group$p_val_adj < 0.05] <- "UP"
markers_group$diffexpressed[markers_group$avg_log2FC < -1.5 & markers_group$p_val_adj < 0.05] <- "Down"
markers_group$delabel <- NA
markers_group$delabel[markers_group$diffexpressed != "no"] <- markers_group$gene[markers_group$diffexpressed != "no"]
markers_group$p_val_plot <- markers_group$p_val_adj
markers_group$p_val_plot[markers_group$p_val_adj == 0] <- 1.452161e-300
g <- ggplot(markers_group, aes(x=avg_log2FC, y=-log10(p_val_plot), col=diffexpressed, label=delabel))+geom_point()+theme(legend.position="None",text = element_text(size = 20)) + geom_text_repel()+facet_wrap(cluster_id~cluster, ncol=length(unique(markers_group$cluster))) +ylim(0,350)
ggsave(paste0(samplename,"_markers_",group.by.col,"_volcano.pdf"), g, width=20, height=20,limitsize=FALSE)

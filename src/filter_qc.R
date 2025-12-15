"Filter the Seurat object based on several conditions and output a file for scrublet

Usage: filter_qc.R --sc_file=<value> --mincount_p_gene=<value> --maxcount_p_gene=<value> --mincell_p_gene=<value> --maxcell_p_gene=<value> --mincount_p_cell=<value> --maxcount_p_cell=<value> --mingene_p_cell=<value> --maxgene_p_cell=<value> --maxpct_mt=<value> --maxpct_rb=<value> --rm_mt=<value> --rm_rb=<value> 

Options:
  -h --help			Show this screen.
  --sc_file=<value>		*.rds file containing Seurat Object
  --mincount_p_gene=<value>	Minimum number of counts for a gene in all cells combined
  --maxcount_p_gene=<value>	Maximum number of counts for a gene in all cells combined
  --mincell_p_gene=<value>	Minimum number of cells in which a gene is detected
  --maxcell_p_gene=<value>	Maximum number of cells in which a gene is detected
  --mincount_p_cell=<value>	Minimum number of counts per cell
  --maxcount_p_cell=<value>	Maximum number of counts per cell
  --mingene_p_cell=<value>	Minimum number of genes detected per cell
  --maxgene_p_cell=<value>	Maximum number of genes detected per cell
  --maxpct_mt=<value>		Maximum mitochondrial-count/total-count allowed per cell.
  --maxpct_rb=<value>		Maximum ribosomal-count/total-count allowed per cell.
  --rm_mt=<value>		TRUE for filtering mitochondrial genes out.
  --rm_rb=<value>		TRUE for filtering ribosomal genes out.

"-> doc

# --- Load all the libraries
suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

# --- Read all the arguments passed 
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters
sc_file <- arguments$sc_file
mincount_p_gene <- as.numeric(arguments$mincount_p_gene)
maxcount_p_gene <- as.numeric(arguments$maxcount_p_gene)
mincell_p_gene <- as.numeric(arguments$mincell_p_gene)
maxcell_p_gene <- as.numeric(arguments$maxcell_p_gene)
mincount_p_cell <- as.numeric(arguments$mincount_p_cell)
maxcount_p_cell <- as.numeric(arguments$maxcount_p_cell)
mingene_p_cell <- as.numeric(arguments$mingene_p_cell)
maxgene_p_cell <- as.numeric(arguments$maxgene_p_cell)
maxpct_mt <- as.numeric(arguments$maxpct_mt)
maxpct_rb <- as.numeric(arguments$maxpct_rb)
rm_mt <- as.logical(arguments$rm_mt)
rm_rb <- as.logical(arguments$rm_rb)

message("sc_file: ", sc_file)
message("mincount_p_gene: ",mincount_p_gene)
message("maxcount_p_gene: ",maxcount_p_gene)
message("mincell_p_gene: ",mincell_p_gene)
message("maxcell_p_gene: ",maxcell_p_gene)
message("mincount_p_cell: ",mincount_p_cell)
message("maxcount_p_cell: ",maxcount_p_cell)
message("mingene_p_cell: ",mingene_p_cell)
message("maxgene_p_cell: ",maxgene_p_cell)
message("maxpct_mt: ",maxpct_mt)
message("maxpct_rb: ",maxpct_rb)
message("rm_mt: ",rm_mt)
message("rm_rb: ",rm_rb)

# --- Run filtering
samplename <- gsub(pattern="_\\w+.rds",replacement="",sc_file) # Get the sample name from filename
sc <- readRDS(sc_file) # Read in the Seurat object

genes <- rownames(sc)
cells <- colnames(sc)


# Remove mitochondrial gene is asked
if(rm_mt) {
	mt_genes <- grep(pattern="^mt-",genes,ignore.case=TRUE,value=TRUE)
	
	if(length(mt_genes) == 0) {
		message("No mitochondrial genes detected to filter out")
	} else {
		message("Mitochondrial genes filtered out:")
		print(mt_genes)
		genes <- genes[!(genes %in% mt_genes)]
	}
}

# Remove ribosomal gene is asked
if(rm_rb) {
	rb_genes <- grep(pattern="^rps|^rpl",genes,value=TRUE,ignore.case=TRUE)

	if(length(rb_genes) == 0) {
		message("No ribosomal genes detected to filter out")
	} else {
		message("Ribosomal genes filtered out:")
		print(rb_genes)
		genes <- genes[!(genes %in% rb_genes)]
	}
}

# subset the Seurat object
sc_subset <- subset(sc,features=genes) #now the Seurat object only contains genes of interest

# figure out genes and cells that do not meet other quality criteria
gene_metadata <- sc[["RNA"]][[]]
genes <- gene_metadata %>%
	dplyr::filter(as.numeric(count_p_gene) > mincount_p_gene & as.numeric(count_p_gene) < maxcount_p_gene & as.numeric(cell_p_gene) > mincell_p_gene & as.numeric(cell_p_gene) < maxcell_p_gene) %>%
	rownames()

cell_metadata <- sc@meta.data
cells <- cell_metadata  %>%
	dplyr::filter(as.numeric(nCount_RNA) > mincount_p_cell & as.numeric(nCount_RNA) < maxcount_p_cell & as.numeric(nFeature_RNA) > mingene_p_cell & as.numeric(nFeature_RNA) < maxgene_p_cell & as.numeric(pct_mt_RAW) < maxpct_mt & as.numeric(pct_rb_RAW) < maxpct_rb) %>%
	rownames()	

# again subset the Seurat object
sc_subset <- sc_subset[genes, cells]


# --- Write filter report plots
report <- data.frame('sc' = c('Before', 'After'),
		     'n_genes' = c(nrow(sc), nrow(sc_subset)),
		     'n_cells' = c(ncol(sc), ncol(sc_subset)))

g1 <- ggplot(report, aes(x = sc, y = n_genes)) +
	geom_bar(stat = "identity", alpha = 0.3) +
	geom_text(aes(label = n_genes), vjust = -1) +
	theme_classic() +
	ggtitle("Features")+xlab("Filtering")+ylab("Total number")

g2 <- ggplot(report, aes(x = sc, y = n_cells)) +
	geom_bar(stat = "identity", color = "#a6bddb", fill = "#a6bddb") +
	geom_text(aes(label = n_cells), vjust = -1) +
	theme_classic() +
	ggtitle("Cells")+xlab("Filtering")+ylab("Total number")

title <- ggdraw() + 
  draw_label(
    paste0(samplename, " filter report (", Sys.time(), ")"),
    fontface = 'bold',
    hjust = 0.5
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 0)
  )

filename <- paste0(samplename,'_','filter_report.pdf')

pdf(filename, height = 5.85, width = 8.27)
    plot_grid(title, plot_grid(g1, g2, ncol = 2), ncol = 1, rel_heights = c(0.06, 1))
dev.off()
message("The filter report has been created")


# --- Write the filtered Seuart object into an rds file
filename <- paste0(samplename, "_filtered.rds") 
saveRDS(sc_subset, file=filename)


# -------------------------------------------
# --- Files for duplets ---------------------
# --- Write matrix for duplets identification
# -------------------------------------------

countmatrix <- GetAssayData(sc_subset, assay="RNA", slot="counts")
filename <- paste0(samplename,".mtx")
Matrix::writeMM(countmatrix, file=filename)

message("The filtered Seurat Object and Matrix for Duplet filtering in the directory: filtered/")

# The following data below does not seem necessary:
if(FALSE){
# --- Write genes for duplets identification pipeline
genes <- data.frame(Reduce(cbind, sub_sc@rowRanges@elementMetadata@listData))
write.table(genes, paste0(outputfolder, '/', jobname, '_', 'filtered_genes_metadata.tsv'),
	    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# --- Write barcodes for duplets identification pipeline
barcodes <- colData(sub_sc)
write.table(barcodes, paste0(outputfolder, '/', jobname, '_', 'filtered_barcodes_metadata.tsv'),
	    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

}

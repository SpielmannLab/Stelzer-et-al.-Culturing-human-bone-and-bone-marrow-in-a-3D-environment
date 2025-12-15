"Make a Scrublet report and filter the seurat object from doublet cells

Usage: filter_doublets.R --sc_file=<value> --threshold=<value> 
Options:
  -h --help		Show this screen.
  --sc_file=<value>	*.rds file containing the Seurat Object
  --threshold=<value>	Doublet score threshold to filter out doublets (score >= threshold).

Author: 
  Cesar Prada
e-mail: 
  prada@molgen.mpg.de

"-> doc

# --- Load all the libraries
suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters
sc_file <- arguments$sc_file
threshold <- as.numeric(arguments$threshold)

message("sc_file: ",sc_file)
message("threshold: ", threshold)

# --- Run
samplename <- gsub(pattern="_\\w+.rds", replacement="", sc_file) # Get sample name from filename

# Read output from run_scrublet.py
# Read single_cell_object
sc <- readRDS(sc_file) # Read in the Seurat object
message("\n","The following object has been read in:")
print(sc)

# Read scrublet outputs
# UMAP
umap_file <- paste0(samplename,"_doublet_umap.csv")
umap <- read.csv(umap_file,header=FALSE)
colnames(umap) <- c("UMAP_scrublet_1", "UMAP_scrublet_2")
rownames(umap) <- Cells(sc)
message("\n","UMAP space from scrublet:")
print(head(umap))
print(dim(umap))

# Duplet observed-score - read an make a histogram
score_file <- paste0(samplename,"_doublet_score.csv")
score <- read.csv(score_file,header=FALSE)
names(score) <- "doublet_score"

g_score <- ggplot(score, aes(doublet_score)) + geom_histogram() +
	geom_vline(xintercept = threshold) + theme_classic()

# Duplet simulated-score
score_sim_file <- paste0(samplename,"_doublet_score_sim.csv")
score_sim <- read.csv(score_sim_file,header=FALSE)
names(score_sim) <- "doublet_score_simulated"

g_score_sim <- ggplot(score_sim, aes(doublet_score_simulated)) + geom_histogram() +
	geom_vline(xintercept = threshold) + theme_classic()

# Apply the desired threhold and classify cells to be doublets or not
scrublet <- cbind(umap, score)
scrublet <- mutate(scrublet, doublet = ifelse(doublet_score >= threshold, "yes", "no"))

g_umap_score <- ggplot(scrublet, aes(x = UMAP_scrublet_1, y = UMAP_scrublet_2)) +
	geom_point(size = 0.3, aes(color = doublet_score), alpha = 0.4) +
	scale_colour_gradient(low = "grey", high = "#67001f") +
	theme_classic()

g_umap_doublet <- ggplot(scrublet, aes(x = UMAP_scrublet_1, y = UMAP_scrublet_2)) +
	geom_point(size=0.3,aes(color = doublet)) +
	ggtitle(paste(table(umap$doublet), collapse = "/")) +
	theme_classic()

# --- Write scrublet report
title <- ggdraw() + 
  draw_label(
    paste("Scrublet analysis and doublet identification on the sample: ", samplename),
    fontface = 'bold', hjust = 0.5
  ) +
  theme(plot.margin = margin(0, 0, 0, 0)) # add margin to align the title 

filename <- paste0(samplename, '_scrublet_report.pdf')
pdf(filename, width = 11.69, height = 8.27)
plot(
    plot_grid(title, plot_grid(g_score, g_score_sim, g_umap_score, g_umap_doublet,ncol = 2),
	      ncol = 1, rel_heights = c(0.1, 0.9))
    )
dev.off()

# --- Filter out the doublet cells 
sc <- AddMetaData(sc, metadata=scrublet) #Adding all the scrublet info as metadata to Seurat
Idents(sc) <- sc$doublet
sc <- subset(sc,ident="no")

# --- Write the filtered Seurat object as rds file

filename <- paste0(samplename,"_doubletFiltered.rds") 
saveRDS(sc, file = filename)
message("The doublet filered Seurat object has been saved as: ", filename)

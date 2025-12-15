#
# Quality Control of Seurat Single cell Object
#
#count_p_gene	Counts for a gene in all cells combined
#cell_p_gene	Number of cells in which a gene is detected
#count_p_cell	Number of counts per cell
#gene_p_cell	Number of genes detected per cell
#mt		Mitochondrial genes out
#rb		Ribosomal genes out
#pct_mt		Percentage mitochondrial-count/total-count per cell
#pct_rb		Percentage ribosomal-count/total-count per cell




"Quality Control of Seurat Single cell Object

Usage: qc_sc.R --scobject=<file>  --dataset=<value>

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --scobject=<file>    *.rds file SeuratObject.
  --dataset=<value>    Either RAW or FILTERED.

"-> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
suppressMessages(library(ggwordcloud))

# --------------
# Main function
# --------------

plot_qc <- function(sc, mt_genes, rb_genes,dataset, samplename) {

    # ---------------------------
    # Feature centered histograms
    # ---------------------------
    counts <- GetAssayData(object = sc, slot = "counts")
    logcount <- log10(rowSums(sc))
    gene_metadata <- data.frame('count_p_gene' = rowSums(sc), # Number of mRNA per gene
                           'log10_count_p_gene' = ifelse(is.infinite(logcount),
							 0.001, logcount), #Log Number of mRNA per gene
			   'cell_p_gene' = Matrix::rowSums(counts > 0), # Number of cells expressing the gene
			   'gene_names' = rownames(sc)
			   )

    genecount <- sum(gene_metadata$count_p_gene)
    gene_metadata$frac_p_gene <- gene_metadata$count_p_gene/genecount

    top20feat <- head(rownames(arrange(gene_metadata, -count_p_gene),20))
    top20exp <- counts[top20feat,]
    count_p_gene_p_totalgenecount <- lapply(seq(ncol(top20exp)), function(i) top20exp[, i]/sc@meta.data$nCount_RNA[i]) # fraction of gene per cell/ total count for gene
    names(count_p_gene_p_totalgenecount) <- colnames(top20exp)

    top20 <- do.call(rbind,
	      lapply(names(count_p_gene_p_totalgenecount), function(cell) {
		  data.frame('count_p_gene_p_totalgenecount'=count_p_gene_p_totalgenecount[[cell]],
          'gene_names'=names(count_p_gene_p_totalgenecount[[cell]]),
			     'cell'=cell)
		  })
	     )

    print(head(top20))


    gtop20 <- ggplot(top20, aes(x = gene_names, y=count_p_gene_p_totalgenecount)) +
	    geom_jitter(aes(group = gene_names), alpha=.2, color = "grey", shape=".") +
	    geom_violin(alpha=.3, fill = "grey") +
	    scale_x_discrete(limits=rev(rownames(head(arrange(gene_metadata, -count_p_gene),20)))) +
	    coord_flip() +
	    theme_classic() +
	    theme(axis.text=element_text(size=5)) +
	    ggtitle(paste(dataset, "Top 20 deepest genes"),
		    subtitle = paste0("They account for ~",
				     ceiling(sum(gene_metadata[top20feat,]$frac_p_gene)*100),
				     "% of total-counts"))

     g1 <- ggplot(`gene_metadata`, aes(`count_p_gene`)) +
	geom_histogram(color = "grey", fill = "grey",
		       bins = ceiling(max(gene_metadata$count_p_gene)/100)) +
	geom_vline(xintercept = ceiling(median(gene_metadata$count_p_gene)),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("median ~",
			       ceiling(median(gene_metadata$count_p_gene))),
		 x = median(gene_metadata$count_p_gene),
		 y = max(hist(gene_metadata$count_p_gene,
			      breaks = seq(from = 0,
					   to = max(gene_metadata$count_p_gene),
					   by = max(gene_metadata$count_p_gene)/ceiling(max(gene_metadata$count_p_gene)/100)),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(dataset, "Feature count depth"))

    g1_1 <- ggplot(gene_metadata, aes(count_p_gene)) +
	  geom_histogram(color = "grey", fill = "grey",
			 bins = (ceiling(max(gene_metadata$count_p_gene)/100))) +
	  theme_classic() +
	  xlim(c(-0.5, ceiling(quantile(gene_metadata$count_p_gene, 0.5)))) +
	  ggtitle(paste(dataset, "Feature count depth"),
		  subtitle = "Zoomed to median")

    g2_wc <- ggplot(head(arrange(gene_metadata, -count_p_gene), nrow(gene_metadata)*0.01),
		    aes(label = gene_names,
				size = count_p_gene,
				color = count_p_gene)) +
		geom_text_wordcloud(rm_outside = TRUE) +
		scale_size_area(max_size = 2.5) +
		scale_color_gradient(low = "lightgrey", high = "black") +
		theme_minimal() +
		ggtitle("Top 1% count-deepest genes")


    g2 <- ggplot(gene_metadata, aes(log10_count_p_gene)) +
	geom_histogram(color = "grey", fill = "grey", bins=100) +
	geom_vline(xintercept = round(quantile(gene_metadata$log10_count_p_gene, 0.5), digits = 2),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("median ~",
			       round(median(gene_metadata$log10_count_p_gene), digits=2)),
		 x = median(gene_metadata$log10_count_p_gene),
		 y = max(hist(gene_metadata$log10_count_p_gene,
			      breaks = seq(from = min(gene_metadata$log10_count_p_gene),
					   to = max(gene_metadata$log10_count_p_gene),
					   by = (max(gene_metadata$log10_count_p_gene)-min(gene_metadata$log10_count_p_gene))/100),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(dataset, "log10(feature count depth)"))

    g1 <- ggplot(`gene_metadata`, aes(`count_p_gene`)) +
	geom_histogram(color = "grey", fill = "grey",
		       bins = ceiling(max(gene_metadata$count_p_gene)/200)) +
	geom_vline(xintercept = ceiling(quantile(gene_metadata$count_p_gene, 0.5)),
		   linetype = "dashed", color = "black") +
	geom_vline(xintercept = ceiling(quantile(gene_metadata$count_p_gene, .99)),
		   linetype = "dotted", color = "grey") +
	ggplot2::annotate("text",
		 label = paste("median ~",
			       ceiling(median(gene_metadata$count_p_gene))),
		 x = median(gene_metadata$count_p_gene),
		 y = max(hist(gene_metadata$count_p_gene,
			      breaks = seq(from = 0,
					   to = max(gene_metadata$count_p_gene),
					   by = max(gene_metadata$count_p_gene)/ceiling(max(gene_metadata$count_p_gene)/200)),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(dataset, "Feature count depth"))


    g2_1 <- ggplot(gene_metadata, aes(log10_count_p_gene)) +
	  geom_histogram(color = "grey", fill = "grey", binwidth = 0.01) +
          xlim(c(-.1, quantile(gene_metadata$log10_count_p_gene, 0.5))) +
	  theme_classic() +
	  ggtitle(paste(dataset, "log10(feature count depth)"),
		  subtitle = "Zoomed to median")

    g2_2 <- ggplot(`gene_metadata`, aes(`cell_p_gene`)) +
	geom_histogram(color = "grey", fill = "grey",
		       bins = ceiling(max(gene_metadata$cell_p_gene)/100)) +
	geom_vline(xintercept = ceiling(quantile(gene_metadata$cell_p_gene, 0.5)),
		   linetype = "dashed", color = "black") +
	geom_vline(xintercept = ceiling(quantile(gene_metadata$cell_p_gene, .99)),
		   linetype = "dotted", color = "grey") +
	ggplot2::annotate("text",
		 label = paste("median ~",
			       ceiling(median(gene_metadata$cell_p_gene))),
		 x = median(gene_metadata$cell_p_gene),
		 y = max(hist(gene_metadata$cell_p_gene,
			      breaks = seq(from = 0,
					   to = max(gene_metadata$cell_p_gene),
					   by = max(gene_metadata$cell_p_gene)/ceiling(max(gene_metadata$cell_p_gene)/100)),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(dataset, "Feature cell depth"))

    g2_2_wc <- ggplot(head(arrange(gene_metadata, -cell_p_gene), nrow(gene_metadata)*0.01),
		    aes(label = gene_names,
				size = cell_p_gene,
				color = cell_p_gene)) +
		geom_text_wordcloud(rm_outside = TRUE) +
		scale_size_area(max_size = 2.5) +
		scale_color_gradient(low = "lightgrey", high = "black") +
		theme_minimal() +
		ggtitle("Top 1% cell-deepest gene_metadata")



    # ------------------
    # Cell based density
    # ------------------

    logcount <- log10(colSums(sc))
    cell_metadata <- data.frame('count_p_cell' = sc@meta.data$nCount_RNA,
                       'log10_count_p_cell' = ifelse(is.infinite(logcount),
							  0.1, logcount),
		       'gene_p_cell' = sc@meta.data$nFeature_RNA
		       )

    cell_metadata[['barcode']] <- rownames(cell_metadata)
    print(head(cell_metadata))

    # ---
    # Mitochondrial gene information
    # ---

    if(length(mt_genes) == 0) {
	    mtc <- setNames(rep(0, ncol(sc)), colnames(sc))
    } else {
	    mtc <- colSums(sc[mt_genes, ])
    }
    mt <- data.frame(mt_counts = mtc)
    mt$pct_mt <- mt$mt_counts/cell_metadata$count_p_cell
    mt$barcode <- rownames(mt)
    print(head(mt))
    cell_metadata <- merge(cell_metadata, mt, by = "barcode")

    print(head(cell_metadata))

    # ---
    # Ribosomal
    # ---

    if(length(rb_genes) == 0) {
	    rbc <- setNames(rep(0, ncol(sc)), colnames(sc))
    } else {
	    rbc <- colSums(sc[rb_genes, ])
    }
    rb <- data.frame(rb_counts = rbc)
    rb$pct_rb <- rb$rb_counts/cell_metadata$count_p_cell
    rb$barcode <- rownames(rb)

    print("Rb genes")
    print(head(cell_metadata))
    cell_metadata <- merge(cell_metadata, rb, by = "barcode")
    print(head(cell_metadata))

    cell_metadata %>%
	dplyr::arrange(-count_p_cell) %>%
	mutate(`rank` = seq(nrow(cell_metadata))) -> cell_metadata
    print(head(cell_metadata,2))

    rownames(cell_metadata) <- cell_metadata$barcode

    g3 <- ggplot(cell_metadata, aes(count_p_cell)) +
	geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 1000) +
	geom_vline(xintercept = ceiling(quantile(cell_metadata$count_p_cell, 0.5)),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("median ~",
			       ceiling(median(cell_metadata$count_p_cell))),
		 x = median(cell_metadata$count_p_cell),
		 y = max(hist(cell_metadata$count_p_cell,
			      breaks = seq(from = 0,
					   to = max(cell_metadata$count_p_cell),
					   by = max(cell_metadata$count_p_cell)/1000),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(dataset, "Count depth"))

    g3_1 <- ggplot(cell_metadata, aes(count_p_cell)) +
	  geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 500) +
	  xlim(c(-.5, ceiling(quantile(cell_metadata$count_p_cell, 0.5)))) +
	  theme_classic() +
	  ggtitle(paste(dataset, "Count depth - Zoomed to median"))

    g4 <- ggplot(cell_metadata, aes(log10_count_p_cell)) +
	geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 100) +
	geom_vline(xintercept = quantile(cell_metadata$log10_count_p_cell, 0.5),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("median ~",
			       round(median(cell_metadata$log10_count_p_cell), digits=2)),
		 x = median(cell_metadata$log10_count_p_cell),
		 y = max(hist(cell_metadata$log10_count_p_cell,
			      breaks = seq(from = 0,
					   to = max(cell_metadata$log10_count_p_cell),
					   by = max(cell_metadata$log10_count_p_cell)/100),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(dataset, "log10(count depth)"))

    g4_1 <- ggplot(cell_metadata, aes(log10_count_p_cell)) +
	  geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 50) +
	  xlim(c(-.1, quantile(cell_metadata$log10_count_p_cell, 0.5))) +
	  theme_classic() +
	  ggtitle(paste(dataset, "log10(count depth) - Zoomed to median"))

    g5 <- ggplot(cell_metadata, aes(gene_p_cell)) +
	geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 100) +
	geom_vline(xintercept = quantile(cell_metadata$gene_p_cell, 0.5),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("median ~",
			       ceiling(median(cell_metadata$gene_p_cell))),
		 x = median(cell_metadata$gene_p_cell),
		 y = max(hist(cell_metadata$gene_p_cell,
			      breaks = seq(from = 0,
					   to = max(cell_metadata$gene_p_cell),
					   by = max(cell_metadata$gene_p_cell)/100),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(dataset, "Gene depth"))

    g5_1 <- ggplot(cell_metadata, aes(gene_p_cell)) +
	  geom_histogram(color = "#a6bddb", fill = "#a6bddb", binwidth = 1) +
	  xlim(c(-.1, quantile(cell_metadata$gene_p_cell, 0.5))) +
	  theme_classic() +
	  ggtitle(paste(dataset, "Gene depth - Zoomed to median"))

    g6 <- ggplot(cell_metadata, aes(x = `rank`, y = count_p_cell)) +
	geom_point(aes(color = pct_mt),
		   stat = "identity", alpha = .2, size = 0.5) +
	scale_color_gradient(low = "black",  high = "red") +
	ylim(c(0, max(cell_metadata$count_p_cell))) +
	theme_classic() +
	ggtitle(paste(dataset, "Depth vs. rank"))

    g6_1 <- ggplot(cell_metadata, aes(x = `rank`, y = count_p_cell)) +
	  geom_point(aes(color = pct_rb),
		     stat = "identity", alpha = .2, size = 0.5) +
	  scale_color_gradient(low = "pink",  high = "black") +
	  ylim(c(0, max(cell_metadata$count_p_cell))) +
	  theme_classic() +
	  ggtitle(paste(dataset, "Depth vs. rank"))

    g7 <- ggplot(cell_metadata, aes(x = count_p_cell, y = gene_p_cell)) +
	geom_point(aes(color = pct_mt), stat = "identity",
		   alpha = 0.2, size = 0.5) +
	scale_color_gradient(low = "black",  high = "red") +
	theme_classic() +
	ggtitle(paste(dataset, "Genes-detected vs. depth"),
		subtitle = "% Mitochondrial-counts")

    g8 <- ggplot(cell_metadata, aes(x = count_p_cell, y = gene_p_cell)) +
	geom_point(aes(color = pct_rb, fill = pct_rb),
		   stat = "identity", alpha = 0.2, size = 0.5) +
	scale_color_gradient(high = "black",  low = "pink") +
	scale_fill_gradient(high = "black",  low = "pink") +
	theme_classic() +
	ggtitle(paste(dataset, "Genes-detected vs. depth"),
	       subtitle = "% Ribosomal-counts")

    g9 <- ggplot(cell_metadata, aes(x = count_p_cell, y = pct_mt)) +
	geom_point(color = "red", alpha = .2, size = 0.5) +
	theme_classic() +
	ggtitle(paste(dataset, "% Mitochondrial-counts vs. depth"))


    if(max(cell_metadata$pct_mt) == 0) {
	    g10 <-  ggplot(cell_metadata, aes(pct_mt)) +
		        geom_density(fill = "red", color = "white", alpha = 0.8) +
			geom_vline(xintercept = median(cell_metadata$pct_mt),
				   linetype = "dashed", color = "black") +
			theme_classic() +
			coord_flip() +
			ggtitle(paste(dataset, "Mitochondrial-counts density"))
    } else {
	    g10 <- ggplot(cell_metadata, aes(pct_mt)) +
			geom_density(fill = "red", color = "white", alpha = 0.8) +
			geom_vline(xintercept = median(cell_metadata$pct_mt),
				   linetype = "dashed", color = "black") +
			ggplot2::annotate("text",
					  label = paste("median ~",
							round(median(cell_metadata$pct_mt), digits=2)),
					  	 x = median(cell_metadata$pct_mt),
						 y = max(hist(cell_metadata$pct_mt,
							      breaks = seq(from = 0,
									   to = max(cell_metadata$pct_mt),
									   by = max(cell_metadata$pct_mt)/1000),
							      plot = FALSE)$density),
					  hjust = 1,
					  vjust = 0
					  ) +
			theme_classic() +
			coord_flip() +
			ggtitle(paste(dataset, "Mitochondrial-counts density"))
    }

    g11 <- ggplot(cell_metadata, aes(x = count_p_cell, y = pct_rb)) +
	geom_point(alpha = 0.4, color = "pink", size = 0.5) +
	theme_classic() +
	ggtitle(paste(dataset, "% Ribosomal-counts vs. depth"))


    if(max(cell_metadata$pct_rb) == 0) {

	    g12 <-  ggplot(cell_metadata, aes(pct_rb)) +
		        geom_density(fill = "pink", color = "white", alpha = 0.8) +
			geom_vline(xintercept = median(cell_metadata$pct_rb),
				   linetype = "dashed", color = "black") +
			theme_classic() +
			coord_flip() +
			ggtitle(paste(dataset, "Ribosomal-counts density"))

    } else {

	    g12 <- ggplot(cell_metadata, aes(pct_rb)) +
			geom_density(fill = "pink", color = "white", alpha = 0.8) +
			geom_vline(xintercept = median(cell_metadata$pct_rb),
				   linetype = "dashed", color = "black") +
			ggplot2::annotate("text",
					  label = paste("median ~",
							round(median(cell_metadata$pct_rb), digits=2)),
					  	 x = median(cell_metadata$pct_rb),
						 y = max(hist(cell_metadata$pct_rb,
							      breaks = seq(from = 0,
									   to = max(cell_metadata$pct_rb),
									   by = max(cell_metadata$pct_rb)/1000),
							      plot = FALSE)$density),
					  hjust = 1,
					  vjust = 0
					  ) +
			theme_classic() +
			coord_flip() +
			ggtitle(paste(dataset, "Ribosomal-counts density"))

    }

    # ------------------------------
    # Save QC metrics in sc object
    # ------------------------------
    if(dataset == "RAW") {
      colnames(cell_metadata) <- paste0(colnames(cell_metadata), '_RAW')
      cell_metadata <- cell_metadata %>% rename('nCount_RAW'='count_p_cell_RAW','nFeature_RAW'= 'gene_p_cell_RAW')

    }else{
     cell_metadata <- cell_metadata %>% select(-c(count_p_cell, gene_p_cell))
    }
    sc <- AddMetaData(sc, metadata=cell_metadata)
    sc[['RNA']] <- AddMetaData(sc[['RNA']], metadata=gene_metadata)

    print(head(sc[['RNA']]))

    return(list("ps1" = list(g1, g2_wc, g2, gtop20, g2_2, g2_2_wc, g3, g3_1, g4, g4_1),
		"ps2" = list(g5, g5_1, g6, g7, g9, g10, g6_1, g8, g11, g12),
		"sc" = sc))
}

# ---------------
# --- Parameters
# ---------------
sc_file <- arguments$scobject
dataset <- arguments$dataset
samplename <- gsub(pattern="_\\w+.rds",replacement="",sc_file)
sc <- readRDS(sc_file)
if(class(sc) != "Seurat") {
    stop(paste("Sorry", print(class(sc)), "objects are not supported!",
	       "Try Seurat object instead!"))
}


	# Mitochondrial genes
	mt_genes <- grep("^mt-", rownames(sc), ignore.case=TRUE, value=TRUE)

	# Ribosomal genes
	rb_genes <- grep("^Rps|^Rpl", rownames(sc), ignore.case=TRUE, value=TRUE)



# --- Run
results <- plot_qc(sc, mt_genes, rb_genes, dataset,samplename)

# --- Write QC plots
title <- ggdraw() +
  draw_label(
    samplename,
    fontface = 'bold',
    #x = 0,
    hjust = 0.5
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 0)
  )

qc_file <- paste0(samplename,'_',dataset,'_qc.pdf')

pdf(qc_file, height = 11.69, width = 8.27)
    plot_grid(title, plot_grid(plotlist = results$ps1, ncol = 2),
	      ncol = 1, rel_heights = c(0.02, 1))
    plot_grid(title, plot_grid(plotlist = results$ps2, ncol = 2),
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()

message(paste("The QC report ", qc_file, "has been created."))

# --- Write cell_data_set object with QC metrics
head(results$sc)
sc_file <- paste0(samplename,'_',dataset,".rds")
saveRDS(results$sc, file = sc_file)

message(paste("The Seurat object",
	      sc_file, "has been created."))

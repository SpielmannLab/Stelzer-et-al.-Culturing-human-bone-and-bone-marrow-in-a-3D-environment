"Covert Seurat object to a monocle object. Choose whether or not to keep the UMAP embeddings. Does the normalization, PCA, and (umap embedding, if asked for). The created monocle file will have the filename with <sc_obj> replaced by <cds>.

Usage: analysis_suerat2monocle.R --file_sc_obj=<file> --assay=<value> --keep_embeddings=<value> --npcs=<value>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		The rds file of a seurat object to convert to monocle object
  --assay=<value>         Which assay to use 
  --keep_embeddings=<value>    To keep the seurat UMAP embeddings or not TRUE or FALSE
  --npcs=<value>        How many pcs to use for the PCA analysis (only if keep_embeddings=FALSE)
" -> doc

# Load libraries
suppressMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressMessages(library(dplyr))
suppressMessages(library(docopt))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters: Read in the output of the sc_multi_sample pipline
file_sc_obj <- arguments$file_sc_obj
keep_embeddings <- arguments$keep_embeddings %>%
    as.logical()
assay <- arguments$assay
if ((assay %in% c("null", "NULL")) || is.null(assay)) {
    assay <- "RNA" #The default assay
}
npcs <- arguments$npcs %>%
    as.numeric()

message("file_sc_obj: ", file_sc_obj)
message("assay: ", assay)
message("keep_embeddings: ", keep_embeddings)
message("npcs: ", npcs)

# read file and apply parameters
sc_obj <- readRDS(file_sc_obj)
DefaultAssay(sc_obj) <- assay

# Conversion to monocle3 - hard coding.
cell_metadata <- sc_obj@meta.data
gene_metadata <- sc_obj@assays$RNA@meta.features %>%
    mutate(gene_id = rownames(.), gene_short_name = rownames(.))
expression_matrix <- GetAssayData(sc_obj, assay = assay, slot = "count")

cds <- new_cell_data_set(expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata)

cds <- preprocess_cds(cds, num_dim = npcs) #does normalization, pca etc.

# Whether to create a new UMAP embedding or to use the one from Seurat
if (keep_embeddings) {
    # replace the pca, as the user wants
    reducedDim(cds, type = "PCA") <- Embeddings(sc_obj, reduction = "pca")
    reducedDim(cds, type = "UMAP") <- Embeddings(sc_obj, reduction = "umap")
} else {
    method_to_use <- "PCA" #alternative is "Aligned"
    seed.use <- 111
    cds <- reduce_dimension(cds, preprocess_method = method_to_use)
}

# save the converted monocle3 cds object
filename <- gsub(file_sc_obj, pattern = "*.rds", replacement = "_cds.rds")
saveRDS(cds, file = filename)

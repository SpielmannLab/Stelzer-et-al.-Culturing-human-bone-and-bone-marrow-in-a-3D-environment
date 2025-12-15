"CellRanger count matrix to Seurat Object

Usage: read_10x.R --infolder=<file> --condition=<value> --treatment=<value> --replicate=<value> --batch=<value> --nickname=<value>
Options:
  -h --help			Show this screen.
  --infolder=<file>		Path to directory containg cell ranger barcodes, cells, and matrix
  --condition=<value>		Name of the sample
  --treatment=<value>		Grouping of the sample
  --replicate=<value>		Replicate info of the sample
  --batch=<value>		Batch of the sample
  --nickname=<value>    Nickname of the sample

" -> doc

# --- Load all the libraries
suppressMessages(library(Seurat))
suppressMessages(library(docopt))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args = TRUE)

# --- Parameters
infolder <- arguments$infolder 
condition <- arguments$condition
treatment <- arguments$treatment
replicate <- arguments$replicate
batch <- arguments$batch
nickname <- arguments$nickname

message("infolder: ", infolder)
message("condition: ", condition)
message("treatment: ", treatment)
message("replicate: ", replicate)
message("batch: ", batch)
message("nickname: ", nickname)

# --- Run
sc_data <- Read10X(data.dir = infolder) # Read the data
project <- paste(condition, treatment, replicate, batch, sep = "-")

# If cellplex is used, sc_data is a list of matrices.
# These include a matrix of gex and another for cell-multiplexing-oligo expression.
# Subset only the gene expression matrix.
if (is.list(sc_data)) {
    sc_data <- sc_data[["Gene Expression"]]
}

sc <- CreateSeuratObject(counts = sc_data, project = project) # Create Seurat Object

# Add metadata
sc$condition <- condition
sc$treatment <- treatment
sc$replicate <- replicate
sc$batch <- batch
sc$nickname <- nickname

print(sc) #preview the created Seurat Object

# --- Save the output rds file
filename <- paste0(project, "_10x.rds")

saveRDS(sc, file = filename)
message("The cell ranger count matrix read and seurat object saved as:", filename)
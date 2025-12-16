" Add spliced and unspliced data from the loom file outputted by Velocyto into Seurat object.

Usage: add_velocyto_to_seurat.R --file_sc_obj=<file> --file_velocyto=<file>

Options:
  -h --help Show this screen.
  --file_sc_obj=<file>        .rds file from 10x read (in which no filtering has been applied).
  --file_velocyto=<file>       loom file containing the the Velocyto output - spliced and unspliced counts.

" -> doc

library(docopt)
arguments <- docopt(doc, quoted_args = TRUE)
print(arguments)

# --- Dependencies
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SeuratWrappers))

# ------------------------------
file_sc_obj <- arguments$file_sc_obj
file_velocyto <- arguments$file_velocyto

sc_obj <- readRDS(file_sc_obj)
velocyto <- ReadVelocity(file_velocyto)

### Store assay in a new variable
spliced <- velocyto[["spliced"]]
colnames(spliced) <- colnames(spliced) %>% 
    gsub(pattern = ".+:", replacement = "") %>% 
    gsub(pattern = "x$", replacement = paste0("-1"))

unspliced <- velocyto[["unspliced"]]
colnames(unspliced) <- colnames(unspliced) %>% 
    gsub(pattern = ".+:", replacement = "") %>% 
    gsub(pattern = "x$", replacement = paste0("-1"))

#confirm if rownames and colnames match
matches <- all(c((rownames(spliced) %in% rownames(sc_obj)),
    (rownames(unspliced) %in% rownames(sc_obj))))
if (!matches) stop("The rownames (genes) do not match between the velocyto & seurat object. Tip - check species.")

matches <- all(c((colnames(spliced) %in% colnames(sc_obj)),
    (colnames(unspliced) %in% colnames(sc_obj))))
if (!matches) stop("The colnames (cell barcodes) do not match between the velocyto & seurat object")

sc_obj[["spliced"]] <- CreateAssayObject(counts = spliced)
sc_obj[["unspliced"]] <- CreateAssayObject(counts = unspliced)

print(sc_obj)

filename <- gsub(file_sc_obj, pattern = ".rds", replacement = "_velocyto.rds")
saveRDS(sc_obj, file = filename)
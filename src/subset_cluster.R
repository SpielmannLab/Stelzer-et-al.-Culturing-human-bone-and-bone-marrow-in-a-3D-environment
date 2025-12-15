" Analyse single sample following the Seurat approach

Usage: subset_cluster.R --scobject=<file> --cluster=<value> --resolution=<value> --integrate=<value> [--integrateBy=<value>]

Options:
  -h --help Show this screen.
  --scobject=<file>        Descriptive name for your experiment.
  --cluster=<value>       Main Cluster number to subcluster
  --resolution=<value>     Clustering resolution of interest
  --integrate=<value>     Boolean to integrate
  --integrateBy=<value>    Required if integrate, Split object metadata
" -> doc

library(docopt)
arguments <- docopt(doc, quoted_args = TRUE)
print(arguments)

# --- Dependencies
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))

#--------------Functions
subcluster <- function(sc, cl, samplename, resolution, integrateBy) {
    sc_sub <- subset(x = sc, idents = cl)
    if(integrate == 'true'){
    sc_list <- SplitObject(sc_sub, split.by = integrateBy)
    for (i in 1:length(sc_list)) {
      fname = paste0(samplename, "-", resolution, "-", cl,"-",i , "_subset.rds")
      saveRDS(file=fname, sc_list[[i]])
      }
    }
    else{
    sc_file <- paste0(samplename, "-", resolution, "-", cl, "_subset.rds")
    saveRDS(sc_sub, file = sc_file)
    }
}


# -----------------
# -----Parameters
# -----------------

sc_file <- arguments$scobject
cluster <- as.numeric(arguments$cluster)
resolution <- arguments$resolution
integrate <- arguments$integrate
integrateBy <- arguments$integrateBy

# ------------
# Reading merge sc file
# ------------
samplename <- gsub(pattern = "_\\w+.rds", replacement = "", sc_file)
sc <- readRDS(sc_file)

if (class(sc) != "Seurat") {
    stop(paste("Sorry", print(class(sc)), "objects are not supported!",
    "Try Seurat object instead!"))
}

print(head(sc@meta.data))

# Set idents to the cluster identity at the required clustering resolution
meta_col <- grep(pattern = paste0("snn_res.", resolution), names(sc@meta.data), value = TRUE)
Idents(sc) <- sc@meta.data[, meta_col]

#Subset
subcluster(sc, cluster, samplename, resolution, integrateBy)

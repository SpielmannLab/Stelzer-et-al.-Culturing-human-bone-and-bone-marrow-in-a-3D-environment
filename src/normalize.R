" Normalize the seurat dataset for a single sample using standard or SCT workflow. Call sex, housekeeping genes as well as covars to regress out.

Usage: normalize.R --scobject=<file> [--method=<value>] [--nhvg=<value>] [--covars=<value>] [--genes_regress=<value>] --ncores=<value>

Options:
  -h --help Show this screen.
  --scobject=<file>        Descriptive name for your experiment.
  --method=<value>         Norm. methods. (SCT, or standard = RNA)
  --nhvg=<value>           Number of highly variable genes to consider.
  --covars=<value>         String separated by , with cell metavariables to regress out.
  --genes_regress=<value>  A csv file containing the list of genes to regress out
  --ncores=<value>         Number of threads to use.
"-> doc

library(docopt)
arguments <- docopt(doc, quoted_args = TRUE)
print(arguments)

# --- Dependencies

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
suppressMessages(library(future))
# ------------------------------

#--------------Functions
norm_sc <- function(sc, nhvg, covars, covar_genes) {
    if (method == "standard") {
        DefaultAssay(sc) <- "RNA"

        sc <- NormalizeData(sc)

        sc <- FindVariableFeatures(sc,
            assay = "RNA",
            selection.method = "vst",
            nfeatures = nhvg)

        sc <- CellCycleScoring(sc,
            s.features = cc.genes$s.genes,
            g2m.features = cc.genes$g2m.genes,
            set.ident = FALSE)

        # Call sex based on total gene expression of all x_linked_genes. Use log-normalized counts
        x_linked_genes <- c("Xist", "Tsix")
        x_linked_genes <- lapply(x_linked_genes,
            FUN = function(gene) {grep(rownames(sc), pattern = gene, value=TRUE, ignore.case = TRUE)}
        ) %>% 
            unlist()
        if (length(x_linked_genes) != 0) {
            x_linked_mtx <- FetchData(sc, vars = x_linked_genes, slot = "data")
            x_linked_sum <- Matrix::rowSums(x_linked_mtx)
            sc$x_linked_sum <- x_linked_sum
            sc@meta.data <- sc@meta.data %>%
                mutate(Sex = case_when((x_linked_sum > 0) ~ "Female",
                (x_linked_sum == 0) ~ "Male"))
        } else {
            sc$x_linked_sum <- NA
            sc$Sex <- NA
        }
        

        hk <- data.frame(t(as.matrix(sc[["RNA"]][grep("^Actb$", rownames(sc[['RNA']]),ignore.case = TRUE), rownames(sc@meta.data)])))
        if (ncol(hk) == 0) {
            sc <- AddMetaData(sc, 0, col.name = 'housekeeping_Lognorm')
        } else {
            sc <- AddMetaData(sc, hk, col.name = 'housekeeping_Lognorm')
            print(head(sc@meta.data, 2))
        }

        if (covar_genes != "") {
            cg <- data.frame(t(as.matrix(sc[["RNA"]][grep(covar_genes, rownames(sc[["RNA"]]),ignore.case = TRUE), rownames(sc@meta.data)])))
            sc <- AddMetaData(sc, cg)
            covars <- append(covars, colnames(cg))
        }

        all.genes <- rownames(sc)
        sc <- ScaleData(sc,
            features = all.genes,
            assay = "RNA",
            vars.to.regress = covars)

        return(sc)

    } else if (method == 'SCT') {
        
        # Call sex based on total gene expression of all x_linked_genes. Use log-normalized counts
        x_linked_genes <- c("Xist", "Tsix")
        x_linked_genes <- lapply(x_linked_genes,
            FUN = function(gene) {grep(rownames(sc), pattern = gene, value=TRUE, ignore.case = TRUE)}
        ) %>% 
            unlist()
        if (length(x_linked_genes) != 0) {
            x_linked_mtx <- FetchData(sc, vars = x_linked_genes, slot = "counts")
            x_linked_sum <- Matrix::rowSums(x_linked_mtx)
            sc$x_linked_sum <- x_linked_sum
            sc@meta.data <- sc@meta.data %>%
                mutate(Sex = case_when((x_linked_sum > 0) ~ "Female",
                (x_linked_sum == 0) ~ "Male"))
        } else {
            sc$x_linked_sum <- NA
            sc$Sex <- NA
        }

        hk <- data.frame(t(as.matrix(sc[["RNA"]][grep("^Actb$", rownames(sc[["RNA"]]),ignore.case=TRUE), rownames(sc@meta.data)])))
        if (ncol(hk) == 0) {
            sc <- AddMetaData(sc, 0, col.name = 'housekeeping_SCT')
        } else {
            sc <- AddMetaData(sc, hk, col.name = 'housekeeping_SCT')
        }
        print(head(sc@meta.data, 2))

        if(covar_genes != ""){
            cg <- data.frame(t(as.matrix(sc[["RNA"]][grep(covar_genes, rownames(sc[["RNA"]]),ignore.case=TRUE), rownames(sc@meta.data)])))
            sc <- AddMetaData(sc,cg)
            covars <- append(covars, colnames(cg))
        }
        sc <- SCTransform(sc,
            assay = "RNA",
            vars.to.regress = covars,
            variable.features.n = nhvg)
        DefaultAssay(sc) <- "SCT"

        return(sc)

    } else {
        print('Enter a method for normalization')
    }
}

# -----------------
# -----Parameters
# -----------------

if (is.null(arguments$covars) || arguments$covars=='NULL' || arguments$covars =='null') {
    covars <- NULL
} else {
    covars <- unlist(strsplit(arguments$covars, ","))
}

if (is.null(arguments$method)|| arguments$method=='NULL' || arguments$method =='null') {
    method <- "standard"
} else {
    method <- arguments$method
}

if (is.null(arguments$genes_regress)|| arguments$genes_regress=='NULL' || arguments$genes_regress =='null') {
    covar_genes <- ""
} else {
    covar_genes <- read.delim(file=arguments$genes_regress, sep=',',header=TRUE)
    covar_genes <- as.vector(covar_genes$genes)
    covar_genes <- paste0("^",covar_genes,"$", collapse='|')
}

sc_file <- arguments$scobject
nhvg <- arguments$nhvg

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
samplename <- gsub(pattern = "_\\w+.rds", replacement="", sc_file)
sc <- readRDS(sc_file)

if (class(sc) != "Seurat") {
    stop(paste("Sorry", print(class(sc)), "objects are not supported!",
    "Try Seurat object instead!"))
}

#Normalizing
norm <- norm_sc(sc, nhvg, covars, covar_genes)
head(norm)
sc_file <- paste0(samplename, '_normalized.rds')
saveRDS(norm, file = sc_file)

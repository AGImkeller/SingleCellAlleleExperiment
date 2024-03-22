## ----options, include=FALSE, echo=FALSE---------------------------------------
library(BiocStyle)

## ----message = FALSE----------------------------------------------------------
library(SingleCellAlleleExperiment)
library(scaeData)
library(scran)
library(scater)
library(patchwork)
library(ggplot2)

## -----------------------------------------------------------------------------
example_data_5k <- scaeData::scaeDataGet(dataset = "pbmc_5k")

## -----------------------------------------------------------------------------
example_data_5k

## ----warning=FALSE------------------------------------------------------------
scae <- read_allele_counts(example_data_5k$dir,
                           sample_names = "example_data",
                           filter = "yes",
                           lookup_file = "lookup_table_HLA_only.csv",
                           barcode_file = example_data_5k$barcodes,
                           gene_file = example_data_5k$features,
                           matrix_file = example_data_5k$matrix,
                           filter_threshold = NULL,
                           example_dataset = TRUE,
                           verbose = TRUE)
scae

## -----------------------------------------------------------------------------
scae_nonimmune__allels_subset <- scae[c(rownames(scae_subset(scae, "nonimmune")), rownames(scae_subset(scae, "alleles"))), ]

scae_nonimmune__allels_subset

## -----------------------------------------------------------------------------
scae_nonimmune__immune <- scae[c(rownames(scae_subset(scae, "nonimmune")), rownames(scae_subset(scae, "immune_genes"))), ]

scae_nonimmune__immune

## -----------------------------------------------------------------------------
scae_nonimmune__functional <- scae[c(rownames(scae_subset(scae, "nonimmune")), rownames(scae_subset(scae, "functional_groups"))), ]

scae_nonimmune__functional

## -----------------------------------------------------------------------------
df_ni_a <- modelGeneVar(scae_nonimmune__allels_subset)
top_ni_a <- getTopHVGs(df_ni_a, prop = 0.1)

## -----------------------------------------------------------------------------
df_ni_g <- modelGeneVar(scae_nonimmune__immune)
top_ni_g <- getTopHVGs(df_ni_g, prop = 0.1)

## -----------------------------------------------------------------------------
df_ni_f <- modelGeneVar(scae_nonimmune__functional)
top_ni_f <- getTopHVGs(df_ni_f, prop = 0.1)

## -----------------------------------------------------------------------------

set.seed(18)
scae <- runPCA(scae, ncomponents = 10, subset_row = top_ni_a, exprs_values = "logcounts", name = "PCA_a")
scae <- runPCA(scae, ncomponents = 10, subset_row = top_ni_g, exprs_values = "logcounts", name = "PCA_g")
scae <- runPCA(scae, ncomponents = 10, subset_row = top_ni_f, exprs_values = "logcounts", name = "PCA_f")

## -----------------------------------------------------------------------------
reducedDimNames(scae)

## -----------------------------------------------------------------------------
set.seed(18)
scae <- runTSNE(scae, dimred= "PCA_a",  name = "TSNE_a")
set.seed(18)
scae <- runTSNE(scae, dimred= "PCA_g",  name = "TSNE_g")
set.seed(18)
scae <- runTSNE(scae, dimred= "PCA_f",  name = "TSNE_f")


## -----------------------------------------------------------------------------
reducedDimNames(scae)

head(reducedDims(scae)$TSNE_a)


## ----fig3, fig.height = 4, fig.width = 12, fig.align = "center", warning = FALSE, message=FALSE----
which_tsne <- "TSNE_a"

#EGA_hla_a_alleles
tsne_g_a  <- plotReducedDim(scae, dimred = which_tsne, colour_by = "HLA-A") + ggtitle("HLA-A gene")
tsne_g_a
tsne_g_a1 <- plotReducedDim(scae, dimred = which_tsne, colour_by = "A*02:01:01:01") + ggtitle("Allele A*02:01:01:01")
tsne_g_a1
tsne_g_a2 <- plotReducedDim(scae, dimred = which_tsne, colour_by = "A*29:02:01:01") + ggtitle("Allele A*29:01:01:01")

p2 <- tsne_g_a + tsne_g_a1 + tsne_g_a2

p2

## -----------------------------------------------------------------------------
sessionInfo()


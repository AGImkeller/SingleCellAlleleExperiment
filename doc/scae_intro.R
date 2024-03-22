## ----options, include=FALSE, echo=FALSE---------------------------------------
library(BiocStyle)

## ----message = FALSE----------------------------------------------------------
library(SingleCellAlleleExperiment)
library(scaeData)
library(scater)


## -----------------------------------------------------------------------------
example_data_5k <- scaeData::scaeDataGet(dataset = "pbmc_5k")

## -----------------------------------------------------------------------------
example_data_5k

## ----warning=FALSE------------------------------------------------------------
scae <- read_allele_counts(example_data_5k$dir,
                           sample_names = "example_data",
                           filter = "no",
                           lookup_file = "lookup_table_HLA_only.csv",
                           barcode_file = example_data_5k$barcodes,
                           gene_file = example_data_5k$features,
                           matrix_file = example_data_5k$matrix,
                           filter_threshold = NULL,
                           example_dataset = TRUE,
                           verbose = FALSE)


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

## ----warning=FALSE------------------------------------------------------------
#this is the object used in the further workflow
scae <- read_allele_counts(example_data_5k$dir,
                           sample_names = "example_data",
                           filter = "custom",
                           lookup_file = "lookup_table_HLA_only.csv",
                           barcode_file = example_data_5k$barcodes,
                           gene_file = example_data_5k$features,
                           matrix_file = example_data_5k$matrix,
                           filter_threshold = 282,
                           example_dataset = TRUE)
scae

## -----------------------------------------------------------------------------
rowData(scae)

## -----------------------------------------------------------------------------
colData(scae)

## -----------------------------------------------------------------------------
scae_nonimmune_subset <- scae_subset(scae, "nonimmune")

head(rownames(scae_nonimmune_subset))

## -----------------------------------------------------------------------------
scae_alleles_subset <- scae_subset(scae, "alleles")

head(rownames(scae_alleles_subset))

## -----------------------------------------------------------------------------
scae_immune_genes_subset <- scae_subset(scae, "immune_genes")

head(rownames(scae_immune_genes_subset))

## -----------------------------------------------------------------------------
scae_functional_groups_subset <- scae_subset(scae, "functional_groups")

head(rownames(scae_functional_groups_subset))

## -----------------------------------------------------------------------------
scae_immune_layers_subset <- c(rownames(scae_subset(scae, "alleles")),
                               rownames(scae_subset(scae, "immune_genes")),
                               rownames(scae_subset(scae, "functional_groups")))

scater::plotExpression(scae, scae_immune_layers_subset)

## -----------------------------------------------------------------------------
sessionInfo()


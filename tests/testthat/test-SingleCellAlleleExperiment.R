#initial tests for the object, not completed yet

library(testthat)
library(SingleCellAlleleExperiment)
library(scaeData)
library(Matrix)

#--------------------read in raw data and objects for tests--------------------#
example_data_5k <- scaeData::scaeDataGet(dataset = "pbmc_5k")

barcode_loc <- file.path(example_data_5k$dir, example_data_5k$barcodes)
feature_loc <- file.path(example_data_5k$dir, example_data_5k$features)
matrix_loc  <- file.path(example_data_5k$dir, example_data_5k$matrix)

feature_info <- utils::read.delim(feature_loc, header = FALSE)
cell_names   <- utils::read.csv(barcode_loc, sep = "", header = FALSE)
mat          <- t(Matrix::readMM(matrix_loc))

sce_filter_raw <- SingleCellExperiment(assays = list(counts = mat),
                                       rowData = feature_info,
                                       colData = cell_names)

# perform zero-filtering on colSums in the example dataset as this is the necessary default performed in the constructor
filtered_sce_raw <- sce_filter_raw[, colSums(counts(sce_filter_raw)) > 0]
cell_names <- cell_names[1:length(colData(filtered_sce_raw)$V1), ]
cell_names$V1 <- colData(filtered_sce_raw)$V1

mat_filtered_zero <- counts(filtered_sce_raw)

scae <- read_allele_counts(example_data_5k$dir,
                           sample_names = "example_data_wta",
                           filter = "custom",
                           exp_type = "WTA",
                           lookup_file = "lookup_table_HLA_only.csv",
                           barcode_file = example_data_5k$barcodes,
                           gene_file = example_data_5k$features,
                           matrix_file = example_data_5k$matrix,
                           filter_threshold = 0,
                           example_dataset = TRUE,
                           verbose = TRUE)

#
test_that("rownames and rowData check", {
  #check the names
  expect_equal(feature_info$V1, rownames(rowData(scae[c(rownames(get_nigenes(scae)), rownames(scae_subset_alleles(scae))),])))

  #check if roWData and object rownames are equal
  expect_equal(rownames(scae[c(rownames(get_nigenes(scae)), rownames(scae_subset_alleles(scae))),]), rownames(rowData(scae[c(rownames(get_nigenes(scae)), rownames(scae_subset_alleles(scae))),])))

  #check the dimension
  expect_equal(length(feature_info$V1), length(rownames(scae[c(rownames(get_nigenes(scae)), rownames(scae_subset_alleles(scae))),])))

})

test_that("colnames and colData check", {
  #check the names
  expect_equal(rownames(colData(scae)), cell_names$V1)

  #check if colData and object colnames are equal
  expect_equal(colnames(scae[c(rownames(get_nigenes(scae)), rownames(scae_subset_alleles(scae))),]), rownames(colData(scae[c(rownames(get_nigenes(scae)), rownames(scae_subset_alleles(scae))),])))

  #check the dimension
  expect_equal(length(colnames(scae[c(rownames(get_nigenes(scae)), rownames(scae_subset_alleles(scae))),])), length(cell_names$V1))
})

test_that("assay check", {
  #dim-check
  expect_equal(dim(counts(scae)[1:dim(mat_filtered_zero)[1],]), dim(mat_filtered_zero))

  scae_no_immune_layers <- scae[1:dim(mat_filtered_zero)[1],]

  random_num1 <- sample(1:dim(mat_filtered_zero)[1], 1)
  random_num2 <- sample(1:dim(mat_filtered_zero)[2], 1)

  #check random if random entry is equal
  expect_equal(counts(scae_no_immune_layers)[random_num1, random_num2], mat_filtered_zero[random_num1, random_num2])

  expect_type(scae[random_num1, random_num2], "S4")
})



#test different filter/no-filter modes
test_that("check input-parameter errors", {

  # filter = "custom" but didnt set a threshold in the filter_threshold param
  expect_error(read_allele_counts(example_data_5k$dir,
                                sample_names = "example_data_wta",
                                filter = "custom",
                                exp_type = "WTA",
                                lookup_file = "lookup_table_HLA_only.csv",
                                barcode_file = example_data_5k$barcodes,
                                gene_file = example_data_5k$features,
                                matrix_file = example_data_5k$matrix,
                                filter_threshold = NULL,
                                example_dataset = TRUE,
                                verbose = FALSE),
            regexp = "")

  expect_message(read_allele_counts(example_data_5k$dir,
                                  sample_names = "example_data_wta",
                                  filter = "yes",
                                  exp_type = "WTA",
                                  lookup_file = "lookup_table_HLA_only.csv",
                                  barcode_file = example_data_5k$barcodes,
                                  gene_file = example_data_5k$features,
                                  matrix_file = example_data_5k$matrix,
                                  tag_feature_mtx = "cells_x_genes.genes.txt",
                                  tag_feature_barcodes = "cells_x_genes.barcodes.txt",
                                  filter_threshold = NULL,
                                  example_dataset = TRUE,
                                  verbose = FALSE),
                 regexp = "Filtering performed based on the inflection point at: 282 UMI counts.")


  expect_message(read_allele_counts(example_data_5k$dir,
                                  sample_names = "example_data_wta",
                                  filter = "no",
                                  exp_type = "WTA",
                                  lookup_file = "lookup_table_HLA_only.csv",
                                  barcode_file = example_data_5k$barcodes,
                                  gene_file = example_data_5k$features,
                                  matrix_file = example_data_5k$matrix,
                                  tag_feature_mtx = "cells_x_genes.genes.txt",
                                  tag_feature_barcodes = "cells_x_genes.barcodes.txt",
                                  filter_threshold = NULL,
                                  example_dataset = TRUE,
                                  verbose = FALSE),
                 regexp = "Suggested threshold based on inflection point is at: 282 UMI counts.")

})


test_that("check rowData extension for WTA and Amplicon", {

 expect_equal(colnames(rowData(scae))[1], "Ensembl_ID")

 expect_equal(colnames(rowData(scae))[2], "Symbol")

 expect_equal(colnames(rowData(scae))[3], "NI_I")

 expect_equal(colnames(rowData(scae))[4], "Quant_type")

})


test_that("sample_names and samples_dir", {

  scae_no_sample <- read_allele_counts(example_data_5k$dir,
                               filter = "custom",
                               exp_type = "WTA",
                               lookup_file = "lookup_table_HLA_only.csv",
                               barcode_file = example_data_5k$barcodes,
                               gene_file = example_data_5k$features,
                               matrix_file = example_data_5k$matrix,
                               filter_threshold = 0,
                               example_dataset = TRUE,
                               verbose = TRUE)

  expect_equal(colData(scae_no_sample)$Sample[1], example_data_5k$dir)
})

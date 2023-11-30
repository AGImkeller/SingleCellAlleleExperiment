library(SingleCellAlleleExperiment)
library(testthat)

dir_path <- system.file("extdata", package = "SingleCellAlleleExperiment")

scae <- readAlleleCounts(dir_path,
                         sample_names = "example_data_wta",
                         filter = "custom",
                         symbols = "orgdb",
                         exp_type = "WTA",
                         lookup_file = "lookup_table_HLA_only.csv",
                         barcode_file = "cells_x_genes.barcodes.txt",
                         gene_file = "cells_x_genes.genes.txt",
                         matrix_file = "cells_x_genes.mtx",
                         tag_feature_mtx = "cells_x_genes.genes.txt",
                         tag_feature_barcodes = "cells_x_genes.barcodes.txt",
                         filter_threshold = 0,
                         verbose = TRUE)


# non_immune genes layer
test_that("non-immune genes getter", {

  expect_equal(get_nigenes(scae), scae[grepl("ENSG", rownames(counts(scae)), fixed = TRUE),])

})


# alleles layer
test_that("alleles getter", {

  expect_equal(get_alleles(scae), scae[grepl("*", rownames(counts(scae)), fixed = TRUE),])

})


# immune gene layer
test_that("immune genes getter", {

  expect_equal(get_agenes(scae), scae[grepl("HLA-", rownames(counts(scae)), fixed = TRUE),])

})


# functional class layer
test_that("functional class getter", {

  expect_equal(get_func(scae), scae[grepl("class", rownames(counts(scae)), fixed = TRUE),])

})


# unknown alleles
test_that("unknown alleles getter", {

  grep_con <- grepl("Unkwn", rownames(counts(scae)), fixed = TRUE)
  expect_equal(get_unknown(scae), scae[grep_con,])

})



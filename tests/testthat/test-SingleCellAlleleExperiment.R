#initial tests for the object, not completed yet

library(SingleCellAlleleExperiment)
library(testthat)

dir_path <- system.file("extdata", package = "SingleCellAlleleExperiment")

barcode_loc <- system.file("extdata", "cells_x_genes.barcodes.txt", package = "SingleCellAlleleExperiment")
feature_loc <- system.file("extdata", "cells_x_genes.genes.txt", package = "SingleCellAlleleExperiment")
matrix_loc <- system.file("extdata", "cells_x_genes.mtx", package = "SingleCellAlleleExperiment")

feature_info <- utils::read.delim(feature_loc, header = FALSE)
cell_names   <- utils::read.csv(barcode_loc, sep = "", header = FALSE)
mat          <- t(Matrix::readMM(matrix_loc))

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

test_that("rownames and rowData check", {
  #check the names
  expect_equal(feature_info$V1, rownames(rowData(scae[c(rownames(get_nigenes(scae)), rownames(get_alleles(scae))),])))

  #check if roWData and object rownames are equal
  expect_equal(rownames(scae[c(rownames(get_nigenes(scae)), rownames(get_alleles(scae))),]), rownames(rowData(scae[c(rownames(get_nigenes(scae)), rownames(get_alleles(scae))),])))

  #check the dimension
  expect_equal(length(feature_info$V1), length(rownames(scae[c(rownames(get_nigenes(scae)), rownames(get_alleles(scae))),])))

})

test_that("colnames and colData check", {
  #check the names
  expect_equal(cell_names$V1, rownames(colData(scae)))

  #check if roWData and object rownames are equal
  expect_equal(colnames(scae[c(rownames(get_nigenes(scae)), rownames(get_alleles(scae))),]), rownames(colData(scae[c(rownames(get_nigenes(scae)), rownames(get_alleles(scae))),])))

  #check the dimension
  expect_equal(length(cell_names$V1), length(colnames(scae[c(rownames(get_nigenes(scae)), rownames(get_alleles(scae))),])))
})


#initial tests for the object, not completed yet

library(SingleCellAlleleExperiment)
library(Matrix)
library(testthat)

#--------------------read in raw data and objects for tests--------------------#
dir_path <- system.file("extdata", package = "SingleCellAlleleExperiment")
barcode_loc <- system.file("extdata", "cells_x_genes.barcodes.txt", package = "SingleCellAlleleExperiment")
feature_loc <- system.file("extdata", "cells_x_genes.genes.txt", package = "SingleCellAlleleExperiment")
matrix_loc  <- system.file("extdata", "cells_x_genes.mtx", package = "SingleCellAlleleExperiment")

feature_info <- utils::read.delim(feature_loc, header = FALSE)
cell_names   <- utils::read.csv(barcode_loc, sep = "", header = FALSE)
mat          <- t(Matrix::readMM(matrix_loc))


#read in using `orgdb`
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
  expect_equal(rownames(colData(scae)), cell_names$V1)

  #check if roWData and object rownames are equal
  expect_equal(colnames(scae[c(rownames(get_nigenes(scae)), rownames(get_alleles(scae))),]), rownames(colData(scae[c(rownames(get_nigenes(scae)), rownames(get_alleles(scae))),])))

  #check the dimension
  expect_equal(length(colnames(scae[c(rownames(get_nigenes(scae)), rownames(get_alleles(scae))),])), length(cell_names$V1))
})

test_that("assay check", {
  #dim-check
  expect_equal(dim(counts(scae)[1:dim(mat)[1],]), dim(mat))

  scae_no_immune_layers <- scae[1:dim(mat)[1],]

  random_num1 <- sample(1:dim(mat)[1], 1)
  random_num2 <- sample(1:dim(mat)[2], 1)

  #check random if random entry is equal
  expect_equal(counts(scae_no_immune_layers)[random_num1, random_num2], mat[random_num1, random_num2])

  expect_type(scae[random_num1, random_num2], "S4")
})



#test different filter/no-filter modes

test_that("check input-parameter errors", {

  # filter = "custom" but didnt set a threshold in the filter_threshold param
  expect_error(readAlleleCounts(dir_path,
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
                                filter_threshold = NULL,
                                verbose = FALSE),
            regexp = "")


  # `symbols` parameter does not equal c("biomart", "orgdb"),
  # also left sample_names param empty which is then automatically set to the dir_path param
  expect_error(readAlleleCounts(dir_path,
                                sample_names = "example_data_wta",
                                filter = "yes",
                                symbols = "notorgdb",
                                exp_type = "WTA",
                                lookup_file = "lookup_table_HLA_only.csv",
                                barcode_file = "cells_x_genes.barcodes.txt",
                                gene_file = "cells_x_genes.genes.txt",
                                matrix_file = "cells_x_genes.mtx",
                                tag_feature_mtx = "cells_x_genes.genes.txt",
                                tag_feature_barcodes = "cells_x_genes.barcodes.txt",
                                filter_threshold = NULL,
                                verbose = FALSE),
               regexp = "Invalid value for symbols parameter. Allowed values are 'biomaRt' and 'orgdb'.")

  expect_message(readAlleleCounts(dir_path,
                                  sample_names = "example_data_wta",
                                  filter = "yes",
                                  symbols = "orgdb",
                                  exp_type = "WTA",
                                  lookup_file = "lookup_table_HLA_only.csv",
                                  barcode_file = "cells_x_genes.barcodes.txt",
                                  gene_file = "cells_x_genes.genes.txt",
                                  matrix_file = "cells_x_genes.mtx",
                                  tag_feature_mtx = "cells_x_genes.genes.txt",
                                  tag_feature_barcodes = "cells_x_genes.barcodes.txt",
                                  filter_threshold = NULL,
                                  verbose = FALSE),
                 regexp = "Filtering performed based on the inflection point at: 105 UMI counts.")


  expect_message(readAlleleCounts(dir_path,
                                  sample_names = "example_data_wta",
                                  filter = "no",
                                  symbols = "orgdb",
                                  exp_type = "WTA",
                                  lookup_file = "lookup_table_HLA_only.csv",
                                  barcode_file = "cells_x_genes.barcodes.txt",
                                  gene_file = "cells_x_genes.genes.txt",
                                  matrix_file = "cells_x_genes.mtx",
                                  tag_feature_mtx = "cells_x_genes.genes.txt",
                                  tag_feature_barcodes = "cells_x_genes.barcodes.txt",
                                  filter_threshold = NULL,
                                  verbose = FALSE),
                 regexp = "Suggested threshold based on inflection point is at: 105 UMI counts.")

})


test_that("test for unknown alleles that have correct nomenclature", {
  path_genes  <- file.path(dir_path, "cells_x_genes_unknwn.genes.txt")
  path_matrix <- file.path(dir_path, "cells_x_genes_unknwn.mtx")


  #generate new feature list, with the unknown alleles integrated
  new_entries <- data.frame(V1 = c("Unkwn*01:01:01:01", "Unkwn*02:02:02:02"))
  feature_info_new <- rbind(feature_info, new_entries)

  #generate two count matrix, integrating two more rows (all 0 counts)
  unknown_counts <- matrix(0, nrow = 2, ncol = dim(mat)[2])
  unknown_counts <- as(unknown_counts, "CsparseMatrix")
  new_full_matrix <- rbind(mat, unknown_counts)
  new_full_matrix <- t(new_full_matrix)

  #overwrite raw data and integrate two unknown alleles and update the matrix as well
  write.table(feature_info_new$V1, file = path_genes, sep = "\t", row.names = FALSE, col.names = FALSE)
  writeMM(new_full_matrix, path_matrix)

  expect_message(scae_unknown <- readAlleleCounts(dir_path,
                                  sample_names = "example_data_wta",
                                  filter = "yes",
                                  symbols = "orgdb",
                                  exp_type = "WTA",
                                  lookup_file = "lookup_table_HLA_only.csv",
                                  barcode_file = "cells_x_genes.barcodes.txt",
                                  gene_file = "cells_x_genes_unknwn.genes.txt",
                                  matrix_file = "cells_x_genes_unknwn.mtx",
                                  tag_feature_mtx = "cells_x_genes.genes.txt",
                                  tag_feature_barcodes = "cells_x_genes.barcodes.txt",
                                  filter_threshold = NULL,
                                  verbose = FALSE),
                 regexp = message("Unkwn*01:01:01:01 can't be found in the lookup table
Unkwn*02:02:02:02 can't be found in the lookup table"))

# test new dimension if integration worked
expect_equal(dim(counts(scae_unknown[c(rownames(get_nigenes(scae_unknown)), rownames(get_alleles(scae_unknown))),]))[1],
             (dim(mat)[1] + dim(unknown_counts)[1]))

# delete files with integrated unknown alleles again
unlink(path_genes)
unlink(path_matrix)

})


test_that("test for unknown alleles that have arbitrary identifier (stop execution)", {

  path_genes  <- file.path(dir_path, "cells_x_genes_unknwn_stop.genes.txt")
  path_matrix <- file.path(dir_path, "cells_x_genes_unknwn_stop.mtx")

  #generate new feature list, with the unknown alleles integrated, both dont have a valid nomenclature
  new_entries_stop <- data.frame(V1 = c("Unkwn_allele1", "Unkwn_allele2"))
  feature_info_new_stop <- rbind(feature_info, new_entries_stop)

  #generate two count matrix, integrating two more rows (all 0 counts)
  unknown_counts_stop <- matrix(0, nrow = 2, ncol = dim(mat)[2])
  unknown_counts_stop <- as(unknown_counts_stop, "CsparseMatrix")
  new_full_matrix_stop <- rbind(mat, unknown_counts_stop)
  new_full_matrix_stop <- t(new_full_matrix_stop)

  #overwrite raw data and integrate two unknown alleles and update the matrix as well
  write.table(feature_info_new_stop$V1, file = path_genes, sep = "\t", row.names = FALSE, col.names = FALSE)
  writeMM(new_full_matrix_stop, path_matrix)

  expect_error(readAlleleCounts(dir_path,
                                sample_names = "example_data_wta",
                                filter = "yes",
                                symbols = "orgdb",
                                exp_type = "WTA",
                                lookup_file = "lookup_table_HLA_only.csv",
                                barcode_file = "cells_x_genes.barcodes.txt",
                                gene_file = "cells_x_genes_unknwn_stop.genes.txt",
                                matrix_file = "cells_x_genes_unknwn_stop.mtx",
                                tag_feature_mtx = "cells_x_genes.genes.txt",
                                tag_feature_barcodes = "cells_x_genes.barcodes.txt",
                                filter_threshold = NULL,
                                verbose = TRUE),
               regexp = message("Error in check_unknowns(sce, allele_ids_lookup) :
Allele information contains unknown identifier.
Please check the data and remove rows of the following
allele features identifiers: `Unkwn_allele1 Unkwn_allele2` or use proper nomenclature."))


# delete files with integrated unknown alleles again
unlink(path_genes)
unlink(path_matrix)

})


test_that("check rowData extension with biomart and orgdb", {

 expect_equal(colnames(rowData(scae))[1], "Ensembl_ID")

 expect_equal(colnames(rowData(scae))[2], "Symbol")

 expect_equal(colnames(rowData(scae))[3], "NI_I")

 expect_equal(colnames(rowData(scae))[4], "Quant_type")

})

test_that("check biomart and orgdb", {

  expect_message(readAlleleCounts(dir_path,
                                  sample_names = "example_data_wta",
                                  filter = "yes",
                                  symbols = "biomart",
                                  exp_type = "WTA",
                                  lookup_file = "lookup_table_HLA_only.csv",
                                  barcode_file = "cells_x_genes.barcodes.txt",
                                  gene_file = "cells_x_genes.genes.txt",
                                  matrix_file = "cells_x_genes.mtx",
                                  tag_feature_mtx = "cells_x_genes.genes.txt",
                                  tag_feature_barcodes = "cells_x_genes.barcodes.txt",
                                  filter_threshold = NULL,
                                  verbose = TRUE))

  #regexp = message("Using biomart to retrieve NCBI gene identifiers.")
  expect_message(readAlleleCounts(dir_path,
                                  sample_names = "example_data_wta",
                                  filter = "yes",
                                  symbols = "orgdb",
                                  exp_type = "WTA",
                                  lookup_file = "lookup_table_HLA_only.csv",
                                  barcode_file = "cells_x_genes.barcodes.txt",
                                  gene_file = "cells_x_genes.genes.txt",
                                  matrix_file = "cells_x_genes.mtx",
                                  tag_feature_mtx = "cells_x_genes.genes.txt",
                                  tag_feature_barcodes = "cells_x_genes.barcodes.txt",
                                  filter_threshold = NULL,
                                  verbose = TRUE),
                 regexp = message("Using org.HS to retrieve NCBI gene identifiers."))
})

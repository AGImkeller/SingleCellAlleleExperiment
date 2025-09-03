library(testthat)
library(SingleCellAlleleExperiment)
library(scaeData)
library(Matrix)
#--------------------read in raw data and objects for tests--------------------#
example_data_5k <- scaeData::scaeDataGet(dataset="pbmc_5k")
lookup <- utils::read.csv(system.file("extdata",
                                      "pbmc_5k_lookup_table.csv",
                                      package="scaeData"))

barcode_loc <- file.path(example_data_5k$dir, example_data_5k$barcodes)
feature_loc <- file.path(example_data_5k$dir, example_data_5k$features)
matrix_loc  <- file.path(example_data_5k$dir, example_data_5k$matrix)

feature_info <- utils::read.delim(feature_loc, header=FALSE)
cell_names   <- utils::read.csv(barcode_loc, sep="", header=FALSE)
mat          <- t(Matrix::readMM(matrix_loc))

sce_filter_raw <- SingleCellExperiment(assays = list(counts=mat),
                                       rowData = feature_info,
                                       colData = cell_names)

## perform zero-filtering on colSums in the example dataset
## this is the necessary default performed in the constructor
filtered_sce_raw <- sce_filter_raw[, colSums(counts(sce_filter_raw)) > 0]
cell_names <- cell_names[1:length(colData(filtered_sce_raw)$V1), ]
cell_names <- as.data.frame(cell_names)
cell_names$V1 <- colData(filtered_sce_raw)$V1


mat_filtered_zero <- counts(filtered_sce_raw)

scae <- read_allele_counts(example_data_5k$dir,
                           filter_mode="custom",
                           lookup_file=lookup,
                           barcode_file=example_data_5k$barcodes,
                           gene_file=example_data_5k$features,
                           matrix_file=example_data_5k$matrix,
                           filter_threshold=0,
                           log=FALSE,
                           gene_symbols=TRUE,
                           verbose=TRUE)
show(scae)

test_that("validty test", {

  expect_true(validObject(scae))

  expect_error(rowData(scae)$NI_I <- "random_value", regexp="")
  expect_error(rowData(scae)$NI_I <- NULL, regexp="")
  expect_error(rowData(scae)$Quant_type <- "random_value", regexp="")
  expect_error(rowData(scae) <- "test", regexp="")

  old_scae <- scae
  rowData(scae) <- NULL
  retain_cols <- c("NI_I", "Quant_type")

  expect_equal(rowData(old_scae)[,retain_cols], rowData(scae))

  ## change 'NI_I' column
  expect_error(colnames(rowData(scae))[2] <- "arbitrary_name")
  ## change 'Quant_type' column
  expect_error(colnames(rowData(scae))[3] <- "arbitrary_name")
})


test_that("rownames and rowData check", {

  rn_ni_genes <- rownames(get_nigenes(scae))
  rn_alleles <- rownames(scae_subset_alleles(scae))
  ## check the names
  expect_equal(feature_info$V1,
               rownames(rowData(scae[c(rn_ni_genes, rn_alleles),])))

  ## check if roWData and object rownames are equal
  expect_equal(rownames(scae[c(rn_ni_genes, rn_alleles),]),
               rownames(rowData(scae[c(rn_ni_genes, rn_alleles),])))

  ## check the dimension
  expect_equal(length(feature_info$V1),
               length(rownames(scae[c(rn_ni_genes, rn_alleles),])))

})

test_that("colnames and colData check", {
  ## check the names
  rn_ni_genes <- rownames(get_nigenes(scae))
  rn_alleles <- rownames(scae_subset_alleles(scae))

  expect_equal(rownames(colData(scae)), cell_names$V1)

  ## check if colData and object colnames are equal
  expect_equal(colnames(scae[c(rn_ni_genes, rn_alleles),]),
               rownames(colData(scae[c(rn_ni_genes, rn_alleles),])))

  ## check the dimension
  expect_equal(length(colnames(scae[c(rn_ni_genes, rn_alleles),])),
               length(cell_names$V1))
})

test_that("assay check", {
  ## dim-check
  expect_equal(dim(counts(scae)[1:dim(mat_filtered_zero)[1],]),
               dim(mat_filtered_zero))

  scae_no_immune_layers <- scae[1:dim(mat_filtered_zero)[1],]

  random_num1 <- sample(1:dim(mat_filtered_zero)[1], 1)
  random_num2 <- sample(1:dim(mat_filtered_zero)[2], 1)

  ## check if random entry is equal
  expect_equal(counts(scae_no_immune_layers)[random_num1, random_num2],
               mat_filtered_zero[random_num1, random_num2])

  expect_type(scae, "S4")
})

## Test different filter/no-filter modes
test_that("check input-parameter errors", {

  ## Testing error message of the filter_mode="custom" if threshold is not set
  expect_error(read_allele_counts(example_data_5k$dir,
                                  sample_names="example_data_wta",
                                  filter_mode="custom",
                                  lookup_file=lookup,
                                  barcode_file=example_data_5k$barcodes,
                                  gene_file=example_data_5k$features,
                                  matrix_file=example_data_5k$matrix,
                                  filter_threshold=NULL,
                                  verbose=FALSE),
               regexp = "")

#
  ## Testing for the correct output message of the filter_mode="yes",
  ## also testing for error message if package is not installed
  if (!requireNamespace("DropletUtils", quietly=TRUE)) {
    expect_error(read_allele_counts(example_data_5k$dir,
                                      sample_names="example_data_wta",
                                      filter_mode="yes",
                                      lookup_file=lookup,
                                      barcode_file=example_data_5k$barcodes,
                                      gene_file=example_data_5k$features,
                                      matrix_file=example_data_5k$matrix,
                                      filter_threshold=NULL,
                                      verbose=FALSE),
                 regexp = "")
  }else {
    expect_message(read_allele_counts(example_data_5k$dir,
                                      sample_names="example_data_wta",
                                      filter_mode="yes",
                                      lookup_file=lookup,
                                      barcode_file=example_data_5k$barcodes,
                                      gene_file=example_data_5k$features,
                                      matrix_file=example_data_5k$matrix,
                                      filter_threshold=NULL,
                                      verbose=FALSE),
    regexp = "Filtering performed based on the inflection point at: [0-9]+(\\.[0-9]+)? UMI counts.")
  }


  ## Testing error/output message for `log=TRUE`
  ## also testing for error message if package is not installed
  if (!requireNamespace("scuttle", quietly=TRUE)) {
    expect_error(read_allele_counts(example_data_5k$dir,
                                    sample_names="example_data_wta",
                                    filter_mode="no",
                                    lookup_file=lookup,
                                    barcode_file=example_data_5k$barcodes,
                                    gene_file=example_data_5k$features,
                                    matrix_file=example_data_5k$matrix,
                                    filter_threshold=NULL,
                                    log=TRUE,
                                    verbose=TRUE),
                 regexp="")
  }else {
    expect_message(read_allele_counts(example_data_5k$dir,
                                      sample_names="example_data_wta",
                                      filter_mode="no",
                                      lookup_file=lookup,
                                      barcode_file=example_data_5k$barcodes,
                                      gene_file=example_data_5k$features,
                                      matrix_file=example_data_5k$matrix,
                                      filter_threshold=NULL,
                                      log=TRUE,
                                      verbose=TRUE),
                   regexp="  Generating SCAE object: Generate logcounts assay using Library Factors")
  }


  ## Testing error/output message for `gene_symbols=TRUE`
  ## also testing for error message if package is not installed
  if (!requireNamespace("org.Hs.eg.db", quietly=TRUE)) {
    expect_error(read_allele_counts(example_data_5k$dir,
                                    sample_names="example_data_wta",
                                    filter_mode="no",
                                    lookup_file=lookup,
                                    barcode_file=example_data_5k$barcodes,
                                    gene_file=example_data_5k$features,
                                    matrix_file=example_data_5k$matrix,
                                    filter_threshold=NULL,
                                    gene_symbols=TRUE,
                                    verbose=TRUE),
                 regexp="")
  }else {
    expect_message(read_allele_counts(example_data_5k$dir,
                                      sample_names="example_data_wta",
                                      filter_mode="no",
                                      lookup_file=lookup,
                                      barcode_file=example_data_5k$barcodes,
                                      gene_file=example_data_5k$features,
                                      matrix_file=example_data_5k$matrix,
                                      filter_threshold=NULL,
                                      gene_symbols=TRUE,
                                      verbose=TRUE),
                   regexp="Using org.Hs to retrieve NCBI gene identifiers.")
  }

  ## Testing for the correct output message of the filter_mode="no"
  expect_message(read_allele_counts(example_data_5k$dir,
                                    sample_names="example_data_wta",
                                    filter_mode="no",
                                    lookup_file=lookup,
                                    barcode_file=example_data_5k$barcodes,
                                    gene_file=example_data_5k$features,
                                    matrix_file=example_data_5k$matrix,
                                    filter_threshold=NULL,
                                    verbose=FALSE),
               regexp = "Filtering performed on default value at 0 UMI counts.")

})

test_that("check rowData extension for WTA and Amplicon", {

  expect_equal(colnames(rowData(scae))[1], "Ensembl_ID")

  expect_equal(colnames(rowData(scae))[2], "NI_I")

  expect_equal(colnames(rowData(scae))[3], "Quant_type")

  expect_equal(colnames(rowData(scae))[4], "Symbol")

})

test_that("sample_names and samples_dir", {

  expect_equal(colData(scae)$Sample[1], example_data_5k$dir)
})

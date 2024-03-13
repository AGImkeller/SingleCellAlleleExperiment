library(testthat)
library(SingleCellAlleleExperiment)
library(scaeData)

example_data_5k <- scaeData::scaeDataGet(dataset = "pbmc_5k")

scae <- read_allele_counts(example_data_5k$dir,
                         sample_names = "example_data_wta",
                         filter = "custom",
                         lookup_file = "pbmc_5k_lookup_table.csv",
                         barcode_file = example_data_5k$barcodes,
                         gene_file = example_data_5k$features,
                         matrix_file = example_data_5k$matrix,
                         filter_threshold = 0,
                         example_dataset = TRUE,
                         verbose = TRUE)


# non_immune genes layer
test_that("non-immune genes getter", {

  expect_equal(get_nigenes(scae), scae[grepl("ENSG", rownames(counts(scae)), fixed = TRUE),])

})


# alleles layer
test_that("alleles getter", {

  expect_equal(scae_subset_alleles(scae), scae[grepl("*", rownames(counts(scae)), fixed = TRUE),])

})


# immune gene layer
test_that("immune genes getter", {

  expect_equal(get_agenes(scae), scae[grepl("HLA-", rownames(counts(scae)), fixed = TRUE),])

})


# functional class layer
test_that("functional class getter", {

  expect_equal(scae_subset_functional(scae), scae[grepl("class", rownames(counts(scae)), fixed = TRUE),])

})


test_that("test wrapper getter", {

  #nonimmune
  expect_equal(scae_subset(scae, "nonimmune"), get_nigenes(scae))
  #allele
  expect_equal(scae_subset(scae, "alleles"), scae_subset_alleles(scae))
  #immune
  expect_equal(scae_subset(scae, "immune_genes"), get_agenes(scae))
  #functional
  expect_equal(scae_subset(scae, "functional_groups"), scae_subset_functional(scae))

  expect_message(scae_subset(scae, "wrong_layer"),
                 regexp = "Invalid layer specified, Choose from `nonimmune`, `alleles`, `immune_genes`, `functional_groups`")

})


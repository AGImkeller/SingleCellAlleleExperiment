
 #---------------------------------getters-------------------------------------#

#' Get allele rows
#'
#' @description
#' Getter function returning subsampled SCAE object with all rows containing raw allele information. These rows
#' are identified by "I" in rowData(scae)$NI_I and "A" in rowData(scae)$Quant_type.
#'
#' @param scae A \code{\link{SingleCellAlleleExperiment}} object.
#'
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#'
#' @return A SingleCellAlleleExperiment object.
#'
#' @examples
#'
#' example_data_5k <- scaeData::scaeDataGet(dataset = "pbmc_5k")
#'
#' scae <- read_allele_counts(example_data_5k$dir,
#'                           sample_names = "example_data_wta",
#'                           filter = "custom",
#'                           exp_type = "WTA",
#'                           lookup_file = "lookup_table_HLA_only.csv",
#'                           barcode_file = example_data_5k$barcodes,
#'                           gene_file = example_data_5k$features,
#'                           matrix_file = example_data_5k$matrix,
#'                           filter_threshold = 0,
#'                           example_dataset = TRUE,
#'                           verbose = TRUE)
#'
#'
#'scae
#'
#'rownames(scae_subset_alleles(scae))
#'
#'scae_alleles <- scae_subset_alleles(scae)
#'
#'scae_alleles
#'
#'rownames(scae_alleles)
#'
#' @export
scae_subset_alleles <- function(scae) {
  subset_rows <- stats::complete.cases(rowData(scae)$NI_I, rowData(scae)$Quant_type)
  # alleles of the genes with extended quantification
  alleles <- scae[subset_rows & rowData(scae)$NI_I == "I" & startsWith(rowData(scae)$Quant_type,"A"), ]
  return(alleles)
}


#' Get immune gene rows
#'
#' @description
#' Getter function returning subsampled SCAE object with all rows containing immune gene information. These rows are
#' identfied by "I" in rowData(scae)$NI_I and "G" in rowData(scae)$Quant_type.
#'
#' @param scae A \code{\link{SingleCellAlleleExperiment}} object.
#'
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#'
#' @return A SingleCellAlleleExperiment object.
#'
#' @examples
#'
#' example_data_5k <- scaeData::scaeDataGet(dataset = "pbmc_5k")
#'
#' scae <- read_allele_counts(example_data_5k$dir,
#'                           sample_names = "example_data_wta",
#'                           filter = "custom",
#'                           exp_type = "WTA",
#'                           lookup_file = "lookup_table_HLA_only.csv",
#'                           barcode_file = example_data_5k$barcodes,
#'                           gene_file = example_data_5k$features,
#'                           matrix_file = example_data_5k$matrix,
#'                           filter_threshold = 0,
#'                           example_dataset = TRUE,
#'                           verbose = TRUE)
#'
#'scae
#'
#'rownames(get_agenes(scae))
#'
#'scae_immune_genes <- get_agenes(scae)
#'
#'scae_immune_genes
#'
#'rownames(scae_immune_genes)
#'
#' @export
get_agenes <- function(scae) {
  subset_rows <- stats::complete.cases(rowData(scae)$NI_I, rowData(scae)$Quant_type)
  # genes with extended quantification
  agenes <- scae[subset_rows & rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "G", ]
  return(agenes)
}


#' Get non-immune rows
#'
#' @description
#' Getter function returning subsampled SCAE object with all rows containing non immune gene information. These rows are
#' identified by "NI" in rowData(scae)$NI_I and "G" in rowData(scae)$Quant_type.
#'
#' @param scae A \code{\link{SingleCellAlleleExperiment}} object.
#'
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#'
#' @return A SingleCellAlleleExperiment object.
#'
#' @examples
#'
#' example_data_5k <- scaeData::scaeDataGet(dataset = "pbmc_5k")
#'
#' scae <- read_allele_counts(example_data_5k$dir,
#'                           sample_names = "example_data_wta",
#'                           filter = "custom",
#'                           exp_type = "WTA",
#'                           lookup_file = "lookup_table_HLA_only.csv",
#'                           barcode_file = example_data_5k$barcodes,
#'                           gene_file = example_data_5k$features,
#'                           matrix_file = example_data_5k$matrix,
#'                           filter_threshold = 0,
#'                           example_dataset = TRUE,
#'                           verbose = TRUE)
#'
#'scae
#'
#'rownames(get_nigenes(scae))
#'
#'scae_non_immune_genes <- get_nigenes(scae)
#'
#'scae_non_immune_genes
#'
#'rownames(scae_non_immune_genes)
#'
#' @export
get_nigenes <- function(scae) {
  subset_rows <- stats::complete.cases(rowData(scae)$NI_I, rowData(scae)$Quant_type)
  # classical genes
  agenes <- scae[subset_rows & rowData(scae)$NI_I == "NI" & rowData(scae)$Quant_type == "G", ]
  return(agenes)
}


#' Get functional class rows
#'
#' @description
#' Getter function returning subsampled SCAE object with all rows containing functional class information. These rows are
#' identified by "I" in rowData(scae)$NI_I and "F" in rowData(scae)$Quant_type.
#'
#' @param scae A \code{\link{SingleCellAlleleExperiment}} object.
#'
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#'
#' @return A SingleCellAlleleExperiment object.
#'
#' @examples
#'
#' example_data_5k <- scaeData::scaeDataGet(dataset = "pbmc_5k")
#'
#' scae <- read_allele_counts(example_data_5k$dir,
#'                           sample_names = "example_data_wta",
#'                           filter = "custom",
#'                           exp_type = "WTA",
#'                           lookup_file = "lookup_table_HLA_only.csv",
#'                           barcode_file = example_data_5k$barcodes,
#'                           gene_file = example_data_5k$features,
#'                           matrix_file = example_data_5k$matrix,
#'                           filter_threshold = 0,
#'                           example_dataset = TRUE,
#'                           verbose = TRUE)
#'
#'scae
#'
#'rownames(scae_subset_functional(scae))
#'
#'scae_functional_class <- scae_subset_functional(scae)
#'
#'scae_functional_class
#'
#'rownames(scae_functional_class)
#'
#' @export
scae_subset_functional <- function(scae) {
  subset_rows <- stats::complete.cases(rowData(scae)$NI_I, rowData(scae)$Quant_type)
  # functional groups of the genes with extended quantification
  func <- scae[subset_rows & rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "F", ]
  return(func)
}


#' Get unknown allele rows
#'
#' @description
#' Getter function returning subsampled SCAE object with all rows containing unknown allele information thats present in the data and in the proper nomenclature.
#' These rows are identified by "I" in rowData(scae)$NI_I and "A_unknown" in rowData(scae)$Quant_type.
#'
#' @param scae A \code{\link{SingleCellAlleleExperiment}} object.
#'
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#'
#' @return A SingleCellAlleleExperiment object.
#'
#' @examples
#'
#' example_data_5k <- scaeData::scaeDataGet(dataset = "pbmc_5k")
#'
#' scae <- read_allele_counts(example_data_5k$dir,
#'                           sample_names = "example_data_wta",
#'                           filter = "custom",
#'                           exp_type = "WTA",
#'                           lookup_file = "lookup_table_HLA_only.csv",
#'                           barcode_file = example_data_5k$barcodes,
#'                           gene_file = example_data_5k$features,
#'                           matrix_file = example_data_5k$matrix,
#'                           filter_threshold = 0,
#'                           example_dataset = TRUE,
#'                           verbose = TRUE)
#'
#'scae
#'
#'rownames(scae_subset_unknown_alleles(scae))
#'
#'scae_unknown_alleles <- scae_subset_unknown_alleles(scae)
#'
#'scae_unknown_alleles
#'
#'rownames(scae_unknown_alleles)
#'
#' @export
scae_subset_unknown_alleles <- function(scae) {
  subset_rows <- stats::complete.cases(rowData(scae)$NI_I, rowData(scae)$Quant_type)
  # unknown alleles of the genes with extended quantification
  unknown <- scae[subset_rows & rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "A_unknown", ]
  return(unknown)
}

#------------------------------getter wrapper----------------------------------#

#' Subset SCAE object
#'
#' @param scae SCAE object
#' @param subset character string specifiying a data layer
#'
#' @return SCAE object
#' @export
scae_subset <- function(scae, subset = c("nonimmune", "alleles", "immune_genes", "functional", "unknown_alleles")){

  scae_sub<- switch(subset,
         "nonimmune" = get_nigenes(scae),
         "alleles" =  scae_subset_alleles(scae),
         "immune_genes" = get_agenes(scae),
         "functional" = scae_subset_functional(scae),
         "unknown_alleles" = scae_subset_unknown_alleles(scae),
         message("Invalid layer specified, Choose from `nonimmune`, `alleles`, `immune_genes`, `functional`, `unknown_alleles`"))
  return(scae_sub)
}


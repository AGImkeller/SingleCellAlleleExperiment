 ###############################################################################
 ##############---getters to use for subsampling the main assay---##############
 ###############################################################################

 #---------------------------------getters-------------------------------------#

 #####
#' Get allele rows
#'
#' @description
#' Getter function returning subsampled SCAE object with all rows containing raw allele information. These rows
#' are
#' identified by "I" in rowData(scae)$NI_I and "A" in rowData(scae)$Quant_type.
#'
#' @param scae A \code{\link{SingleCellAlleleExperiment}} object.
#'
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#'
#' @return A SingleCellAlleleExperiment object.
#'
#' @examples
#' library(SingleCellAlleleExperiment)
#'
#' example_data <- system.file("extdata", package = "SingleCellAlleleExperiment")
#'
#' scae <- readAlleleCounts(example_data,
#'                         sample_names = "example_data",
#'                         filter = "yes",
#'                         symbols = "orgdb",
#'                         exp_type = "WTA",
#'                         lookup_file = "lookup_table_HLA_only.csv",
#'                         barcode_file = "cells_x_genes.barcodes.txt",
#'                         gene_file = "cells_x_genes.genes.txt",
#'                         matrix_file = "cells_x_genes.mtx",
#'                         tag_feature_mtx = "cells_x_genes.genes.txt",
#'                         tag_feature_barcodes = "cells_x_genes.barcodes.txt",
#'                         filter_threshold = NULL)
#'
#'
#'scae
#'
#'rownames(get_alleles(scae))
#'
#'scae_alleles <- get_alleles(scae)
#'
#'scae_alleles
#'
#'rownames(scae_alleles)
#'
#' @export
get_alleles <- function(scae) {
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
#' library(SingleCellAlleleExperiment)
#'
#' example_data <- system.file("extdata", package = "SingleCellAlleleExperiment")
#'
#' scae <- readAlleleCounts(example_data,
#'                         sample_names = "example_data",
#'                         filter = "yes",
#'                         symbols = "orgdb",
#'                         exp_type = "WTA",
#'                         lookup_file = "lookup_table_HLA_only.csv",
#'                         barcode_file = "cells_x_genes.barcodes.txt",
#'                         gene_file = "cells_x_genes.genes.txt",
#'                         matrix_file = "cells_x_genes.mtx",
#'                         tag_feature_mtx = "cells_x_genes.genes.txt",
#'                         tag_feature_barcodes = "cells_x_genes.barcodes.txt",
#'                         filter_threshold = NULL)
#'
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
#' library(SingleCellAlleleExperiment)
#'
#' example_data <- system.file("extdata", package = "SingleCellAlleleExperiment")
#'
#' scae <- readAlleleCounts(example_data,
#'                         sample_names = "example_data",
#'                         filter = "yes",
#'                         symbols = "orgdb",
#'                         exp_type = "WTA",
#'                         lookup_file = "lookup_table_HLA_only.csv",
#'                         barcode_file = "cells_x_genes.barcodes.txt",
#'                         gene_file = "cells_x_genes.genes.txt",
#'                         matrix_file = "cells_x_genes.mtx",
#'                         tag_feature_mtx = "cells_x_genes.genes.txt",
#'                         tag_feature_barcodes = "cells_x_genes.barcodes.txt",
#'                         filter_threshold = NULL)
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
#' library(SingleCellAlleleExperiment)
#'
#' example_data <- system.file("extdata", package = "SingleCellAlleleExperiment")
#'
#' scae <- readAlleleCounts(example_data,
#'                         sample_names = "example_data",
#'                         filter = "yes",
#'                         symbols = "orgdb",
#'                         exp_type = "WTA",
#'                         lookup_file = "lookup_table_HLA_only.csv",
#'                         barcode_file = "cells_x_genes.barcodes.txt",
#'                         gene_file = "cells_x_genes.genes.txt",
#'                         matrix_file = "cells_x_genes.mtx",
#'                         tag_feature_mtx = "cells_x_genes.genes.txt",
#'                         tag_feature_barcodes = "cells_x_genes.barcodes.txt",
#'                         filter_threshold = NULL)
#'
#'scae
#'
#'rownames(get_func(scae))
#'
#'scae_functional_class <- get_func(scae)
#'
#'scae_functional_class
#'
#'rownames(scae_functional_class)
#'
#'
#' @export
get_func <- function(scae) {
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
#' library(SingleCellAlleleExperiment)
#'
#' example_data <- system.file("extdata", package = "SingleCellAlleleExperiment")
#'
#' scae <- readAlleleCounts(example_data,
#'                         sample_names = "example_data",
#'                         filter = "yes",
#'                         symbols = "orgdb",
#'                         exp_type = "WTA",
#'                         lookup_file = "lookup_table_HLA_only.csv",
#'                         barcode_file = "cells_x_genes.barcodes.txt",
#'                         gene_file = "cells_x_genes.genes.txt",
#'                         matrix_file = "cells_x_genes.mtx",
#'                         tag_feature_mtx = "cells_x_genes.genes.txt",
#'                         tag_feature_barcodes = "cells_x_genes.barcodes.txt",
#'                         filter_threshold = NULL)
#'
#'scae
#'
#'rownames(get_unknown(scae))
#'
#'scae_unknown_alleles <- get_unknown(scae)
#'
#'scae_unknown_alleles
#'
#'rownames(scae_unknown_alleles)
#'
#' @export
get_unknown <- function(scae) {
  subset_rows <- stats::complete.cases(rowData(scae)$NI_I, rowData(scae)$Quant_type)
  # unknown alleles of the genes with extended quantification
  unknown <- scae[subset_rows & rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "A_unknown", ]
  return(unknown)
}
 #####

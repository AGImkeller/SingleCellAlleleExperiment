 ###############################################################################
 ##############---getters to use for subsampling the main assay---##############
 ###############################################################################

 #---------------------------------getters-------------------------------------#

 #####
#' Get allele rows
#'
#' @description
#' Getter function returning subsampled SCAE-object with all rows containing raw allele information. These rows are
#' identified by "I" in rowData(scae)$NI_I and "A" in rowData(scae)$Quant_type.
#'
#' @param scae SingleCellAlleleExperiment object
#'
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#'
#' @return subsampled SingleCellAlleleExperiment object
#' @export
get_alleles <- function(scae) {
  subset_rows <- stats::complete.cases(rowData(scae)$NI_I, rowData(scae)$Quant_type)
  alleles <- scae[subset_rows & rowData(scae)$NI_I == "I" & startsWith(rowData(scae)$Quant_type,"A"), ]
  return(alleles)
}

#' Get allele gene rows
#'
#' @description
#' Getter function returning subsampled SCAE-object with all rows containing allele gene information. These rows are
#' identfied by "I" in rowData(scae)$NI_I and "G" in rowData(scae)$Quant_type.
#'
#' @param scae SingleCellAlleleExperiment object
#'
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#'
#' @return subsampled SingleCellAlleleExperiment object
#' @export
get_agenes <- function(scae) {
  subset_rows <- stats::complete.cases(rowData(scae)$NI_I, rowData(scae)$Quant_type)
  agenes <- scae[subset_rows & rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "G", ]
  return(agenes)
}

#' Get non immune rows
#'
#' @description
#' Getter function returning subsampled SCAE-object with all rows containing non immune gene information. These rows are
#' identified by "NI" in rowData(scae)$NI_I and "G" in rowData(scae)$Quant_type.
#'
#' @param scae SingleCellAlleleExperiment object
#'
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#'
#' @return subsampled SingleCellAlleleExperiment object
#' @export
get_nigenes <- function(scae) {
  subset_rows <- stats::complete.cases(rowData(scae)$NI_I, rowData(scae)$Quant_type)
  agenes <- scae[subset_rows & rowData(scae)$NI_I == "NI" & rowData(scae)$Quant_type == "G", ]
  return(agenes)
}

#' Get functional class rows
#'
#' @description
#' Getter function returning subsampled SCAE-object with all rows containing functional class information. These rows are
#' identified by "I" in rowData(scae)$NI_I and "F" in rowData(scae)$Quant_type.
#'
#' @param scae SingleCellAlleleExperiment object
#'
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#'
#' @return subsampled SingleCellAlleleExperiment object
#' @export
get_func <- function(scae) {
  subset_rows <- stats::complete.cases(rowData(scae)$NI_I, rowData(scae)$Quant_type)
  func <- scae[subset_rows & rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "F", ]
  return(func)
}

#' Get unknown allele rows
#'
#' @description
#' Getter function returning subsampled SCAE-object with all rows containing unknown allele information thats present in the data and in the proper nomenclature.
#' These rows are identified by "I" in rowData(scae)$NI_I and "A_unknown" in rowData(scae)$Quant_type.
#'
#' @param scae SingleCellAlleleExperiment object
#'
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#'
#' @return subsampled SingleCellAlleleExperiment object
#' @export
get_unknown <- function(scae) {
  subset_rows <- stats::complete.cases(rowData(scae)$NI_I, rowData(scae)$Quant_type)
  unknown <- scae[subset_rows & rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "A_unknown", ]
  return(unknown)
}
 #####
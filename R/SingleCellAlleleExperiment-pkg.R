#' SingleCellAlleleExperiment
#'
#' SingleCellAlleleExperiment: Defines a S4 class that is based on SingleCellExperiment.
#' In addition to the usual gene layer, SingleCellAlleleExperiment can also
#' store data for immune genes such as HLAs, Immunoglobulins and KIRs
#' at the allele level and at the level of functionally similar groups of immune genes.
#'
#' We have developed a workflow, that allows quantification of expression and
#' interactive exploration of donor-specific alleles of different immune genes.
#' The workflow is divided into two software packages and one additional
#' data package.
#'
#' - The scIGD software package consist of a Snakemake workflow designed to
#' automate and streamline the genotyping process for immune genes,
#' focusing on key targets such as HLAs and KIRs,
#' and enabling allele-specific quantification
#' from single-cell RNA-sequencing (scRNA-seq) data using donor-specific references.
#' For detailed information of the performed steps and
#' how to utilize this workflow, please refer to its documentation.
#'
#' - To harness the full analytical potential of the results, we have developed
#' a dedicated R package, SingleCellAlleleExperiment presented in this repository.
#' This package provides a comprehensive multi-layer data structure,
#' enabling the representation of immune genes at specific levels,
#' including alleles, genes, and groups of functionally similar genes and thus,
#' allows data analysis across these immunologically relevant, different layers of annotation.
#'
#' - The scaeData is an R/ExperimentHub data package providing datasets generated
#' and processed by the scIGD software package which can be used to explore the data
#' and potential downstream analysis workflows using the here presented
#' novel SingleCellAlleleExperiment data structure.
#' Refer to scaeData for more information regarding the available datasets
#' and source of raw data.
#'
#'
#' @seealso https://github.com/AGImkeller/scIGD/ for the definition of the
#' quantification workflow.
#' @seealso https://github.com/AGImkeller/scaeData for information regarding
#' the example data sets.
#'
#' @name SingleCellAlleleExperiment-pkg
#' @docType package
#' @keywords internal
"_PACKAGE"

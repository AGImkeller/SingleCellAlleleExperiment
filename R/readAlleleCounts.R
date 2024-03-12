
#-----------------------------read_allele_counts-------------------------------#

#####
#' Reading in allele quantification data into SingleCellAlleleExperiment object
#'
#' @description
#' Main read in function for reading in given allele quantification data and
#' loading the data into an `SingleCellAlleleExperiment` object. Input data are stored in a shared folder.
#' Expected naming scheme of the files from the data generating method:
#'
#'    * quantification matrix: `cells_x_genes.mtx`
#'    * barcode information: `cells_x_genes.barcodes.txt`
#'    * feature information: `cells_x_genes.genes.txt`
#'    * allele lookup table: `lookup_table_HLA_only`
#'
#' File identifiers can be specifically stated if the identifiers are different.
#'
#' @param samples_dir A character string determining the path to one directory containing all input files.
#' @param sample_names A character string for a sample identifier. Can be used to describe the used dataset or sample.
#' @param filter_mode A vector containing three character strings that describe different options for filtering. The value `"yes"` uses the inflection point of the knee plot to filter out low-quality cells.
#' The value `"no"` computes the knee plot and stops function execution. This mode serves as a preflight mode to observe the knee plot before filtering. The value `"custom"` allows for setting a custom threshold in the `filter_threshold` parameter.
#' @param BPPARAM A BiocParallelParam object specifying how loading should be parallelized for multiple samples.
#' @param lookup_file A character string determining the name of the lookup table file.
#' @param barcode_file A character string determining the name of the file containing the barcode identifiers.
#' @param gene_file A character string determining the name of the file containing the feature identifiers.
#' @param matrix_file A character string determining the name of the file containing the count matrix.
#' @param filter_threshold An integer value used as a threshold for filtering low-quality barcodes/cells. Standard value is `NULL` when using `filter = c("yes", "no")`. Value must be provided when using `filter = "custom"`.
#' @param example_dataset A logical parameter used when reading in example datasets from `scaeData`. `FALSE` is default. Set to `TRUE` if you want to generate an `SCAE object` with one of the available datasets in `scaeData`.
#' @param verbose A logical parameter to decide if runtime-messages should be shown during function execution. Use `FALSE` if no info runtime-messages should be shown (default), and `TRUE` for showing runtime-messages.
#'
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors DataFrame ROWNAMES
#'
#' @return A SingleCellAlleleExperiment object.
#'
#' @examples
#'
#' example_data_5k <- scaeData::scaeDataGet(dataset="pbmc_5k")
#'
#'
#' # preflight mode, not generating an SCAE object
#' # used for quality-assessment by plotting the knee plot
#' scae_preflight <- read_allele_counts(example_data_5k$dir,
#'                         sample_names = "example_data",
#'                         filter = "no",
#'                         lookup_file = "lookup_table_HLA_only.csv",
#'                         barcode_file = example_data_5k$barcodes,
#'                         gene_file = example_data_5k$features,
#'                         matrix_file = example_data_5k$matrix,
#'                         filter_threshold = NULL,
#'                         example_dataset = TRUE)
#'
#'
#' # automatic filtering mode, filtering out low-quality cells on the inflection point of the knee plot
#' scae_filtered <- read_allele_counts(example_data_5k$dir,
#'                         sample_names = "example_data",
#'                         filter = "yes",
#'                         lookup_file = "lookup_table_HLA_only.csv",
#'                         barcode_file = example_data_5k$barcodes,
#'                         gene_file = example_data_5k$features,
#'                         matrix_file = example_data_5k$matrix,
#'                         filter_threshold = NULL,
#'                         example_dataset = TRUE,
#'                         verbose = TRUE)
#'
#' scae_filtered
#'
#'
#' # custom filtering mode, setting up a custom filter threshold for filtering out
#' # low-quality cells (e.g. after using the preflight mode and assessing the knee plot)
#' scae_custom_filter <- read_allele_counts(example_data_5k$dir,
#'                         sample_names = "example_data",
#'                         filter = "custom",
#'                         lookup_file = "lookup_table_HLA_only.csv",
#'                         barcode_file = example_data_5k$barcodes,
#'                         gene_file = example_data_5k$features,
#'                         matrix_file = example_data_5k$matrix,
#'                         filter_threshold = 105,
#'                         example_dataset = TRUE)
#'
#' scae_custom_filter
#'
#'
#' @export
read_allele_counts <- function(samples_dir,
                               sample_names = names(samples_dir),
                               filter_mode = c("yes", "no", "custom"),
                               lookup_file = "lookup_table_HLA_only.csv",
                               barcode_file = "cells_x_genes.barcodes.txt",
                               gene_file = "cells_x_genes.genes.txt",
                               matrix_file = "cells_x_genes.mtx",
                               filter_threshold = NULL,
                               example_dataset = FALSE,
                               verbose = FALSE,
                               BPPARAM = BiocParallel::SerialParam()){

  rt_one_readin_start <- Sys.time()
  if (is.null(sample_names)) {
    sample_names <- samples_dir
  }

  if (filter_mode == "custom" & is.null(filter_threshold)) {
    stop("For custom filtering you need to state a integer value >0 in the 'filter_threshold' parameter.")
  }

  #reading in files
  load_out <- BiocParallel::bplapply(samples_dir,
                                     FUN = read_from_sparse_allele,
                                     barcode_file = barcode_file,
                                     gene_file = gene_file,
                                     matrix_file = matrix_file,
                                     BPPARAM = BPPARAM)

  current <- load_out[[1]]
  full_data <- current$mat
  feature_info <- current$feature_info
  cell_names <- current$cell_names
  exp_type <- current$exp_type

  #prepare colData
  cell_info_list <- S4Vectors::DataFrame(Sample = rep(sample_names,
                                                      length(cell_names)),
                                         Barcode = cell_names$V1,
                                         row.names = NULL)
  #prepare rowData
  rownames(feature_info) <- feature_info[,1]

  cnames <- cell_info_list$Barcode
  colnames(full_data) <- cnames

  full_data <- as(full_data, "CsparseMatrix")


  if (example_dataset){
    integrated_lookup_dir <- system.file("extdata", package = "SingleCellAlleleExperiment")
    lookup <- read_Lookup(integrated_lookup_dir, lookup_file)
  }else{
    lookup <- read_Lookup(samples_dir, lookup_file)
  }

  #preflight mode, only for plotting the knee plot
  if (filter_mode == "no"){
    inflection_threshold <- plot_knee(full_data, feature_info, cell_names)
    message("Suggested threshold based on inflection point is at: ", inflection_threshold, " UMI counts.")
    return()
  }

  #filtering on the inflection point shown in the advanced knee plot
  if (filter_mode == "yes"){
    inflection_threshold <- plot_knee(full_data, feature_info, cell_names)
    message("Filtering performed based on the inflection point at: ", inflection_threshold, " UMI counts.")
  }

  #putting a custom filter threshold
  if (filter_mode == "custom"){
    inflection_threshold <- filter_threshold
  }

  if (verbose){
    rt_one_readin_end <- Sys.time()
    diff_rt_one <- round(rt_one_readin_end - rt_one_readin_start, digits = 2)
    message("Runtime check (1/2) Read_in: ",      diff_rt_one, " seconds")
  }

  rt_two_scae_start <- Sys.time()
  sce <- SingleCellAlleleExperiment(assays = list(counts = full_data),
                                    rowData = feature_info,
                                    colData = cell_info_list,
                                    threshold = inflection_threshold,
                                    exp_type = exp_type,
                                    lookup = lookup,
                                    verbose = verbose)

  if (verbose){
    rt_two_scae_end <- Sys.time()
    diff_rt_two <- round(rt_two_scae_end - rt_two_scae_start, digits = 2)
    message("Runtime check (2/2) Generating SCAE completed: ",    diff_rt_two, " seconds")
    diff_rt_total <- rt_two_scae_end - rt_one_readin_start
    message("Total runtime, completed read_in, filtering and normalization and generating scae object ",       ceiling(diff_rt_total), " seconds")
  }

  return(sce)
}


# Inspired from https://github.com/MarioniLab/DropletUtils/blob/devel/R/read10xCounts.R
#' Reading in allele-aware quantification data
#'
#' @description
#' Internal function used in `read_allele_counts()` that reads in the data stated in the given directory path.
#'
#' @param path A character string determining the path to the directory containing the input files.
#' @param barcode_file A character string determining the name of the file containing the sample-tag quantification data.
#' @param gene_file A character string determining the name of the file containing the feature identifiers.
#' @param matrix_file A character string determining the name of the file containing the count matrix.
#'
#' @importFrom utils read.delim read.csv
#' @importFrom Matrix readMM t
#'
#' @return A list with three data.frames containing the input data information.
read_from_sparse_allele <- function(path,
                                    barcode_file,
                                    gene_file,
                                    matrix_file){

  barcode_loc <- file.path(path, barcode_file)
  feature_loc <- file.path(path, gene_file)
  matrix_loc  <- file.path(path, matrix_file)

  feature_info <- utils::read.delim(feature_loc, header = FALSE)
  cell_names   <- utils::read.csv(barcode_loc, sep = "", header = FALSE)
  mat          <- Matrix::readMM(matrix_loc)

  possible_names <- c("Ensembl_ID", "Symbol")

  if (grepl("ENS", feature_info$V1[1])){
    exp_type = "ENS"
    colnames(feature_info) <- possible_names[1]
  }else{
    exp_type = "noENS"
    colnames(feature_info) <- possible_names[2]
  }

  list(mat =  Matrix::t(mat),
       cell_names =  cell_names,
       feature_info = feature_info,
       exp_type = exp_type)
}

#' Read in allele lookup
#'
#' @description
#' Internal function used in `read_allele_counts()` to read in the allele lookup table.
#'
#' @param path A character string determining the path to the directory containing the input files.
#' @param lookup_file A character string determining the name of the lookup table file.
#'
#' @importFrom utils read.csv
#'
#' @return A data.frame containing a representation of the lookup table.
read_Lookup <- function(path, lookup_file){
  lookup_loc <- file.path(path, lookup_file)
  lookup <- utils::read.csv(lookup_loc)
  lookup
}

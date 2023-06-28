################################################################################
#######---Read in function creating SingleCellAlleleExperiment object---########
###############---file contains main function and its helpers---###############
################################################################################

#-----------------------------readAlleleCounts---------------------------------#

#####
#' Reading in allele quantification data into SingleCellAlleleExperiment object
#'
#' @description
#' Main read in function for reading in given allele quantification data and
#' loading the data into an SingleCellAlleleExperiment object. Input data are stored in a shared folder.
#' Expected naming scheme of the files: quantification matrix: "matrix.mtx"
#'                                      barcode information: "barcodes.txt"
#'                                      feature information: "features.txt"
#'                                      allele lookup table: "lookup_table_HLA_only"
#'
#' @param samples character string input containing the path to the directory containing the
#'   input files
#' @param sample.names character string for a sample_name identifier
#' @param col.names binary variable indicating whether quantification assay should contain column names
#' @param compressed binary variable whether the input files are .gz compressed
#' @param filter_threshold count threshold for filtering barcodes/cells
#' @param BPPARAM A BiocParallelParam object specifying how loading should be parallelized for multiple samples
#'
#' @return SingleCellAlleleExperiment object
#'
# #example
#'
#' @export
readAlleleCounts <- function (samples,
                              sample.names = names(samples),
                              col.names = TRUE,
                              compressed = NULL,
                              filter_threshold = 0,
                              BPPARAM = BiocParallel::SerialParam()){

  rt_one_readin_start <- Sys.time()
  if (is.null(sample.names)) {
    sample.names <- samples
  }

  load.out <- BiocParallel::bplapply(samples,
                                     FUN = read_from_sparse_allele,
                                     compressed = compressed,
                                     BPPARAM = BPPARAM)

  current <- load.out[[1]]
  full_data <- current$mat
  feature_info_list <- current$feature.info
  cell.names <- current$cell.names
  cell_info_list <- S4Vectors::DataFrame(Sample = rep(sample.names,
                                                      length(cell.names)),
                                         Barcode = cell.names$V1,
                                         row.names = NULL)
  feature_info <- feature_info_list
  S4Vectors::ROWNAMES(feature_info) <- feature_info$ID

  if (col.names) {
    cnames <- cell_info_list$Barcode
    colnames(full_data) <- cnames
  }

  full_data <- as(full_data, "CsparseMatrix")
  lookup <- readLookup(samples)

  #####
  rt_one_readin_end <- Sys.time()
  diff_rt_one <- rt_one_readin_end - rt_one_readin_start
  print(paste("Runtime check (1/2) Read_in:",      diff_rt_one))
  #####

  rt_two_scae_start <- Sys.time()
  sce <- SingleCellAlleleExperiment(assays = list(counts = full_data),
                                    rowData = feature_info,
                                    colData = cell_info_list,
                                    threshold = filter_threshold,
                                    lookup = lookup)
  #####
  rt_two_scae_end <- Sys.time()
  diff_rt_two <- rt_two_scae_end - rt_two_scae_start
  print(paste("Runtime check (2/2) Generating SCAE completed:",     diff_rt_two))
  diff_rt_total <- rt_two_scae_end - rt_one_readin_start
  print(paste("Total runtime, completed read_in, filtering and normalization and generating scae object",     diff_rt_total))
  #####
  return(sce)
}

# Inspired from https://github.com/MarioniLab/DropletUtils/blob/devel/R/read10xCounts.R
# TODO: das hier stimmt nicht, check Server-version fuer fertige Dokumentation
#
#' `read_from_sparse_allele` read in function for reading in raw data.
#'
#'
#' @param path character string input containing the path to the directory containing the
#' input files
#' @param compressed binary classification if the input data are .gz compressed
#'
#' @import utils
#' @importFrom Matrix readMM
#' @importFrom Matrix t
# @import SingleCellExperiment
#'
#' @return list with the read_in data sorted into different slots
read_from_sparse_allele <- function(path, compressed = NULL){
  #defining path
  barcode.loc <- file.path(path, "barcodes.txt")
  feature.loc <- file.path(path, "features.txt")
  matrix.loc  <- file.path(path, "matrix.mtx")
  #check for compression
  barcode.loc <- check_for_compressed(barcode.loc, compressed)
  feature.loc <- check_for_compressed(feature.loc, compressed)
  matrix.loc  <- check_for_compressed(matrix.loc, compressed)
  #read in the files
  feature.info <- utils::read.delim(feature.loc, header = FALSE)
  cell.names   <- utils::read.csv(barcode.loc, sep = "", header = FALSE)
  mat <- Matrix::readMM(matrix.loc)

  possible.names <- c("ID", "Symbol", "Type", "Chromosome", "Start", "End")
  colnames(feature.info) <- utils::head(possible.names, ncol(feature.info))

  list(mat = Matrix::t(mat),
       cell.names = cell.names,
       feature.info = feature.info)
}

# Taken from https://github.com/MarioniLab/DropletUtils/blob/devel/R/read10xCounts.R
#' Check if files in given directory are gzipped
#'
#' @description
#' Internal function used in `read_from_sparse_allele()` that checks whether the input files are gzipped.
#'
#' @param path file path of the directory containing the input files as character string; taken from
#' readAllleleCounts()
#' @param compressed binary variable taken from readAlleleCounts() indicating if the user files are gzipped
#' @param error binary variable taken from readAlleleCounts() indicating if there should be an error message, when the gzipped files
#' cannot be found in the given directory
#'
#' @return character string adding ".gz" to the filenames
check_for_compressed <- function(path, compressed, error=TRUE) {
  original <- path
  if (isTRUE(compressed)) {
    path <- paste0(path, ".gz")
  }else if (is.null(compressed) && !file.exists(path)) {
    path <- paste0(path, ".gz")
    if (error && !file.exists(path)) {
      # Explicit error here to avoid users getting confused.
      stop(sprintf("cannot find '%s' or its gzip-compressed form", original))
    }
  }
  path
}

#' Read in allele lookup
#'
#' @description
#' Internal function used in `readAlleleCounts()` to read in the allele lookup table.
#'
#' @param path file path of the directory containing the input files as character string
#'
#' @return lookup table
readLookup <- function(path){
  lookup.loc <- file.path(path, "lookup_table_HLA_only.csv")
  lookup <- utils::read.csv(lookup.loc)
  lookup
}
#####

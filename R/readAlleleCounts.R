################################################################################
#######---Read in function creating SingleCellAlleleExperiment object---########
###############---file contains main function and its helpers---###############
################################################################################

#-----------------------------readAlleleCounts---------------------------------#

#####
#' Reading in allele quantification data into a SingleCellAlleleExperiment object
#'
#' @description
#' Main read in function for reading in given allele quantification data and
#' loading the data into a `SingleCellAlleleExperiment` object. Input data are stored in a shared folder with expected file names.
#' Expected naming scheme of the files:
#'
#'    * quantification matrix: `cells_x_genes.matrix.mtx`
#'    * barcode information:   `cells_x_genes.barcodes.txt`
#'    * feature information:   `cells_x_genes.genes.txt`
#'    * lookup table:          `lookup_table_HLA_only`
#'
#' @param samples character string of the file path to the directory containing the expected input files
#' @param sample_names character string for a sample_name identifier. If left empty it will take the path to the files as a sample identifier.
#' @param filter_threshold positive integer value used as count threshold for filtering barcodes/cells. Default value is set at `filter_threshold = 0`.
#' @param BPPARAM A BiocParallelParam object specifying how loading should be parallelized for multiple samples. Default is `BiocParallel::SerialParam()`.
#' @param exp_type character string to pass information about the used experimental approach. Either use `exp_type = "WTA"` or `exp_type = "Amplicon"`.
#' @param symbols character string to determine which database-function to use to retrieve the NCBI gene names if `exp_type = "WTA"`. Default is `symbols = "biomart"`. You can choose between `exp_type = c("biomart", "orgdb")`. `orgdb` parameter is suggested for offline-usage.
#'
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors DataFrame ROWNAMES
#'
#' @return SingleCellAlleleExperiment object
#'
#' @export
readAlleleCounts <- function (samples,
                              sample_names = names(samples),
                              filter_threshold = 0,
                              exp_type = c("WTA", "Amplicon"),
                              symbols = NULL,
                              BPPARAM = BiocParallel::SerialParam()){

  rt_one_readin_start <- Sys.time()
  if (is.null(sample_names)) {
    sample_names <- samples
  }

  if (is.null(symbols)) {
    symbols <- "biomart"
  } else {
    if (!symbols %in% c("biomart", "orgdb")) {
      stop("Invalid value for symbols parameter. Allowed values are 'biomaRt' and 'orgdb'.")
    }
  }

  #reading in files
  load_out <- BiocParallel::bplapply(samples,
                                     FUN = read_from_sparse_allele,
                                     exp_type = exp_type,
                                     BPPARAM = BPPARAM)

  current <- load_out[[1]]
  full_data <- current$mat
  feature_info <- current$feature_info
  cell_names <- current$cell_names

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
  lookup <- readLookup(samples, exp_type)

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
                                    exp_type = exp_type,
                                    symbols = symbols,
                                    lookup = lookup)
  if (exp_type == "Amplicon"){
    rt_six_scae_start <- Sys.time()
    sce <- add_sample_tags(samples, sce)
    #####
    rt_six_scae_end <- Sys.time()
    diff_rt_six <- rt_six_scae_end - rt_six_scae_start
    print(paste("     Generating SCAE (6/X) adding sample tags:", diff_rt_six))
    #####
  }

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
#' Reading in allele-aware quantification data
#'
#' @description
#' Internal function used in `readAlleleCounts()` that reads in the data stated in the given directory path.
#'
#'
#' @param path character string of the file path to the directory containing the expected input files
#' @param exp_type character string to pass information about the used experimental approach. Either use `exp_type = "WTA"` or `exp_type = "Amplicon"`.
#'
#' @importFrom utils read.delim read.csv
#' @importFrom Matrix readMM t
#'
#' @return list with the read_in data sorted into different slots
read_from_sparse_allele <- function(path, exp_type = exp_type){
  barcode_loc <- file.path(path, "cells_x_genes.barcodes.txt")
  feature_loc <- file.path(path, "cells_x_genes.genes.txt")
  matrix_loc  <- file.path(path, "cells_x_genes.mtx")

  feature_info <- utils::read.delim(feature_loc, header = FALSE)
  cell_names   <- utils::read.csv(barcode_loc, sep = "", header = FALSE)
  mat          <- Matrix::readMM(matrix_loc)

  possible_names <- c("Ensembl.ID", "Symbol")

  if (exp_type == "WTA"){
    colnames(feature_info) <- possible_names[1]
  }else {
    colnames(feature_info) <- possible_names[2]
  }

  list(mat = Matrix::t(mat),
       cell_names = cell_names,
       feature_info = feature_info)
}

#' Read in lookup table
#'
#' @description
#' Internal function used in `readAlleleCounts()` to read in the lookup table.
#'
#' @param path character string of the file path to the directory containing the expected input files
#' @param exp_type character string to pass information about the used experimental approach. Either use `exp_type = "WTA"` or `exp_type = "Amplicon"`.
#'
#' @importFrom utils read.csv
#'
#' @return lookup table
readLookup <- function(path, exp_type){
  # check for different naming schemes of the lookup table for the different experimental approaches
  if (exp_type == "WTA"){
    lookup_loc <- file.path(path, "lookup_table_HLA_only.csv")
    lookup <- utils::read.csv(lookup_loc)
  }else if (exp_type == "Amplicon"){
    lookup_loc <- file.path(path, "lookup_table_HLA_amplicon.csv")
    lookup <- utils::read.csv(lookup_loc)
  }
  lookup
}
#####

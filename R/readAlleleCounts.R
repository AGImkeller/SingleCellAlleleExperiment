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
#' loading the data into an `SingleCellAlleleExperiment` object. Input data are stored in a shared folder.
#' Expected naming scheme of the files:
#'
#'    * quantification matrix: `matrix.mtx`
#'    * barcode information: `barcodes.txt`
#'    * feature information: `features.txt`
#'    * allele lookup table: `lookup_table_HLA_only`
#'
#' @param samples character string input containing the path to the directory containing the
#'   input files
#' @param sample.names character string for a sample_name identifier
#' @param filter_threshold count threshold for filtering barcodes/cells
#' @param BPPARAM A BiocParallelParam object specifying how loading should be parallelized for multiple samples
#' @param exp_type either `WTA` or `Amplicon` depending on the used experiments technology
#' @param symbols identifier used to choose which database-function to use to retrieve the ncbi gene names
#'
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors DataFrame ROWNAMES
#'
#' @return SingleCellAlleleExperiment object
#'
#' @export
readAlleleCounts <- function (samples,
                              sample.names = names(samples),
                              filter_threshold = 0,
                              exp_type = c("WTA", "Amplicon"),
                              symbols = NULL,
                              lookup_file = "lookup_table_HLA_only.txt",
                              barcode_file = "cells_x_genes.barcodes.txt",
                              gene_file = "cells_x_genes.genes.txt",
                              matrix_file = "cells_x_genes.mtx",
                              tag_feature_mtx = "cells_x_features.mtx", 
                              tag_feature_barcodes = "cells_x_features.barcodes.txt",
                              BPPARAM = BiocParallel::SerialParam()){

  rt_one_readin_start <- Sys.time()
  if (is.null(sample.names)) {
    sample.names <- samples
  }

  if (is.null(symbols)) {
    symbols <- "biomart"
  } else {
    if (!symbols %in% c("biomart", "orgdb")) {
      stop("Invalid value for symbols parameter. Allowed values are 'biomaRt' and 'orgdb'.")
    }
  }

  #reading in files
  load.out <- BiocParallel::bplapply(samples,
                                     FUN = read_from_sparse_allele,
                                     exp_type = exp_type,
                                     barcode_file = barcode_file,
                                     gene_file = gene_file,
                                     matrix_file = matrix_file,
                                     BPPARAM = BPPARAM)

  current <- load.out[[1]]
  full_data <- current$mat
  feature_info <- current$feature.info
  cell.names <- current$cell.names

  #prepare colData
  cell_info_list <- S4Vectors::DataFrame(Sample = rep(sample.names,
                                                      length(cell.names)),
                                         Barcode = cell.names$V1,
                                         row.names = NULL)
  #prepare rowData
  rownames(feature_info) <- feature_info[,1]

  cnames <- cell_info_list$Barcode
  colnames(full_data) <- cnames

  full_data <- as(full_data, "CsparseMatrix")
  lookup <- readLookup(samples, exp_type, lookup_file)

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
    sce <- add_sample_tags(samples, sce, tag_feature_mtx, tag_feature_barcodes)
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
#' @param path character string input containing the path to the directory containing the
#' input files
#' @param exp_type either `WTA` or `Amplicon` depending on the used experiments technology
#'
#' @importFrom utils read.delim read.csv
#' @importFrom Matrix readMM t
#'
#' @return list with the read_in data sorted into different slots
read_from_sparse_allele <- function(path, exp_type = exp_type, 
                                    barcode_file,
                                    gene_file,
                                    matrix_file){
  # this needs to be provided as an input argument, not hardcoded
  barcode.loc <- file.path(path, barcode_file)
  feature.loc <- file.path(path, gene_file)
  matrix.loc  <- file.path(path, matrix_file)

  feature.info <- utils::read.delim(feature.loc, header = FALSE)
  cell.names   <- utils::read.csv(barcode.loc, sep = "", header = FALSE)
  mat          <- Matrix::readMM(matrix.loc)
  
  
  # call the kneeplot function somewhere here. the chosen kneepoint can be selected automatically
  possible.names <- c("Ensembl.ID", "Symbol")

  if (exp_type == "WTA"){
    colnames(feature.info) <- possible.names[1]
  }else {
    colnames(feature.info) <- possible.names[2]
  }

  list(mat = Matrix::t(mat),
       cell.names = cell.names,
       feature.info = feature.info)
}

#' Read in allele lookup
#'
#' @description
#' Internal function used in `readAlleleCounts()` to read in the allele lookup table.
#'
#' @param path file path of the directory containing the input files as character string
#' @param exp_type either "WTA" or "Amplicon" depending on the used experiments technology
#'
#' @importFrom utils read.csv
#'
#' @return lookup table
readLookup <- function(path, exp_type, lookup_file){
    lookup.loc <- file.path(path, lookup_file)
    lookup <- utils::read.csv(lookup.loc)
  lookup
}
#####

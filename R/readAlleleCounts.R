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
#' @param sample_names character string for a sample_name identifier
#' @param filter count threshold for filtering barcodes/cells
#' @param BPPARAM A BiocParallelParam object specifying how loading should be parallelized for multiple samples
#' @param exp_type either `WTA` or `Amplicon` depending on the used experiments technology
#' @param symbols identifier used to choose which database-function to use to retrieve the ncbi gene names
#' @param lookup_file add text here
#' @param barcode_file add text here
#' @param gene_file add text here
#' @param matrix_file add text here
#' @param tag_feature_mtx add text here
#' @param tag_feature_barcodes add text here
#' @param filter_threshold add text here
#'
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors DataFrame ROWNAMES
#'
#' @return SingleCellAlleleExperiment object
#'
#' @export
readAlleleCounts <- function (samples,
                              sample_names = names(samples),
                              filter = c("yes", "no", "custom" ),
                              exp_type = c("WTA", "Amplicon"),
                              symbols = NULL,
                              lookup_file = "lookup_table_HLA_only.txt",
                              barcode_file = "cells_x_genes.barcodes.txt",
                              gene_file = "cells_x_genes.genes.txt",
                              matrix_file = "cells_x_genes.mtx",
                              tag_feature_mtx = "cells_x_features.mtx",
                              tag_feature_barcodes = "cells_x_features.barcodes.txt",
                              filter_threshold = NULL,
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

  if (filter == "custom" & is.null(filter_threshold)) {
    stop("For custom filtering you need to state a integer value >0 in the 'filter_threshold' parameter.")
  }




  #reading in files
  load_out <- BiocParallel::bplapply(samples,
                                     FUN = read_from_sparse_allele,
                                     exp_type = exp_type,
                                     barcode_file = barcode_file,
                                     gene_file = gene_file,
                                     matrix_file = matrix_file,
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
  lookup <- readLookup(samples, exp_type, lookup_file)



  #preflight mode, only for plotting the knee plots
  if (filter == "no"){
    inflection_threshold <- plotKnee(full_data, feature_info, cell_names)
    cat("suggested threshold based on inflection point is at: ", inflection_threshold, " UMI counts.\n")
    stop()
  }

  #filtering on the inflection point of the shown in the advanced knee plot
  if (filter == "yes"){
    inflection_threshold <- plotKnee(full_data, feature_info, cell_names)
    cat("Filtering performed based on the inflection point at: ", inflection_threshold, " UMI counts.\n")

  }

  #putting a custom filtering threshold
  if (filter == "custom"){
    inflection_threshold <- filter_threshold
  }


  #####
  rt_one_readin_end <- Sys.time()
  diff_rt_one <- rt_one_readin_end - rt_one_readin_start
  print(paste("Runtime check (1/2) Read_in:",      diff_rt_one))
  #####

  rt_two_scae_start <- Sys.time()
  sce <- SingleCellAlleleExperiment(assays = list(counts = full_data),
                                    rowData = feature_info,
                                    colData = cell_info_list,
                                    threshold = inflection_threshold,
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
#' @param barcode_file add text here
#' @param gene_file add text here
#' @param matrix_file add text here
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
  barcode_loc <- file.path(path, barcode_file)
  feature_loc <- file.path(path, gene_file)
  matrix_loc  <- file.path(path, matrix_file)

  feature_info <- utils::read.delim(feature_loc, header = FALSE)
  cell_names   <- utils::read.csv(barcode_loc, sep = "", header = FALSE)
  mat          <- Matrix::readMM(matrix_loc)


  # call the kneeplot function somewhere here. the chosen kneepoint can be selected automatically
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

#' Read in allele lookup
#'
#' @description
#' Internal function used in `readAlleleCounts()` to read in the allele lookup table.
#'
#' @param path file path of the directory containing the input files as character string
#' @param exp_type either "WTA" or "Amplicon" depending on the used experiments technology
#' @param lookup_file add text here
#'
#' @importFrom utils read.csv
#'
#' @return lookup table
readLookup <- function(path, exp_type, lookup_file){
    lookup_loc <- file.path(path, lookup_file)
    lookup <- utils::read.csv(lookup_loc)
  lookup
}
#####

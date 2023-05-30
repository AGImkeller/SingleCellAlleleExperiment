#library(BiocParallel)
#library(S4Vectors)
#library(SingleCellExperiment)
#library(rhdf5)
#library(HDF5Array)
#library(SingleCellAlleleExperiment)
#library(Matrix)
#library(tools)
#library(DropletUtils)
#library(BiocGenerics)
#library(png)
#library(spatialLIBD)
#library(SpatialExperiment


#----------------------------Helper_functions-----------------------------------
#####
#' `.check_for_compressed()` internal function that checks whether the input files are gzipped
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
  } else if (is.null(compressed) && !file.exists(path)) {
    path <- paste0(path, ".gz")
    if (error && !file.exists(path)) {
      # Explicit error here to avoid users getting confused.
      stop(sprintf("cannot find '%s' or its gzip-compressed form", original))
    }
  }
  path
}

#' `check_txt_tsv` internal function that checks if the input files are .txt or .tsv files; might be needed when
#' spatial data get incorporated
#'
#' @param path file path of the directory containing the input files as character string; taken from
#' readAllleleCounts()
#' @param filename character string containing the names of the files to be checked
#'
#' @return file name either appended with .txt or .tsv - depending on what could be found
#' in the input directory
check_txt_tsv <- function(path, filename){
  files <- list.files(path, pattern = paste0("^", filename, "\\."), full.names = TRUE)

  if (length(files) == 0) {
    print(paste0("No file with name ", filename, " found"))
  } else {
    for (file in files) {
      ext <- file_ext(file)
      if (ext == "txt") {
        #print(paste0(file, " is a text file"))
        filename <- paste0(filename, ".txt")
      } else if (ext == "tsv") {
        #print(paste0(file, " is a TSV file"))
        filename <- paste0(filename, ".tsv")
      } else {
        print(paste0(file, " is not a text or TSV file"))
      }
    }
  }
  return(filename)
}
#####
#-------------------------------------------------------------------------------
#-----------------------------Main_functions------------------------------------
#####

#' `allele_loader` internal wrapper function for choosing the right read-in function
#'
#' @param path character string input containing the path to the directory containing the
#' input files; helper function used in readAlleleCounts()
#' @param type character vector containing options for the input data_type
#' @param compressed binary classification if the input data are .gz compressed
#'
#' @return depending on the chosen type, fitting read_in function will be used
allele_loader <- function(path, type, compressed){
  temp_type <- type
  if (temp_type == "auto"){
    temp_type <- if (grepl("\\.h5", path)) "HDF5" else "sparse"
  }
  if (temp_type == "sparse"){
    read_from_sparse_allele(path, compressed = compressed)
  }else if (temp_type == "prefix"){
    read_from_sparse_allele(path, is.prefix = TRUE, compressed = compressed)
  }
}

#' `readLookup` internal function to read in the lookup table
#'
#' @param path character string input containing the path to the directory containing the
#' input files; helper function used in readAlleleCounts()
#' @importFrom utils read.csv
#'
#' @return lookup table
readLookup <- function(path){
  lookup.loc <- file.path(path, "lookup_table_HLA_only.csv")
  lookup <- read.csv(lookup.loc)
  lookup
}


#' `readAlleleCounts` main read_in function for reading in the allele quantification data and
#' loading the data into an SingleCellAlleleExperiment object. Input data are stored in a shared folder.
#' Expected naming scheme of the files: quantification matrix: "matrix.mtx"
#'                                      barcode information: "barcodes.txt"
#'                                      feature information: "features.txt"
#'                                      allele lookup table: "lookup_table_HLA_only"
#'
#' @param samples character string input containing the path to the directory containing the
#' input files
#' @param sample.names character string for a sample_name identifier
#' @param col.names binary variable indicating whether quantification assay should contain column names
#' @param type vector of character strings containing options for the input data_type
#' @param compressed binary variable whether the input files are .gz compressed
#' @param BPPARAM A BiocParallelParam object specifying how loading should be parallelized for multiple samples
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @import BiocParallel
#' @import HDF5Array
#' @import tools
#' @import methods
#' @import DelayedArray
#'
#' @return SingleCellAlleleExperiment object
#' @export
readAlleleCounts <- function (samples,
                              sample.names = names(samples),
                              col.names = TRUE,
                              type = c("auto", "sparse", "HDF5", "prefix"),
                              compressed = NULL,
                              BPPARAM = SerialParam()){
  type = match.arg(type)
  if (is.null(sample.names)) {
    sample.names <- samples
  }
  rt_one_readin_start <- Sys.time()
  load.out <- bplapply(samples, FUN = allele_loader,
                       type = type, compressed = compressed,
                       BPPARAM = BPPARAM)
  rt_one_readin_end <- Sys.time()
  diff_rt_one <- rt_one_readin_end - rt_one_readin_start
  print(paste("Runtime check (1/2) Read_in:",      diff_rt_one))

  #determining how many samples are contained in "samples"
  #predefine the dataslots in the corresponding length
  nsets <- length(samples)
  full_data <- vector("list", nsets)
  feature_info_list <- vector("list", nsets)
  cell_info_list <- vector("list", nsets)
  #load the data from load.out into the predefined slots
  for (i in 1:nsets) {
    current <- load.out[[i]]
    full_data[[i]] <- current$mat
    feature_info_list[[i]] <- current$feature.info
    cell.names <- current$cell.names
    #creates DataFrame with info to put in colData of sce
    cell_info_list[[i]] <- S4Vectors::DataFrame(Sample = rep(sample.names[i],length(cell.names)),
                                     Barcode = cell.names$V1, row.names = NULL)
  }
  if (nsets > 1 && length(unique(feature_info_list)) != 1L) {
    stop("gene information differs between runs")
  }
  #assignment of variables for object assignment, after variables contain data for all given samples
  feature_info <- feature_info_list[[1]]
  S4Vectors::ROWNAMES(feature_info) <- feature_info$ID
  full_data <- do.call(cbind, full_data)
  cell_info <- do.call(rbind, cell_info_list)

  #add the barcodes as columnames to the assay_data
  if (col.names) {
    if (nsets == 1L) {
      cnames <- cell_info$Barcode
    }
    else {
      sid <- rep(seq_along(cell_info_list), vapply(cell_info_list, nrow, 1L))
      cnames <- paste0(sid, "_", cell_info$Barcode)
    }
    colnames(full_data) <- cnames
  }
  full_data <- DelayedArray(full_data)
  lookup <-  readLookup(samples)

  rt_two_scae_start <- Sys.time()
  sce <- SingleCellAlleleExperiment(assays=list(counts = full_data),
                                    rowData = feature_info,
                                    colData = cell_info, lookup = lookup)
  rt_two_scae_end <- Sys.time()
  diff_rt_two <- rt_two_scae_end - rt_two_scae_start
  print(paste("Runtime check (2/2) Generating SCAE completed:",     diff_rt_two))
  diff_rt_total <- rt_two_scae_end - rt_one_readin_start
  print(paste("Total runtime, completed read_in and generating scae object",     diff_rt_total))

  return(sce)
}


#' `read_from_sparse_allele` read in function for reading in raw data.
#'
#'
#' @param path character string input containing the path to the directory containing the
#' input files
#' @param is.prefix binary if the input data contains sample_prefix information
#' @param compressed binary classification if the input data are .gz compressed
#'
#' @import SingleCellExperiment
#' @import utils
#' @importFrom Matrix readMM
#'
#' @return list with the read_in data sorted into different slots
read_from_sparse_allele <- function(path, is.prefix = FALSE, compressed = NULL){
  Fun <- if (is.prefix) paste0 else file.path
  #check if barcodes/features are .txt or .tsv files
  barcode <- check_txt_tsv(path, "barcodes")
  features <- check_txt_tsv(path, "features")

  barcode.loc <- Fun(path, barcode)
  feature.loc <- Fun(path, features)
  matrix.loc <- Fun(path, "matrix.mtx")

  barcode.loc <- check_for_compressed(barcode.loc, compressed)
  feature.loc <- check_for_compressed(feature.loc, compressed)
  matrix.loc <- check_for_compressed(matrix.loc, compressed)

  feature.info <- read.delim(feature.loc, header = FALSE)
  cell.names <- read.csv(barcode.loc, sep = "", header = FALSE)

  possible.names <- c("ID", "Symbol", "Type", "Chromosome", "Start", "End")
  colnames(feature.info) <- head(possible.names, ncol(feature.info))

  mat = Matrix::readMM(matrix.loc)

  if (barcode == "barcodes.txt"){
    mat = t(mat)
  }
  list(
    mat = mat,
    cell.names = cell.names,
    feature.info = feature.info
  )
}
#####

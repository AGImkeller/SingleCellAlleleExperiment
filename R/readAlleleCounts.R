library(BiocParallel)
library(Matrix)
library(S4Vectors)
library(SingleCellExperiment)
library(rhdf5)
library(HDF5Array)
library(SingleCellAlleleExperiment)
library(Matrix)
library(tools)
#library(DropletUtils)
#library(BiocGenerics)
#library(png)
#library(spatialLIBD)
#library(SpatialExperiment

#this version uses the SingleCellAlleleExperiment Constructor
#----------------------------Helper_functions-----------------------------------
#####
#if the input.txts are saved as .gz
.check_for_compressed <- function(path, compressed, error=TRUE) {
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

#more or less for single/spatial distinction // spatial data from 10X Visium are saved
#as .tsv, not as txt
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
#wrapper choosing the right read_in function
#raw or hdf5 data
.allele_loader <- function(run, type, compressed){
  temp_type <- type
  if (temp_type == "auto"){
    temp_type <- if (grepl("\\.h5", run)) "HDF5" else "sparse"
  }
  if (temp_type == "sparse"){
    read_from_sparse_allele(run, compressed = compressed)
  }else if (temp_type == "prefix"){
    read_from_sparse_allele(run, is.prefix = TRUE, compressed = compressed)
  }else {
    read_from_hdf5_allele(run)
  }
}

#main function -use this-
readAlleleCounts <- function (samples,
                              sample.names = names(samples),
                              col.names = TRUE,
                              type = c("auto", "sparse", "HDF5", "prefix"),
                              compressed = NULL,
                              images = "lowres", load = TRUE,
                              BPPARAM = SerialParam(),
                              lookup){
  type = match.arg(type)
  if (is.null(sample.names)) {
    sample.names <- samples
  }
  #load.out for loading the different data
  load.out <- bplapply(samples, FUN = .allele_loader,
                       type = type, compressed = compressed,
                       BPPARAM = BPPARAM)
  #determining how many samples are contained in samples
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
    cell_info_list[[i]] <- DataFrame(Sample = rep(sample.names[i],length(cell.names)),
                                     Barcode = cell.names$V1, row.names = NULL)
  }
  if (nsets > 1 && length(unique(feature_info_list)) != 1L) {
    stop("gene information differs between runs")
  }
  #assignment of variables for object assignment, after variables contain data for all given samples
  feature_info <- feature_info_list[[1]]
  ROWNAMES(feature_info) <- feature_info$ID
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
  #write assay as an DelayedArray-Object
  full_data <- DelayedArray(full_data)

  sce <- SingleCellAlleleExperiment(assays=list(counts = full_data),
                                    rowData = feature_info,
                                    colData = cell_info, lookup = lookup)
  return(sce)
}

read_from_sparse_allele <- function(path, is.prefix = FALSE, compressed = NULL){
  Fun <- if (is.prefix) paste0 else file.path
  #check if barcodes/features
  barcode <- check_txt_tsv(path, "barcodes")
  features <- check_txt_tsv(path, "features")

  barcode.loc <- Fun(path, barcode)
  feature.loc <- Fun(path, features)
  matrix.loc <- Fun(path, "matrix.mtx")

  barcode.loc <- .check_for_compressed(barcode.loc, compressed)
  feature.loc <- .check_for_compressed(feature.loc, compressed)
  matrix.loc <- .check_for_compressed(matrix.loc, compressed)

  feature.info <- read.delim(feature.loc, header = FALSE)
  cell.names <- read.csv(barcode.loc, sep = "", header = FALSE)

  possible.names <- c("ID", "Symbol", "Type", "Chromosome", "Start", "End")
  colnames(feature.info) <- head(possible.names, ncol(feature.info))

  mat = readMM(matrix.loc)

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
#-------------------------------------------------------------------------------

####------------------------------test-case-------------------------------------
lookup_hla <- read.csv("C:/Users/Jonas/OneDrive/Desktop/input_test_files/counts_unfiltered___/lookup_table_HLA_only.csv")
sample_x7_auto_sparse_txt <- "C:/Users/Jonas/OneDrive/Desktop/input_test_files/counts_unfiltered___"

sce_test7_sparse_auto_txt <- readAlleleCounts(sample_x7_auto_sparse_txt, sample.names = "sample_x",
                                              type = "auto", col.names = TRUE, lookup=lookup_hla)

sce_test7_sparse_auto_txt

#counts(sce_test7_sparse_auto_txt)
rowData(sce_test7_sparse_auto_txt)

get_alleles(sce_test7_sparse_auto_txt)
rowData(get_alleles(sce_test7_sparse_auto_txt))
counts(get_alleles(sce_test7_sparse_auto_txt))
rownames(get_alleles(sce_test7_sparse_auto_txt))

get_agenes(sce_test7_sparse_auto_txt)
rowData(get_agenes(sce_test7_sparse_auto_txt))
counts(get_agenes(sce_test7_sparse_auto_txt))

get_func(sce_test7_sparse_auto_txt)
rowData(get_func(sce_test7_sparse_auto_txt))
counts(get_func(sce_test7_sparse_auto_txt))


counts(sce_test7_sparse_auto_txt)
colData(sce_test7_sparse_auto_txt)


#####



lookup_hla <- read.csv("C:/Users/Jonas/OneDrive/Desktop/input_test_files/counts_unfiltered___/lookup_table_HLA_only.csv")
lookup_hla

lookup_hla[grepl("DRB1", lookup_hla$Allele, fixed = TRUE),]



################################################################################
##############---Class definition of SingleCellAlleleExperiment---##############
############---file contains all helpers to transform sce to scae---############
################################################################################

#---------SingleCellAlleleExperiment class definition and constructor----------#

#####
#' SingleCellAlleleExperiment-class definition
#'
#' @description
#' Defining the SingleCellAlleleExperiment-class derived from SingleCellExperiment-class.
#'
#' @import methods
# @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @import SingleCellExperiment
#'
#' @return definition for the scae class
.scae <- setClass("SingleCellAlleleExperiment", contains = "SingleCellExperiment")

#' Constructor SingleCellAlleleExperiment-class
#'
#' @description
#' Constructor for the SingleCellAllelExperiment-class.
#' Constructor is used in the read in function `readAlleleCounts()`. Performing all necessary steps to transform
#' a SingleCellExperiment object into the extended SingleCellAlleleExperiment object. SCAE objects
#' contain an extended count-assay aswell as extended rowData.
#'
#' @param ... parameters to pass to SingleCellExperiment constructor
#' @param lookup allele lookup file
#' @param threshold count threshold for filtering barcodes/cells
#'
#' @return SingleCellAlleleExperiment object
#' @export
SingleCellAlleleExperiment <- function(..., lookup, threshold){
  sce <- SingleCellExperiment::SingleCellExperiment(...)

  rt_scae_lookup_start <- Sys.time()
  sce_add_look <- lookup_in(sce)
  #####
  rt_scae_lookup_end <- Sys.time()
  diff_rt_scae_lookup <- rt_scae_lookup_end - rt_scae_lookup_start
  print(paste("     Generating SCAE (1/5) extending rowData:", diff_rt_scae_lookup))
  #####

  rt_scae_filt_norm_start <- Sys.time()
  sce_filter_norm <- filter_norm(sce_add_look, threshold)
  #####
  rt_scae_filt_norm_end <- Sys.time()
  diff_rt_scae_filt_norm <- rt_scae_filt_norm_end - rt_scae_filt_norm_start
  print(paste("     Generating SCAE (2/5) filtering and normalization:", diff_rt_scae_filt_norm))
  #####

  rt_scae_a2g_start <- Sys.time()
  scae <- alleles2genes(sce_filter_norm, lookup)
  #####
  rt_scae_a2g_end <- Sys.time()
  diff_rt_scae_a2g <- rt_scae_a2g_end - rt_scae_a2g_start
  print(paste("     Generating SCAE (3/5) alleles2genes:", diff_rt_scae_a2g))
  #####

  rt_scae_g2f_start <- Sys.time()
  scae <- genes2functional(scae, lookup)
  #####
  rt_scae_g2f_end <- Sys.time()
  diff_rt_scae_g2f <- rt_scae_g2f_end - rt_scae_g2f_start
  print(paste("     Generating SCAE (4/5) genes2functional:", diff_rt_scae_g2f))
  #####

  rt_scae_log_start <- Sys.time()
  scae <- log_transform(scae)
  #####
  rt_scae_log_end <- Sys.time()
  diff_rt_scae_log <- rt_scae_log_end - rt_scae_log_start
  print(paste("     Generating SCAE (5/5) log_transform:", diff_rt_scae_g2f))
  #####

  .scae(scae)
}
#####

################################################################################
#--------------------Functions used in the SCAE-Constructor--------------------#
################################################################################

#-1-------------------------------lookup_in------------------------------------#

#####
#' Extending rowData
#'
#' @description
#' Internal function used in the `SingleCellAlleleExperiment()` constructor adding information to the SingleCellAlleleExperiment object by
#' extending the rowData by two columns. "NI_I" is a classifier for each feature_row if its considered a
#' non-immune (NI) or immune (I) gene. "Quant_type" is a classifier for determinig which row is related to which
#' subassay of the extended main assay in the SingleCellAlleleExperiment. "A" corresponds to allele, "G" to allele gene and
#' ""F" to functional allele class.
#'
#' @param sce SingleCellExperiment object
#'
#' @return SingleCellAlleleExperiment object with extended rowData
lookup_in <- function(sce){
  new_sce <- sce
  allele_names_all <- find_allele_ids(new_sce)

  SummarizedExperiment::rowData(new_sce[allele_names_all,])$NI_I <- "I"
  SummarizedExperiment::rowData(new_sce[allele_names_all,])$Quant_type <- "A"

  SummarizedExperiment::rowData(new_sce)[!(rownames(new_sce) %in% allele_names_all), ]$NI_I <- "NI"
  SummarizedExperiment::rowData(new_sce)[!(rownames(new_sce) %in% allele_names_all), ]$Quant_type <- "G"
  new_sce
}
#####

#-2------------------barcode filtering and normalization-----------------------#

#####
#' Computing size factors and normalize the filtered counts
#'
#' @description
#' Internal function used in `filter_norm()`
#'
#' @param sce SingleCellExperiment object
#' @param num_of_bins number of bins for the histograms of sizeFactors; default set to 75
#'
#' @return normalized SingleCellExperiment object
sfactors <- function(sce, num_of_bins = 75){
  working_copy <- sce
  df_scales <- scuttle::computeLibraryFactors(working_copy)
  normed_counts <- scuttle::normalizeCounts(df_scales, size_factors = SingleCellExperiment::sizeFactors(df_scales), transform = "none")
  SingleCellExperiment::counts(df_scales) <- normed_counts
  df_scales
}

#' Preprocessing
#'
#' @description
#' Internal function used in `SingleCellAlleleExperiment()` constructor as a preprocessing step for
#' filtering the barcodes and normalizing the count values.
#'
#'
#' @param sce SingleCellExperiment object
#' @param threshold Counts threshold for barcodes
#'
#' @return filtered and normalized SingleCellExperiment object
filter_norm <- function(sce, threshold = 0){
  working_copy <- sce
  filtered <- working_copy[,Matrix::colSums(SingleCellExperiment::counts(working_copy)) > threshold]

  normed <- sfactors(filtered)
  normed
}
#####

#-3-----------------------------allele2genes-----------------------------------#

#####
#' Identify rows containing allele information
#'
#' @description
#' Internal function used in `get_allelecounts()` to subsample the quantification assay and only
#' return the rows specifying allele-quantification information.
#'
#' @param sce SingleCellExperiment object
#'
#' @return subsample of the scae containing all rows with allele_counts
find_allele_ids <- function(sce){
  a <- !grepl("ENS", rownames(SingleCellExperiment::counts(sce)), fixed = TRUE)
  allele_names_all <- rownames(SingleCellExperiment::counts(sce)[a,])
  allele_names_all
}

#' Internal Error handler
#'
#' @description
#' Internal function used in `get_allelecounts()` to check if an allele identifier cannot be found in the lookup table AND does not have proper nomenclature
#' form, then the execution of further step and thus the generation of an SingleCellAlleleExperiment object stopped.
#'
#' @param sce SingleCellExperiment object
#' @param find_allele_ids list containing all allele identifiers present in the raw data
#'
#' @return stops code execution if condition not met
check_unknowns <- function(sce, find_allele_ids){
  names <- find_allele_ids
  #checks if all the identifiers of find_allele_ids have a "*" (nomenclature)
  check_star <- sum(grepl("*", names, fixed = TRUE))
  check_length <- length(names)

  if (check_star != check_length){
    star <- !grepl("*", names, fixed = TRUE)
    unknown_info <- rownames(sce[names[star],])
    stop("Allele information contains unknown identifier. Please check the data and remove rows of the following allele features identifiers: `", unknown_info, "` or use proper nomenclature.")
  }
}

#' Find not yet known allele identifiers
#'
#' @description
#' Internal function used in `get_allelecounts()`to find allele identifier that arent present in the lookup table.
#'
#' @param scae SingleCellAlleleExperimentobject
#' @param agene_names list of allele gene names
#'
#' @return list of identifiers that can not be found in the allele lookup table
find_not_ident <- function(scae, agene_names){
  scae_copy <- scae
  scae_copy <- get_alleles(scae_copy)
  scae_copy_counts <- SingleCellExperiment::counts(scae_copy)
  rownames(scae_copy_counts) <- agene_names
  not_ids <- rownames(scae_copy_counts[!grepl(c("^HLA|^TCR|^KIR"), rownames(scae_copy_counts)), , drop = FALSE])
  not_ids
}

#' Build new substring
#'
#' @description
#' Internal function used in `get_allelecounts()`. Function is used to cut character string at "*" character and return it.
#'
#' @param allele_id list of unidentified allele_names that are still in a proper nomenclature form
#'
#' @return new identifier with a cut off name to be present in the list of gene_names
cutname <- function(allele_id){
  id <- allele_id
  id <- strsplit(id, "\\*")[[1]][1]
  id
}

#' Get Subassay with allele gene names and raw allele quantification
#'
#' @description
#' Internal function used to build a subassay containing counts from raw alleles
#' The rownames  of this subassay are already translated to the corresponding allele gene identifier, which
#' are extracted from the allele lookup table
#'
#' @param sce SingleCellExperiment object
#' @param lookup allele lookup table
#'
#' @importFrom Matrix sparseMatrix
#'
#' @return subsample of the sce containing all allele counts with the allele gene identifier as rownames
get_allelecounts <- function(sce, lookup){
  wor_copy <- sce

  allele_ids_lookup <- find_allele_ids(sce)
  check_unknowns(wor_copy, allele_ids_lookup)
  unknown <- FALSE

  list_alid <- list()
  for (i in 1:length(allele_ids_lookup)){
    if (allele_ids_lookup[i] %in% lookup[grepl(allele_ids_lookup[i], lookup$Allele, fixed = TRUE),]){
      new_ids <- list(lookup[grepl(allele_ids_lookup[i], lookup$Allele, fixed = TRUE),]$Gene)
      list_alid[[length(list_alid) + 1]] <- new_ids
    }else{
      print(paste(allele_ids_lookup[i], "cant be found in the lookup table"))
      new_ids <- list(cutname(allele_ids_lookup[i]))
      list_alid[[length(list_alid) + 1]] <- new_ids
      unknown <- TRUE
    }
  }
  alid_gene_names <- unlist(list_alid)

  alleletogene_counts <- SingleCellExperiment::counts(wor_copy)[allele_ids_lookup,]
  rownames(alleletogene_counts) <- alid_gene_names

  if (unknown){
    not_ids <- find_not_ident(wor_copy, alid_gene_names)
    return_unknown <- c(alleletogene_counts, not_ids)
    return(return_unknown)
  }else {
    return_known <- c(alleletogene_counts)
    return(return_known)
  }
}

#' Building first new subassay for SingleCellAllelexperiment object
#'
#' @description
#' Internal function for the first assay extension used in the `SingleCellAlleleExperiment()` constructor
#' computing the first of the two new subassays that get appended to the
#' quantification assay. This subassay contains the allele gene identifiers instead of the allelen identifiers and
#' sums up the expression counts of alleles that have the same allele gene identifiers.
#'
#' @param sce SingleCellExperiment object
#' @param lookup allele lookup table
#'
#' @importFrom Matrix sparseMatrix
# @importFrom DelayedArray DelayedArray
#'
#' @return adds allele gene-subassay containing summarized count information for the allele genes
alleles2genes <- function(sce, lookup){
  w_copy <- sce
  unknown <- FALSE

  v_acounts <- get_allelecounts(w_copy, lookup)
  if (length(v_acounts) < 2){
    alleletogene_counts <- v_acounts[1][[1]]
  }else {
    alleletogene_counts <- v_acounts[1][[1]]
    not_ids <- v_acounts[2][[1]]
    unknown <- TRUE
  }

  uniqs <- unique(rownames(alleletogene_counts))
  al_gene <- matrix(0,
                    nrow = length(uniqs),
                    ncol = ncol(alleletogene_counts))
  rownames(al_gene) <- uniqs

  for (i in 1:length(uniqs)){
    uniq_sum <- Matrix::colSums(alleletogene_counts[rownames(alleletogene_counts) %in% uniqs[i], , drop = FALSE])
    al_gene[i,] <- uniq_sum
  }

  if (unknown){
    al_gene <- al_gene[!(rownames(al_gene) %in% not_ids), , drop = FALSE]
    SummarizedExperiment::rowData(w_copy[startsWith(rownames(SingleCellExperiment::rowData(w_copy)), not_ids)])$Quant_type <- "A_unknown"
  }
  al_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = al_gene))
  SummarizedExperiment::colData(al_sce) <- SingleCellExperiment::colData(w_copy)
  SummarizedExperiment::rowData(al_sce)$ID <- rownames(al_gene)
  new_sce <- rbind(w_copy, al_sce)

  SummarizedExperiment::rowData(new_sce[rownames(new_sce) %in% uniqs])$NI_I <- "I"
  SummarizedExperiment::rowData(new_sce[rownames(new_sce) %in% uniqs])$Quant_type <- "G"

  new_sce
}
#####

#-4------------------------------genes2func------------------------------------#

#####
#' Building second new subassay for the SingleCellAlleleExperiment object
#'
#' @description
#' Internal function for the second assay extension used in the `SingleCellAlleleExperiment()` constructor
#' computing the second of the two new subassays that get appended to the
#' quantification assay. This subassay contains the functional allele classes and
#' sums up the expression counts of the allele genes that are in the same functional group.
#'
#' @param sce SingleCellExperiment object
#' @param lookup allele lookup table
#'
#' @return adds functional subassay containing summarized count information for the functional allele classes
genes2functional <- function(sce, lookup){
  func_copy <- sce

  gene_names <- rownames(get_agenes(func_copy))
  list_func  <- list()
  for (i in 1:length(gene_names)){
    func_classes <- lookup$Function[lookup$Gene %in% gene_names[i]][1]
    list_func[[length(list_func) + 1]] <- func_classes
  }
  mtx <- matrix(0, nrow=length(gene_names), ncol = 2)
  mtx[,1] <- gene_names
  mtx[,2] <- unlist(list_func)

  func_group_genes <- list()
  for (i in unique(mtx[,2])){
    func_genes <- mtx[mtx[,2] == i, 1]
    func_group_genes[[length(func_group_genes) + 1]] <- func_genes
  }


  gene_func <- matrix(0,
                      nrow = length(func_group_genes),
                      ncol = ncol(func_copy[1,]))
  rownames(gene_func) <- unique(mtx[,2])

  agene_counts <- SingleCellExperiment::counts(get_agenes(func_copy))
  #für jede klasse kombinieren wir die gencounts von allen genen, dieser klasse
  for (i in 1:length(func_group_genes)){
    gene_colsums <- Matrix::colSums(agene_counts[func_group_genes[[i]],])
    gene_func[i,] <- gene_colsums
  }

  func_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = gene_func))
  SummarizedExperiment::colData(func_sce) <- SingleCellExperiment::colData(func_copy)
  SummarizedExperiment::rowData(func_sce)$ID <- unique(mtx[,2])
  final_scae <- rbind(func_copy, func_sce)

  SummarizedExperiment::rowData(final_scae[unique(mtx[,2])])$NI_I <-  "I"
  SummarizedExperiment::rowData(final_scae[unique(mtx[,2])])$Quant_type <-  "F"

  final_scae
}
#####

#-5-------------------------log transform counts-------------------------------#

#####
#' log transform normalized counts
#'
#' @description
#' Internal function used in in
#'
#'
#' @param sce SingleCellExperiment object
#'
#' @return SingleCellAlleleExperiment object with an additional assay containing
#' log normalized counts
log_transform <- function(sce){
  working_copy <- sce
  log_counts <- log1p(SingleCellExperiment::counts(working_copy))
  SummarizedExperiment::assays(working_copy)$logcounts <- log_counts

  SingleCellExperiment::counts(working_copy) <- DelayedArray::DelayedArray(SingleCellExperiment::counts(working_copy))
  SingleCellExperiment::logcounts(working_copy) <- DelayedArray::DelayedArray(SingleCellExperiment::logcounts(working_copy))

  working_copy
}
######

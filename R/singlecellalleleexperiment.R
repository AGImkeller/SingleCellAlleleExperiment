#library(SingleCellExperiment)

#' `.scae` SingleCellAlleleExperiment-Class derived from SingleCellExperiment-class
#'
#' @import SingleCellExperiment
#' @return definition for the scae class
.scae <- setClass("SingleCellAlleleExperiment", contains = "SingleCellExperiment")


#' `SingleCellAlleleExperiment` Constructor for the scae class
#'
#' @param ... parameters to fill in the sce constructor
#' @param lookup allele lookup file
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @import tools
#'
#' @return SingleCellAlleleExperiment object
#' @export
SingleCellAlleleExperiment <- function(..., lookup){
  sce <- SingleCellExperiment(...)

  sce_add_look <- lookup_in(sce)

  rt_scae_a2g_start <- Sys.time()
  scae <- alleles2genes(sce_add_look, lookup)
  rt_scae_a2g_end <- Sys.time()
  diff_rt_scae_a2g <- rt_scae_a2g_end - rt_scae_a2g_start
  print(paste("     Generating SCAE (1/2) alleles2genes:", diff_rt_scae_a2g))

  rt_scae_g2f_start <- Sys.time()
  scae <- genes2functional(scae, lookup)
  rt_scae_g2f_end <- Sys.time()
  diff_rt_scae_g2f <- rt_scae_g2f_end - rt_scae_g2f_start
  print(paste("     Generating SCAE (2/2) genes2functional:", diff_rt_scae_g2f))

  .scae(scae)
}

#-------------------------allele2genes + helpers--------------------------------
#####
#' `find_allele_ids` internal function used to subsample the quantification assay and only
#' return the rows specifying allele-quantification information
#'
#' @param sce SingleCellExperiment object
#'
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#'
#' @return subsample of the sce containing all rows with allele_counts
find_allele_ids <- function(sce){
  allele_char = "*"
  all_genes = rownames(counts(sce))
  a <- grepl(allele_char, all_genes, fixed = TRUE)
  allele_names_all <- rownames(counts(sce)[a,])
  allele_names_all
}

#' `get_allelecounts` internal function used to build a subassay containing only the allelecounts.
#' The rownames  of this subassay are already translated to the corresponding allele_gene names, which
#' are extracted from the allele lookup table
#'
#' @param sce SingleCellExperiment object
#' @param lookup allele lookup table
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @importFrom Matrix sparseMatrix
#' @return subsample of the sce containing all allele_counts with the allele_gen_names as rownames
get_allelecounts <- function(sce, lookup){
  wor_copy <- sce

  allele_ids_lookup <- find_allele_ids(wor_copy)
  list_alid <- list()
  for (i in 1:length(allele_ids_lookup)){
    new_ids <- list(lookup[grepl(allele_ids_lookup[i], lookup$Allele, fixed = TRUE),]$Gene)
    list_alid[[length(list_alid) + 1]] <- new_ids
  }
  alid_gene_names <- unlist(list_alid)

  alleletogene_counts <- counts(wor_copy)[allele_ids_lookup,]
  alleletogene_counts <- head(alleletogene_counts, n=length(alid_gene_names))
  rownames(alleletogene_counts) <- alid_gene_names
  alleletogene_counts <- as(alleletogene_counts, "sparseMatrix")
  alleletogene_counts
}

#' `allele2genes` first internal translation function used in the SingleCellAlleleExperiment constructor
#' computing the one of the two new subassay that get appended to the
#' quantification assay. This subassay contains the allele-gene names instead of the allele-identifiers and
#' sums up the expression counts of alleles that have the same allele-gene name.
#'
#' @param sce SingleCellExperiment object
#' @param lookup allele lookup table
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import HDF5Array
#' @importFrom Matrix sparseMatrix
#' @return adds allele-gene subassay containing summarized count-information for the allele-genes to the main assay
#'
alleles2genes <- function(sce, lookup){
  w_copy <- sce
  #allele counts with allele_gene names as rownames
  alleletogene_counts <- get_allelecounts(w_copy, lookup)
  uniqs <- unique(rownames(alleletogene_counts))
  #build new matrix for the colSum
  al_gene <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                          dims = c(length(uniqs), ncol(alleletogene_counts)),
                          dimnames = list(uniqs, NULL))

  for (i in 1:length(uniqs)){
    uniq_sum <- colSums(alleletogene_counts[rownames(alleletogene_counts) %in% uniqs[i],])
    al_gene[i,] <- uniq_sum
  }
  al_gene <- DelayedArray(al_gene)

  al_sce <- SingleCellExperiment(assays = list(counts = al_gene))
  colData(al_sce) <- colData(w_copy)
  rowData(al_sce)$ID <- uniqs

  #Merge metadata and colData from the original sce object to al_sce
  new_sce <- rbind(sce, al_sce)
  new_sce <- lookup_gene(new_sce, lookup)
  new_sce

  #----test to compare the counts----
  #test_counts <- sum(alleletogene_counts[rownames(alleletogene_counts) %in% uniqs[1],])
  #test_colSum <- sum(al_gene[rownames(al_gene) %in% uniqs[1],])
  #cat(paste("Counts of", uniqs[1] ,"before colSum: ",test_counts,"\n"))
  #cat(paste("Counts of", uniqs[1], "after colSum: ",test_colSum, "\n"))
  #cat(paste("---------------------------------------", "\n"))

  #if (test_counts == test_colSum){
  #  return(new_sce)
  #}else{
  #  print("Counts and ColSum do not match", "\n")
  #}
}

#####

#--------------------------genes2func + helpers---------------------------------
#####

#' `genes2functional` second internal translation function used in the SingleCellAlleleExperiment constructor
#' computing the second of the two new subassays that get appended to the
#' quantification assay. This subassay contains the functional allele classes and
#' sums up the expression counts of the allele genes that are in the same functional group
#'
#' @param sce SingleCellExperiment object
#' @param lookup allele lookup table
#'
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import HDF5Array
#'
#' @return adds functional subassay containing summarized count information for the functional allele classes
genes2functional <- function(sce, lookup){
  func_copy <- sce
  #alle hla_names zu gene_names benennen
  gene_names <- rownames(get_agenes(func_copy))
  list_func <- list()
  for (i in 1:length(gene_names)){
    func_classes <- lookup$Function[lookup$Gene %in% gene_names[i]][1]
    list_func[[length(list_func) + 1]] <- func_classes
  }
  mtx <- matrix(0, nrow=length(gene_names), ncol = 2)
  mtx[,1] <- gene_names
  mtx[,2] <- unlist(list_func)

  #getting a list that contains the agene_names of each functional class
  func_group_genes <- list()
  for (i in unique(mtx[,2])){
    func_genes <- mtx[mtx[,2] == i, 1]
    func_group_genes[[length(func_group_genes) + 1]] <- func_genes
  }
  #create new matrix
  gene_func <- matrix(0, nrow = length(func_group_genes), ncol = length(counts(func_copy[1,])), )
  gene_func <- as(gene_func, "sparseMatrix")
  rownames(gene_func) <- unique(mtx[,2])

  agene_counts <- counts(get_agenes(func_copy))
  agene_counts <- as(agene_counts, "sparseMatrix")

  #für jede klasse kombinieren wir die gencounts von allen genen, dieser klasse
  for (i in 1:length(func_group_genes)){
    gene_colsums <- colSums(agene_counts[func_group_genes[[i]],])
    gene_func[i,] <- gene_colsums
  }
  gene_func <- DelayedArray(gene_func)

  func_sce <- SingleCellExperiment(assays = list(counts = gene_func))
  colData(func_sce) <- colData(func_copy)
  rowData(func_sce)$ID <- unique(mtx[,2])

  final_scae <- rbind(sce, func_sce)

  #updating rowData
  rowData(final_scae[unique(mtx[,2])])$NI_I <-  "I"
  rowData(final_scae[unique(mtx[,2])])$Quant_type <-  "F"

  #----testing the counts before and after colSums------
  #func_one <- sum(counts(func_copy[func_group_genes[[1]],]))
  #func_two <- sum(counts(func_copy[func_group_genes[[2]],]))

  #cat(paste("Counts of", func_group_genes[1] ,"before colSum: ", func_one,"\n"))
  #cat(paste("Counts of", func_group_genes[2] ,"before colSum: ", func_two,"\n"))

  #cat(paste("Counts of", rownames(gene_func)[1], "after colSum: ", sum(gene_func[1,]), "\n"))
  #cat(paste("Counts of", rownames(gene_func)[2], "after colSum: ", sum(gene_func[2,])))

  final_scae
}
#####
#-----------------------------common utils--------------------------------------
#####

#' `lookup_in` internal function adding information to the SingleCellAlleleExperiment object by
#' extending the rowData by two columns. "NI_I" is a classifier for each feature_row if its considered a
#' non-immune (NI) or immune (I) gene. "Quant_type" is a classifier for determinig which row is related to which
#' subassay of the extended main assay in the SingleCellAlleleExperiment
#'
#' @param sce SingleCellExperiment object
#'
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#'
#' @return SingleCellExperiment object with extended rowData
#'
lookup_in <- function(sce){
  new_sce <- sce
  allele_names_all <- find_allele_ids(new_sce)

  rowData(new_sce[allele_names_all,])$NI_I <- "I"
  rowData(new_sce[allele_names_all,])$Quant_type <- "A"

  rowData(new_sce)[!(rownames(new_sce) %in% allele_names_all), ]$NI_I <- "NI"
  rowData(new_sce)[!(rownames(new_sce) %in% allele_names_all), ]$Quant_type <- "G"
  new_sce
}

#' `lookup_gene` internal function adding information to the SingleCellAlleleExperiment object by
#' filling the rowData columns ("NI_I" and "Quant_type") with the correct information for the rows of the newly added
#' subassays after the
#'
#' @param sce singlecellexperiment object
#' @param lookup allele lookup table
#'
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#'
#' @return singlecellexperiment object
lookup_gene <- function(sce, lookup){
  new_sce_comb <- sce
  names_check <- get_allelecounts(new_sce_comb, lookup)
  uniqs_comb <- unique(rownames(names_check))

  rowData(new_sce_comb[rownames(new_sce_comb) %in% uniqs_comb])$NI_I <- "I"
  rowData(new_sce_comb[rownames(new_sce_comb) %in% uniqs_comb,])$Quant_type <- "G"
  new_sce_comb
}

#' `get_alleles` getter function which creates scae-object for the subassay containing allele information
#'
#' @param scae SingleCellAlleleExperiment object
#'
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#'
#' @return subsampled SingleCellAlleleExperiment object
#' @export
get_alleles <- function(scae){
  alleles <- scae[rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "A"]
  alleles
}

#' `get_agenes` getter function which creates scae-object for the subassay containing allele_gene information
#'
#' @param scae SingleCellAlleleExperiment object
#'
#' @import SummarizedExperiment
#'
#' @return subsampled SingleCellAlleleExperiment object
#' @export
get_agenes <- function(scae){
  agenes <- scae[rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "G"]
  agenes
}

#' `get_func` getter function which creates scae-object for the subassay containing functional class information
#'
#' @param scae SingleCellAlleleExperiment object
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#'
#' @return subsampled SingleCellAlleleExperiment object
#' @export
get_func <- function(scae){
  func <- scae[rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "F"]
  func
}

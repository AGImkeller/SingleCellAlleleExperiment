
#---------SingleCellAlleleExperiment class definition and constructor----------#

#' SingleCellAlleleExperiment class definition
#'
#' @description
#' Defining the `SingleCellAlleleExperiment` class derived from `SingleCellExperiment` class.
#'
#' @importFrom methods new
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @return definition for the scae class
.scae <- setClass("SingleCellAlleleExperiment", contains = "SingleCellExperiment")

#' Constructor SingleCellAlleleExperiment class
#'
#' @description
#' Constructor for the `SingleCellAllelExperiment` (SCAE) class.
#' Constructor is used in the read in function `read_allele_counts()`. Performing all necessary steps to transform
#' a `SingleCellExperiment` object into the extended `SingleCellAlleleExperiment` object. SCAE objects
#' contain allele, gene and functional level quantification results. The additional layers are stored as additional
#' rows in the count assays as well as in extended rowData.
#'
#' @param ... Arguments passed to the \code{\link{SingleCellExperiment}} constructor to fill the slots of the SCE-class.
#' @param lookup A data.frame object containing the lookup table.
#' @param threshold An integer value used as a threshold for filtering low-quality barcodes/cells.
#' @param exp_type A vector containing two character strings. Either `"WTA"` or `"Amplicon"` are valid inputs. Choose one depending on the used transcriptomics approach.
#' @param verbose A logical parameter to decide if runtime-messages should be shown during function execution.
#'  Use `FALSE` if no info runtime-messages should be shown (default), and `TRUE` for showing runtime-messages.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @return A SingleCellAlleleExperiment object.
SingleCellAlleleExperiment <- function(..., threshold, exp_type, lookup, verbose = FALSE){
  sce <- SingleCellExperiment(...)

  rt_scae_lookup_start <- Sys.time()
  sce_add_look <- ext_rd(sce, exp_type, verbose = verbose)

  if (verbose){
    rt_scae_lookup_end <- Sys.time()
    diff_rt_scae_lookup <- round(rt_scae_lookup_end - rt_scae_lookup_start, digits = 2)
    message("     Generating SCAE (1/5) extending rowData: ", diff_rt_scae_lookup, " seconds")
  }


  rt_scae_filt_norm_start <- Sys.time()
  sce_filter_norm <- filter_norm(sce_add_look, threshold)
  if (verbose){
    rt_scae_filt_norm_end <- Sys.time()
    diff_rt_scae_filt_norm <- round(rt_scae_filt_norm_end - rt_scae_filt_norm_start, digits = 2)
    message("     Generating SCAE (2/5) filtering and normalization: ", diff_rt_scae_filt_norm, " seconds")
  }


  rt_scae_a2g_start <- Sys.time()
  scae <- alleles2genes(sce_filter_norm, lookup, exp_type)
  if (verbose){
    rt_scae_a2g_end <- Sys.time()
    diff_rt_scae_a2g <- round(rt_scae_a2g_end - rt_scae_a2g_start, digits = 2)
    message("     Generating SCAE (3/5) alleles2genes: ", diff_rt_scae_a2g, " seconds")
  }


  rt_scae_g2f_start <- Sys.time()
  scae <- genes2functional(scae, lookup, exp_type)
  if (verbose){
    rt_scae_g2f_end <- Sys.time()
    diff_rt_scae_g2f <- round(rt_scae_g2f_end - rt_scae_g2f_start, digits = 2)
    message("     Generating SCAE (4/5) genes2functional: ", diff_rt_scae_g2f, " seconds")
  }


  rt_scae_log_start <- Sys.time()
  scae <- log_transform(scae)
  if (verbose){
    rt_scae_log_end <- Sys.time()
    diff_rt_scae_log <- round(rt_scae_log_end - rt_scae_log_start, digits = 2)
    message("     Generating SCAE (5/5) log_transform: ", diff_rt_scae_g2f, " seconds")
  }

  .scae(scae)
}


#--------------------Functions used in the SCAE-Constructor--------------------#

#-1--------------------------------ext_rd--------------------------------------#

#' Extending rowData
#'
#' @description
#' Internal function used in the `SingleCellAlleleExperiment()` constructor adding information to the SingleCellAlleleExperiment object by
#' extending the rowData by two columns. `NI_I` is a classifier for each feature_row if its considered a
#' non-immune (NI) or immune (I) gene. `Quant_type` is a classifier for determining which row is related to which
#' subassay of the extended main assay in the `SingleCellAlleleExperiment`. "A" corresponds to allele, "G" to allele gene and
#' "F" to functional allele class.
#'
#' @param sce A \code{\link{SingleCellExperiment}} object. Object is initally constructed in the `SingleCellAlleleExperiment` constructor.
#' @param exp_type A vector containing two character strings. Either `"WTA"` or `"Amplicon"` are valid inputs. Choose one depending on the used transcriptomics approach.
#' @param verbose A logical parameter to decide if runtime-messages should be shown during function execution.
#'  Use `FALSE` if no info runtime-messages should be shown (default), and `TRUE` for showing runtime-messages.
#'
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom SingleCellExperiment rowData
#'
#' @return A SingleCellExperiment object.
ext_rd <- function(sce, exp_type, verbose = FALSE){
  if (exp_type == "ENS"){
    if (verbose){
      message("Using org.Hs to retrieve NCBI gene identifiers.")
    }
    gene_symbols <- get_ncbi_org(sce)
    rowData(sce)$Symbol <- gene_symbols
  }

  allele_names_all <- find_allele_ids(sce)

  # Group of genes for which extended informaton is stored
  rowData(sce[allele_names_all,])$NI_I <- "I"

  # Allele level
  rowData(sce[allele_names_all,])$Quant_type <- "A"

  # Group of genes for which classical (gene level) informaton is stored
  rowData(sce)[!(rownames(sce) %in% allele_names_all), ]$NI_I <- "NI"

  # Gene level
  rowData(sce)[!(rownames(sce) %in% allele_names_all), ]$Quant_type <- "G"

  rowData(sce)[rownames(rowData(scae_subset_alleles(sce))),]$Symbol <- rownames(rowData(scae_subset_alleles(sce)))

  sce
}

# Code provided by Ahmad Al Ajami
#' Get NCBI genes using the org.HS.db package
#'
#' @description
#' This internal function is not as accurate (does not retrieve as many gene names as `biomaRt`) but can be used without
#' internet connection.
#'
#' @param sce A \code{\link{SingleCellExperiment}} object.
#'
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL org.Hs.egENSEMBL
#' @importFrom AnnotationDbi mappedkeys
#' @importFrom methods as
#' @importFrom SingleCellExperiment rowData
#'
#' @return A list of character strings for gene names.
get_ncbi_org <- function(sce){
  ensembl_ids <- rowData(sce)$Ensembl_ID
  ensembl_ids <- sub("\\..*", "", ensembl_ids)

  Hs_symbol  <- org.Hs.eg.db::org.Hs.egSYMBOL
  Hs_ensembl <- org.Hs.eg.db::org.Hs.egENSEMBL
  mapped_Hs_genes_symbol  <- AnnotationDbi::mappedkeys(Hs_symbol)
  mapped_Hs_genes_ensembl <- AnnotationDbi::mappedkeys(Hs_ensembl)
  Hs_symbol_df  <- as.data.frame(Hs_symbol[mapped_Hs_genes_symbol])
  Hs_ensembl_df <- as.data.frame(Hs_ensembl[mapped_Hs_genes_ensembl])

  Hs_mapping <- merge(Hs_symbol_df, Hs_ensembl_df)

  indic <- match(ensembl_ids, Hs_mapping$ensembl_id)
  ncbi_symbols <- Hs_mapping$symbol[match(ensembl_ids, Hs_mapping$ensembl_id)]

  return(ncbi_symbols)
}


#-2------------------barcode filtering and normalization-----------------------#

#' Preprocessing
#'
#' @description
#' Internal function used in `SingleCellAlleleExperiment()` constructor as a preprocessing step for
#' filtering the barcodes and normalizing the count values.
#'
#' @param sce A \code{\link{SingleCellExperiment}} object.
#' @param threshold An integer value used as a threshold for filtering low-quality barcodes/cells.
#'
#' @importFrom Matrix colSums
#' @importFrom SingleCellExperiment counts
#' @importFrom scuttle computeLibraryFactors
#'
#' @return A SingleCellExperiment object.
filter_norm <- function(sce, threshold = 0){

  filtered  <- sce[, colSums(counts(sce)) > threshold]
  df_scales <- computeLibraryFactors(filtered)
  df_scales
}


#-3-----------------------------allele2genes-----------------------------------#

#' Identify rows containing allele information for WTA
#'
#' @description
#' Internal function used in `get_allelecounts()` to subsample the quantification assay and only
#' return the rows specifying allele-quantification information.
#'
#' @param sce A \code{\link{SingleCellExperiment}} object.
#' @importFrom SingleCellExperiment counts
#'
#' @return A SingleCellExperiment object.
find_allele_ids <- function(sce){
  a <- grepl("*", rownames(counts(sce)), fixed = TRUE)
  if (sum(a) == 0){
    a <- grepl("HLA-", rownames(counts(sce)), fixed = TRUE)
  }

  allele_names_all <- rownames(counts(sce)[a,])
  allele_names_all
}

#' Get Subassay with allele gene names and raw allele quantification
#'
#' @description
#' Internal function used to build a subassay containing counts from raw alleles
#' The rownames  of this subassay are already translated to the corresponding allele gene identifier, which
#' are extracted from the allele lookup table
#'
#' @param sce A \code{\link{SingleCellExperiment}} object.
#' @param lookup A data.frame object containing the lookup table.
#'
#' @importFrom SingleCellExperiment counts
#'
#' @return A SingleCellExperiment object.
get_allelecounts <- function(sce, lookup){

  allele_ids_lookup <- find_allele_ids(sce)
  list_alid <- list()

  for (i in seq_along(allele_ids_lookup)){
    new_ids <- list(lookup[grepl(allele_ids_lookup[i], lookup$Allele, fixed = TRUE),]$Gene)
    list_alid[[length(list_alid) + 1]] <- new_ids
  }
  alid_gene_names <- unlist(list_alid)

  alleletogene_counts <- counts(scae_subset_alleles(sce))
  rownames(alleletogene_counts) <- alid_gene_names

  return_known <- c(alleletogene_counts)
  return(return_known)
}

#' Building first new subassay for SingleCellAllelexperiment object
#'
#' @description
#' Internal function for the first assay extension used in the `SingleCellAlleleExperiment()` constructor
#' computing the first of the two new subassays that get appended to the
#' quantification assay. This subassay contains the allele gene identifiers instead of the allelen identifiers and
#' sums up the expression counts of alleles that have the same allele gene identifiers.
#'
#' @param sce A \code{\link{SingleCellExperiment}} object.
#' @param lookup A data.frame object containing the lookup table.
#' @param exp_type A character string determining wether the gene symbols in the input data are Ensemble identifiers or ncbi identifiers. Only used internally, not related to input done by the user.
#'
#' @importFrom Matrix colSums
#' @importFrom SummarizedExperiment rowData<- colData<-
#' @importFrom SingleCellExperiment rowData colData SingleCellExperiment
#' @importFrom BiocGenerics rbind
#'
#' @return A SingleCellExperiment object.
alleles2genes <- function(sce, lookup, exp_type){

  v_acounts <- get_allelecounts(sce, lookup)

  alleletogene_counts <- v_acounts[1][[1]]

  uniqs   <- unique(rownames(alleletogene_counts))
  al_gene <- matrix(0, nrow = length(uniqs), ncol = ncol(alleletogene_counts))
  rownames(al_gene) <- uniqs

  for (i in seq_along(uniqs)){
    uniq_sum <- colSums(alleletogene_counts[rownames(alleletogene_counts) %in% uniqs[i], , drop = FALSE])
    al_gene[i,] <- uniq_sum
  }

  al_sce <- SingleCellExperiment(assays = list(counts = al_gene),
                                 colData = colData(sce))
  rowData(al_sce)$Symbol <- rownames(al_gene)


  if (exp_type == "ENS"){
    rowData(al_sce)$Ensembl_ID <- rownames(al_gene)
  }

  new_sce <- BiocGenerics::rbind(sce, al_sce)

  rowData(new_sce[rownames(new_sce) %in% uniqs])$NI_I <- "I"
  rowData(new_sce[rownames(new_sce) %in% uniqs])$Quant_type <- "G"

  new_sce
}


#-4------------------------------genes2func------------------------------------#

#' Building second new subassay for the SingleCellAlleleExperiment object
#'
#' @description
#' Internal function for the second assay extension used in the `SingleCellAlleleExperiment()` constructor
#' computing the second of the two new subassays that get appended to the
#' quantification assay. This subassay contains the functional allele classes and
#' sums up the expression counts of the allele genes that are in the same functional group.
#'
#' @param sce A \code{\link{SingleCellExperiment}} object.
#' @param lookup A data.frame object containing the lookup table.
#' @param exp_type A character string determining wether the gene symbols in the input data are Ensemble identifiers or ncbi identifiers. Only used internally, not related to input done by the user.
#'
#' @importFrom SingleCellExperiment colData counts SingleCellExperiment
#' @importFrom SummarizedExperiment colData<- rowData<-
#' @importFrom Matrix colSums
#' @importFrom BiocGenerics rbind
#'
#' @return A SingleCellExperiment object.
genes2functional <- function(sce, lookup, exp_type){

  #find functional classes for each gene
  gene_names <- rownames(get_agenes(sce))
  list_func  <- list()
  for (i in seq_along(gene_names)){
    func_classes <- lookup$Function[lookup$Gene %in% gene_names[i]][1]
    list_func[[length(list_func) + 1]] <- func_classes
  }
  gene_func_names <- unlist(list_func)

  genetofunc_counts <- counts(get_agenes(sce))
  rownames(genetofunc_counts) <- gene_func_names

  uniqs     <- unique(rownames(genetofunc_counts))
  gene_func <- matrix(0,
                      nrow = length(uniqs),
                      ncol = ncol(sce[1,]))
  rownames(gene_func) <- uniqs

  for (i in seq_along(uniqs)){
    gene_colsums  <- colSums(genetofunc_counts[rownames(genetofunc_counts) %in% uniqs[i], , drop = FALSE])
    gene_func[i,] <- gene_colsums
  }

  func_sce <- SingleCellExperiment(assays = list(counts = gene_func),
                                   colData = colData(sce))
  rowData(func_sce)$Symbol <- rownames(func_sce)

  if (exp_type == "ENS"){
    rowData(func_sce)$Ensembl_ID <- rownames(func_sce)
  }

  final_scae <- BiocGenerics::rbind(sce, func_sce)

  # Genes with extended quantification
  rowData(final_scae[rownames(final_scae) %in% uniqs])$NI_I <- "I"
  # Functional level
  rowData(final_scae[rownames(final_scae) %in% uniqs])$Quant_type <- "F"

  final_scae
}


#-5-------------------------log transform counts-------------------------------#

#' Log-transform normalized counts
#'
#' @description
#' Internal function used in the `SingleCellAlleleExperiment()` constructor to log-normalize the raw counts and add them to the \code{\link{logcounts}} assay.
#'
#' @param sce A \code{\link{SingleCellExperiment}} object.
#'
#' @importFrom scuttle normalizeCounts
#' @importFrom SingleCellExperiment sizeFactors counts logcounts counts<- logcounts<-
#' @importFrom SummarizedExperiment assays<- assays
#' @importFrom DelayedArray DelayedArray
#'
#' @return A SingleCellExperiment object.
log_transform <- function(sce){

  normed_counts <- normalizeCounts(sce,
                                   size_factors = sizeFactors(sce),
                                   transform = "log")

  assays(sce)$logcounts  <- normed_counts
  counts(sce)    <- DelayedArray(counts(sce))
  logcounts(sce) <- DelayedArray(logcounts(sce))

  sce
}

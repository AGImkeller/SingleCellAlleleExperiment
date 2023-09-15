################################################################################
##############---Class definition of SingleCellAlleleExperiment---##############
############---file contains all helpers to transform sce to scae---############
################################################################################

#---------SingleCellAlleleExperiment class definition and constructor----------#

#####
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
#' Constructor is used in the read in function `readAlleleCounts()`. Performing all necessary steps to transform
#' a `SingleCellExperiment` object into the extended `SingleCellAlleleExperiment` object. SCAE objects
#' contain allele, gene and functional level quantification results. The additional layers are stored as additional
#' rows in the count assays as well as in extended rowData.
#'
#' @param ... parameters to pass to SingleCellExperiment constructor
#' @param lookup allele lookup file
#' @param threshold count threshold for filtering barcodes/cells
#' @param exp_type either `"WTA"` or `"Amplicon"` depending on the used experiments technology
#' @param symbols identifier used to choose which database-function to use to retrieve the ncbi gene names
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @return SingleCellAlleleExperiment object
#' @export
SingleCellAlleleExperiment <- function(..., threshold, exp_type, symbols, lookup){
  sce <- SingleCellExperiment(...)

  rt_scae_lookup_start <- Sys.time()
  sce_add_look <- ext_rd(sce, exp_type, symbols)
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
  scae <- alleles2genes(sce_filter_norm, lookup, exp_type)
  #####
  rt_scae_a2g_end <- Sys.time()
  diff_rt_scae_a2g <- rt_scae_a2g_end - rt_scae_a2g_start
  print(paste("     Generating SCAE (3/5) alleles2genes:", diff_rt_scae_a2g))
  #####

  rt_scae_g2f_start <- Sys.time()
  scae <- genes2functional(scae, lookup, exp_type)
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

#-1--------------------------------ext_rd--------------------------------------#

#####
#' Extending rowData
#'
#' @description
#' Internal function used in the `SingleCellAlleleExperiment()` constructor adding information to the SingleCellAlleleExperiment object by
#' extending the rowData by two columns. `NI_I` is a classifier for each feature_row if its considered a
#' non-immune (NI) or immune (I) gene. `Quant_type` is a classifier for determining which row is related to which
#' subassay of the extended main assay in the `SingleCellAlleleExperiment`. "A" corresponds to allele, "G" to allele gene and
#' "F" to functional allele class.
#'
#' @param sce SingleCellExperiment object
#' @param exp_type either `"WTA"` or `"Amplicon"` depending on the used experiments technology
#' @param symbols identifier used to choose which database-function to use to retrieve the NCBI gene names
#'
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom SingleCellExperiment rowData
#'
#' @return SingleCellAlleleExperiment object with extended rowData
ext_rd <- function(sce, exp_type, symbols){

  if (exp_type == "WTA"){
    if (symbols == "biomart"){
      ensembl_ids  <- unlist(rowData(sce)$Ensembl.ID)
      gene_symbols <- get_ncbi_gene_names(ensembl_ids)
    }
    if (symbols == "orgdb"){
      gene_symbols <- get_ncbi_org(sce)
    }
    rowData(sce)$Symbol <- gene_symbols
  }
  #extend rowData with new classification columns
  allele_names_all <- find_allele_ids(sce, exp_type)

  # Group of genes for which extended informaton is stored
  rowData(sce[allele_names_all,])$NI_I <- "I"
  # Allele level
  rowData(sce[allele_names_all,])$Quant_type <- "A"

  # Group of genes for which classical (gene level) informaton is stored
  rowData(sce)[!(rownames(sce) %in% allele_names_all), ]$NI_I <- "NI"
  # Gene level
  rowData(sce)[!(rownames(sce) %in% allele_names_all), ]$Quant_type <- "G"

  rowData(sce)[rownames(rowData(get_alleles(sce))),]$Symbol <- rownames(rowData(get_alleles(sce)))

  sce
}

#' Get Ncbi genes using biomaRt
#'
#' @description
#' This internal function is used to retrieve the gene-symbol names to the corresponding ENSG accession numbers in the WTA experiment approach.
#' Internet connection is mandatory, as its retrieving the newest possible dataset every time. If you have to work offline then use the
#' get_ncbi_org by specifying `symbols = "orgdb"` in the corresponding `symbols` parameter of the `readAlleleCounts()` function.
#'
#' @param ensembl_ids vector containing ensembl.ids
#'
#' @importFrom biomaRt useMart getBM
#'
#' @return vector containing the ncbi names matching the right rows for the ensembl.ids
get_ncbi_gene_names <- function(ensembl_ids) {
  ensembl_ids <- sub("\\..*", "", ensembl_ids)
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")

  results <- biomaRt::getBM(
    attributes = attributes,
    filters = "ensembl_gene_id",
    values  = ensembl_ids,
    mart    = ensembl
  )

  ncbi_gene_names <- rep(NA_character_, length(ensembl_ids))
  matching_indices <- match(results$ensembl_gene_id, ensembl_ids)
  ncbi_gene_names[matching_indices[!is.na(matching_indices)]] <- results$external_gene_name

  ncbi_gene_names[ncbi_gene_names == ""] <- NA
  return(ncbi_gene_names)
}

#' Get NCBI genes using the org.HS.db package
#'
#' @description
#' This internal function is not as accurate (does not retrieve as many ncbi gene names as `biomaRt`) but can be used without
#' internet connection.
#'
#' @param scae SingleCellAlleleExperiment object
#'
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL org.Hs.egENSEMBL
#' @importFrom AnnotationDbi mappedkeys
#' @importFrom methods as
#' @importFrom SingleCellExperiment rowData
#'
#' @return list of gene names
get_ncbi_org <- function(scae){
  ensembl_ids <- rowData(scae)$Ensembl.ID
  ensembl_ids <- sub("\\..*", "", ensembl_ids)

  Hs_symbol  <- org.Hs.eg.db::org.Hs.egSYMBOL
  Hs_ensembl <- org.Hs.eg.db::org.Hs.egENSEMBL
  mapped_Hs_genes.symbol  <- AnnotationDbi::mappedkeys(Hs_symbol)
  mapped_Hs_genes.ensembl <- AnnotationDbi::mappedkeys(Hs_ensembl)
  Hs_symbol.df  <- as.data.frame(Hs_symbol[mapped_Hs_genes.symbol])
  Hs_ensembl.df <- as.data.frame(Hs_ensembl[mapped_Hs_genes.ensembl])

  Hs_mapping <- merge(Hs_symbol.df, Hs_ensembl.df)

  indic <- match(ensembl_ids, Hs_mapping$ensembl_id)
  ncbi_symbols <- Hs_mapping$symbol[match(ensembl_ids, Hs_mapping$ensembl_id)]

  return(ncbi_symbols)
}
#####

#-2------------------barcode filtering and normalization-----------------------#

#####

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
#' @importFrom Matrix colSums
#' @importFrom SingleCellExperiment counts
#' @importFrom scuttle computeLibraryFactors
#'
#' @return filtered and normalized SingleCellExperiment object
filter_norm <- function(sce, threshold = 0){
  
  filtered  <- sce[, colSums(counts(sce)) > threshold]
  df_scales <- computeLibraryFactors(filtered)
  df_scales
}
#####

#-3-----------------------------allele2genes-----------------------------------#

#####
#' Identify rows containing allele information for WTA
#'
#' @description
#' Internal function used in `get_allelecounts()` to subsample the quantification assay and only
#' return the rows specifying allele-quantification information.
#'
#' @param sce SingleCellExperiment object
#' @param exp_type either `"WTA"` or `"Amplicon"` depending on the used experiments technology
#'
#' @importFrom SingleCellExperiment counts
#'
#' @return subsample of the scae containing all rows with allele_counts
find_allele_ids <- function(sce, exp_type){
  a <- switch(exp_type,
              "WTA"      = !grepl("ENS", rownames(counts(sce)), fixed = TRUE),
              # This needs to be fixed on the long run, because not all genes with extended information
              # match HLA
              "Amplicon" =  grepl("HLA-", rownames(counts(sce)), fixed = TRUE),
              NA)
  allele_names_all <- rownames(counts(sce)[a,])
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
  check_star   <- sum(grepl("*", names, fixed = TRUE))
  check_length <- length(names)

  if (check_star != check_length){
    star <- !grepl("*", names, fixed = TRUE)
    unknown_info <- rownames(sce[names[star],])
    stop("Allele information contains unknown identifier.
         Please check the data and remove rows of the following
         allele features identifiers: `",unknown_info, " ` or use proper nomenclature.")
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
#' @importFrom SingleCellExperiment counts
#'
#' @return list of identifiers that can not be found in the allele lookup table
find_not_ident <- function(scae, agene_names){
  #return allele genes that dont start with HLA (not found in lookup table)
  scae_counts <- counts(get_alleles(scae))
  rownames(scae_counts) <- agene_names
  # This needs to be fixed on the long run, because not all genes with extended information
  # match HLA
  not_ids <- rownames(scae_counts[!grepl(c("^HLA"), rownames(scae_counts)), , drop = FALSE])

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
  # does this work for all alleles?
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
#' @param exp_type either `"WTA"` or `"Amplicon"` depending on the used experiments technology
#'
#' @importFrom SingleCellExperiment counts
#'
#' @return subsample of the sce containing all allele counts with the allele gene identifier as rownames
get_allelecounts <- function(sce, lookup, exp_type){

  allele_ids_lookup <- find_allele_ids(sce, exp_type)
  if (exp_type == "WTA"){
    check_unknowns(sce, allele_ids_lookup)
  }
  unknown <- FALSE

  list_alid <- list()
  for (i in seq_along(allele_ids_lookup)){
    if (allele_ids_lookup[i] %in% lookup[grepl(allele_ids_lookup[i], lookup$Allele, fixed = TRUE),]){
      new_ids <- list(lookup[grepl(allele_ids_lookup[i], lookup$Allele, fixed = TRUE),]$Gene)
      list_alid[[length(list_alid) + 1]] <- new_ids
    }else{
      print(paste(allele_ids_lookup[i], "can't be found in the lookup table"))
      new_ids <- list(cutname(allele_ids_lookup[i]))
      list_alid[[length(list_alid) + 1]] <- new_ids
      unknown <- TRUE
    }
  }
  alid_gene_names <- unlist(list_alid)

  alleletogene_counts <- counts(get_alleles(sce))
  rownames(alleletogene_counts) <- alid_gene_names

  if (unknown){
    not_ids <- find_not_ident(sce, alid_gene_names)
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
#' @param exp_type either `"WTA"` or `"Amplicon"` depending on the used experiments technology
#'
#' @importFrom Matrix colSums
#' @importFrom SummarizedExperiment rowData<- colData<-
#' @importFrom SingleCellExperiment rowData colData SingleCellExperiment
#' @importFrom BiocGenerics rbind
#'
#' @return adds allele gene-subassay containing summarized count information for the allele genes
alleles2genes <- function(sce, lookup, exp_type){
  unknown <- FALSE

  v_acounts <- get_allelecounts(sce, lookup, exp_type)

  if (length(v_acounts) < 2){
    alleletogene_counts <- v_acounts[1][[1]]
  }else {
    alleletogene_counts <- v_acounts[1][[1]]
    not_ids <- unlist(v_acounts[2:length(v_acounts)])
    unknown <- TRUE
  }

  uniqs   <- unique(rownames(alleletogene_counts))
  al_gene <- matrix(0, nrow = length(uniqs), ncol = ncol(alleletogene_counts))
  rownames(al_gene) <- uniqs

  for (i in seq_along(uniqs)){
    uniq_sum <- colSums(alleletogene_counts[rownames(alleletogene_counts) %in% uniqs[i], , drop = FALSE])
    al_gene[i,] <- uniq_sum
  }

  if (unknown){
    al_gene <- al_gene[!(rownames(al_gene) %in% not_ids), , drop = FALSE]
    filtered_rows <- sapply(rownames(rowData(sce)), function(rowname){
      any(startsWith(rowname, not_ids))
    })
    rowData(sce)[filtered_rows, "Quant_type"] <- "A_unknown"
  }

  al_sce <- SingleCellExperiment(assays = list(counts = al_gene),
                                 colData = colData(sce))
  rowData(al_sce)$Symbol <- rownames(al_gene)
  if (exp_type == "WTA"){
    rowData(al_sce)$Ensembl.ID <- rownames(al_gene)
  }

  new_sce <- BiocGenerics::rbind(sce, al_sce)

  rowData(new_sce[rownames(new_sce) %in% uniqs])$NI_I <- "I"
  rowData(new_sce[rownames(new_sce) %in% uniqs])$Quant_type <- "G"

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
#' @param exp_type either `"WTA"` or `"Amplicon"` depending on the used experiments technology
#'
#' @importFrom SingleCellExperiment colData counts SingleCellExperiment
#' @importFrom SummarizedExperiment colData<- rowData<-
#' @importFrom Matrix colSums
#' @importFrom BiocGenerics rbind
#'
#' @return adds functional subassay containing summarized count information for the functional allele classes
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
  if (exp_type == "WTA"){
    rowData(func_sce)$Ensembl.ID <- rownames(func_sce)
  }

  final_scae <- BiocGenerics::rbind(sce, func_sce)
  
  # Genes with extended quantification
  rowData(final_scae[rownames(final_scae) %in% uniqs])$NI_I <- "I"
  # Functional level
  rowData(final_scae[rownames(final_scae) %in% uniqs])$Quant_type <- "F"

  final_scae
}
#####

#-5-------------------------log transform counts-------------------------------#

#####
#' Log-transform normalized counts
#'
#' @description
#' Internal function used in the `SingleCellAlleleExperiment()` constructor to log-normalize the raw counts and add them to `logcounts` assay.
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom scuttle normalizeCounts
#' @importFrom SingleCellExperiment sizeFactors counts logcounts counts<- logcounts<-
#' @importFrom SummarizedExperiment assays<- assays
#' @importFrom DelayedArray DelayedArray
#'
#' @return SingleCellAlleleExperiment object with an additional assay containing
log_transform <- function(sce){

  normed_counts <- normalizeCounts(sce,
                                   size_factors = sizeFactors(sce),
                                   transform = "log")

  assays(sce)$logcounts   <- normed_counts

  counts(sce)    <- DelayedArray(counts(sce))
  logcounts(sce) <- DelayedArray(logcounts(sce))

  sce
}
######

#-6-------------------------add sample tags------------------------------------#

#####
#' Adding sample tag information to colData
#'
#' @description
#' Internal function used in 'readAlleleCounts()'. Stated here because its supposed to be a transformation step of the SCAE object.
#' Adding sample tag information to colData
#'
#'
#'
#' @param path character string input containing the path to the directory containing the
#'   input files
#' @param scae SingleCellAlleleExperiment object
#'
#' @importFrom methods as
#' @importFrom utils read.table
#' @importFrom Matrix readMM rowSums
#' @importFrom MatrixGenerics rowMaxs
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SingleCellExperiment colData
#'
#' @return updated SingleCellAlleleExperiment object
add_sample_tags <- function(path, scae, tag_feature_mtx, tag_feature_barcodes){

  tags  <- Matrix::readMM(paste0(path, tag_feature_mtx, sep = ""))
  cells <- utils::read.table(paste0(path, tag_feature_barcodes, sep = ""), header = FALSE)
  rownames(tags) <- cells$V1
  colnames(tags) <- paste("ST_", seq_len(12), sep = "")
  tags <- methods::as(tags, "CsparseMatrix")
  dom.tags <- data.frame(Cellular.Barcode = rownames(tags),
                         Sample.Tag = ifelse(rowMaxs(tags) >= 0.75 * rowSums(tags),
                                             colnames(tags)[max.col(tags)],
                                             "Multiplet"))
  dom.tags <- dom.tags[!dom.tags$Sample.Tag == 'Multiplet', ]

  inter <- intersect(rownames(colData(scae)), dom.tags$Cellular.Barcode)
  colData(scae)$sample.tags <- NA
  colData(scae[, inter])$sample.tags <- dom.tags[inter, "Sample.Tag"]
  scae <- scae[, !is.na(colData(scae)$sample.tags)]
  scae
}
#####

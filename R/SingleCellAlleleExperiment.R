
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
#' Constructor is used in the read in function `readAlleleCounts()`. Performing all necessary steps to transform
#' a `SingleCellExperiment` object into the extended `SingleCellAlleleExperiment` object. SCAE objects
#' contain allele, gene and functional level quantification results. The additional layers are stored as additional
#' rows in the count assays as well as in extended rowData.
#'
#' @param ... Arguments passed to the \code{\link{SingleCellExperiment}} constructor to fill the slots of the SCE-class.
#' @param lookup A data.frame object containing the lookup table.
#' @param threshold An integer value used as a threshold for filtering low-quality barcodes/cells.
#' @param exp_type A vector containing two character strings. Either `"WTA"` or `"Amplicon"` are valid inputs. Choose one depending on the used transcriptomics approach.
#' @param symbols A character string used to determine which database-funtion to use to retrieve NCBI gene names. The value `"orgdb"` uses the \code{\link{org.Hs.eg.db}} package.
#' The value `"biomart"` uses the `biomaRt` package. Standard value is set to `NULL` and is updated to `"biomaRt"` during runtime if not specified.
#' @param verbose A logical parameter to decide if runtime-messages should be shown during function execution.
#'  Use `FALSE` if no info runtime-messages should be shown (default), and `TRUE` for showing runtime-messages.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @return A SingleCellAlleleExperiment object.
SingleCellAlleleExperiment <- function(..., threshold, exp_type, symbols, lookup, verbose = FALSE){
  sce <- SingleCellExperiment(...)

  rt_scae_lookup_start <- Sys.time()
  sce_add_look <- ext_rd(sce, exp_type, symbols, verbose = verbose)
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
#' @param symbols A character string used to determine which database-funtion to use to retrieve NCBI gene names. The value `"orgdb"` uses the \code{\link{org.Hs.eg.db}} package.
#' The value `"biomart"` [biomaRt](\code{\link{biomaRt}}) package. Standard value is set to `NULL` and is updated to `"biomaRt"` during runtime if not specified.
#' @param verbose A logical parameter to decide if runtime-messages should be shown during function execution.
#'  Use `FALSE` if no info runtime-messages should be shown (default), and `TRUE` for showing runtime-messages.
#'
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom SingleCellExperiment rowData
#'
#' @return A SingleCellExperiment object.
ext_rd <- function(sce, exp_type, symbols, verbose = FALSE){
  test <- symbols

  if (exp_type == "WTA"){
    if (symbols == "biomart"){
      test <- tryCatch({
        gene_symbols <- get_ncbi_gene_names(sce)
        if (verbose){
          message("Using biomart to retrieve NCBI gene identifiers.")
        }
        test <- "biomart"
      }, message = function(e){
        if (grepl("unavailable", e$message)) {
          message("Ensembl service is currently unavailable, using org.Hs.db instead")
          symbols <- "orgdb"
        } else if (grepl("unresponsive", e$message)){
          message("Ensembl service is currently unavailable, using org.Hs.db instead")
          symbols <- "orgdb"
        } else {
          symbols <- "biomart"
        }
        return(symbols)
      })
    }
    symbols = test

    if (symbols == "orgdb"){
      gene_symbols <- get_ncbi_org(sce)
      if (verbose){
        message("Using org.Hs to retrieve NCBI gene identifiers.")
      }
    }

    rowData(sce)$Symbol <- gene_symbols

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
  }
  sce
}



#' Get NCBI genes using biomaRt
#'
#' @description
#' This internal function is used to retrieve the gene-symbol names to the corresponding ENSG accession numbers in the WTA experiment approach.
#' Internet connection is mandatory, as its retrieving the newest possible dataset every time. If you have to work offline then use the
#' get_ncbi_org by specifying `symbols = "orgdb"` in the corresponding `symbols` parameter of the `readAlleleCounts()` function.
#'
#' @param sce A \code{\link{SingleCellExperiment}} object.
#'
#' @importFrom biomaRt useMart getBM
#'
#' @return A vector containing character strings for NCBI gene names.
get_ncbi_gene_names <- function(sce) {

  ensembl_ids_sce  <- unlist(rowData(sce)$Ensembl_ID)

  ensembl_ids <- sub("\\..*", "", ensembl_ids_sce)

  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")

  results <- biomaRt::getBM(
    attributes <-  attributes,
    filters <- "ensembl_gene_id",
    values  <- ensembl_ids,
    mart    <- ensembl
  )
  ncbi_gene_names  <- rep(NA_character_, length(ensembl_ids))
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
#' @param exp_type symbols A character string used to determine which database-funtion to use to retrieve NCBI gene names. The value `"orgdb"` uses the \code{\link{org.Hs.eg.db}} package.
#' The value `"biomart"` [biomaRt](\code{\link{biomaRt}}) package. Standard value is set to `NULL` and is updated to `"biomaRt"` during runtime if not specified.
#'
#' @importFrom SingleCellExperiment counts
#'
#' @return A SingleCellExperiment object.
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
#' @param sce A \code{\link{SingleCellExperiment}} object.
#' @param find_allele_ids A list containing allele identifiers present in the raw data.
#'
#' @return Error message if condition is not met.
check_unknowns <- function(sce, find_allele_ids){
  names <- find_allele_ids
  # checks if all the identifiers of find_allele_ids have a "*" (nomenclature)
  check_star   <- sum(grepl("*", names, fixed = TRUE))
  check_length <- length(names)

  if (check_star != check_length){
    star <- !grepl("*", names, fixed = TRUE)
    unknown_info <- rownames(sce[names[star],])
    unknown_info_sep <- paste(unknown_info, collapse = " ")
    stop("Allele information contains unknown identifier.
         Please check the data and remove rows of the following
         allele features identifiers: `",unknown_info_sep, "` or use proper nomenclature.")
  }
}

#' Find not yet known allele identifiers
#'
#' @description
#' Internal function used in `get_allelecounts()`to find allele identifier that are not present in the lookup table.
#'
#' @param sce A \code{\link{SingleCellExperiment}} object.
#' @param agene_names A list of allele gene names.
#'
#' @importFrom SingleCellExperiment counts
#'
#' @return A list of character strings of identifiers that can not be found in the allele lookup table.
find_not_ident <- function(sce, agene_names){
  # return allele genes that do not start with HLA (not found in lookup table)
  scae_counts <- counts(get_alleles(sce))
  rownames(scae_counts) <- agene_names
  # This needs to be fixed on the long run, because not all genes with extended information
  # match HLA
  not_ids <- rownames(scae_counts[!grepl(c("^HLA"), rownames(scae_counts)), , drop = FALSE])

  not_ids
}

#' Build new substring
#'
#' @description
#' Internal function used in `get_allelecounts()`. Function is used to cut character string at "*" character and return it. Only used to cut the names
#' if unknown alleles, as the lookup table does not provide any corresponding immune gene name.
#'
#' @param allele_id A list of character strings for unidentified allele_names that have proper [allele-nomenclature](https://hla.alleles.org/nomenclature/index.html).
#'
#' @return A list of charcter strings resembling HLA gene names.
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
#' @param sce A \code{\link{SingleCellExperiment}} object.
#' @param lookup A data.frame object containing the lookup table.
#' @param exp_type A vector containing two character strings. Either `"WTA"` or `"Amplicon"` are valid inputs. Choose one depending on the used transcriptomics approach.
#'
#' @importFrom SingleCellExperiment counts
#'
#' @return A SingleCellExperiment object.
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
      message(allele_ids_lookup[i], " can't be found in the lookup table")
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
    return_known   <- c(alleletogene_counts)
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
#' @param sce A \code{\link{SingleCellExperiment}} object.
#' @param lookup A data.frame object containing the lookup table.
#' @param exp_type A vector containing two character strings. Either `"WTA"` or `"Amplicon"` are valid inputs. Choose one depending on the used transcriptomics approach.
#'
#' @importFrom Matrix colSums
#' @importFrom SummarizedExperiment rowData<- colData<-
#' @importFrom SingleCellExperiment rowData colData SingleCellExperiment
#' @importFrom BiocGenerics rbind
#'
#' @return A SingleCellExperiment object.
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

  if (unknown) {
    al_gene <- al_gene[!(rownames(al_gene) %in% not_ids), , drop = FALSE]
    filtered_rows <- vapply(rownames(rowData(sce)), function(rowname) {
      any(startsWith(rowname, not_ids))
    }, logical(1))
    rowData(sce)[filtered_rows, "Quant_type"] <- "A_unknown"
  }

  al_sce <- SingleCellExperiment(assays = list(counts = al_gene),
                                 colData = colData(sce))
  rowData(al_sce)$Symbol <- rownames(al_gene)
  if (exp_type == "WTA"){
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
#' @param exp_type A vector containing two character strings. Either `"WTA"` or `"Amplicon"` are valid inputs. Choose one depending on the used transcriptomics approach.
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
  if (exp_type == "WTA"){
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


#-6-------------------------add sample tags------------------------------------#

#' Adding sample tag information to colData
#'
#' @description
#' Internal function used in 'readAlleleCounts()'. Stated here because its supposed to be a transformation step of the SCAE object.
#' Adding sample tag information to colData
#'
#' @param path character string input containing the path to the directory containing the
#'   input files
#' @param scae A \code{\link{SingleCellAlleleExperiment}} object.
#' @param tag_feature_mtx A character string determining the name of the file containing the sample-tag quantification data.
#' @param tag_feature_barcodes A character string determining the name of the file containing the sample-tag barcode identifiers.
#'
#' @importFrom methods as
#' @importFrom utils read.table
#' @importFrom Matrix readMM rowSums
#' @importFrom MatrixGenerics rowMaxs
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SingleCellExperiment colData
#'
#' @return A SingleCellAlleleExperiment object.
add_sample_tags <- function(path, scae, tag_feature_mtx, tag_feature_barcodes){

  tags  <- Matrix::readMM(paste0(path, tag_feature_mtx, sep = ""))
  cells <- utils::read.table(paste0(path, tag_feature_barcodes, sep = ""), header = FALSE)
  rownames(tags) <- cells$V1
  colnames(tags) <- paste("ST_", seq_len(12), sep = "")
  tags <- methods::as(tags, "CsparseMatrix")
  dom_tags <- data.frame(Cellular_Barcode = rownames(tags),
                         Sample_Tag = ifelse(rowMaxs(tags) >= 0.75 * rowSums(tags),
                                             colnames(tags)[max.col(tags)],
                                             "Multiplet"))
  dom_tags <- dom_tags[!dom_tags$Sample_Tag == 'Multiplet', ]

  inter <- intersect(rownames(colData(scae)), dom_tags$Cellular_Barcode)
  colData(scae)$sample_tags <- NA
  colData(scae[, inter])$sample_tags <- dom_tags[inter, "Sample_Tag"]
  scae <- scae[, !is.na(colData(scae)$sample_tags)]
  scae
}

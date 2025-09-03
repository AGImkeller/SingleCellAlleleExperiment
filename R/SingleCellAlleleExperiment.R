
#' The SingleCellAlleleExperiment class
#'
#' The SingleCellAlleleExperiment class is a comprehensive multi-layer data
#' structure, enabling the representatino of immune genes at specific levels,
#' including alleles, genes and groups of functionally similar genes. This data
#' representation allows data handling and data analysis across these
#' immunological relevant, different layers of annotation.
#'
#' The SingleCellAlleleExperiment class builds upon and extends the data
#' representation that can be facilitated using a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object.
#'
#' The Constructor `SingleCellAlleleExperiment()` can be used on its own,
#' if raw data is processed accordingly (see examples) OR in a more
#' convenient way using this packages read in function `read_allele_counts()`
#'
#' A getter function `scae_subset()` allows to subset the object according to
#' the newly implemented layers.
#'
#'
#' @seealso [read_allele_counts()]
#' @seealso [scae_subset()]
#'
#' @param ... Arguments passed to the \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' constructor to fill the slots of the SCE-class.
#' @param lookup A data.frame object containing the lookup table.
#' @param metadata A list containing a dataframe and two integer values of
#' information regarding plotting a knee plot for quality control. This parameter
#' is linked to `filter_mode="yes"` in the `read_allele_counts()` function.
#' @param threshold An integer value used as a threshold for filtering
#' low-quality barcodes/cells.
#' @param exp_type Internal character string parameter that determines in which
#' format the gene symbols in the input data are. Can be `c("ENS","noENS")`
#' @param log A logical parameter which determines if the user wants to
#' compute the `logcounts` assay.
#' @param gene_symbols A logical parameter to decide whether to compute additional
#' gene gene symbols in case the raw data only contains ENSEMBL gene identifiers.
#' @param verbose A logical parameter to decide if runtime-messages should be
#' shown during function execution. Use `FALSE` if no info runtime-messages
#' should be shown (default), and `TRUE` for showing runtime-messages.
#'
#' @details
#' In this class, similar to the \code{\link[SingleCellExperiment]{SingleCellExperiment}} class,
#' rows should represent genomic features (including immune genes, represented
#' as allele information), while columns represent single cells/barcodes.
#'
#' The SingleCellAlleleExperiment data structure serves as a data representation
#' for data generated with the `scIGD` workflow.
#' This workflow allows for the quantification of expression and interactive
#' exploration of donor-specific alleles of different immune genes and its
#'
## @seealso
## \code{\link{https://github.com/AGImkeller/scIGD/}},
## for information about the "single-cell ImmunoGenomic Diversity" **scIGD** workflow.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment sizeFactors counts
#' logcounts counts<- logcounts<-
#' @importFrom SummarizedExperiment assays<- assays metadata<-
#' @importFrom DelayedArray DelayedArray
#' @importFrom methods new
#'
#' @examples
#' ##-If you want to use the Constructor on its own, some preprocessing is
#' ##-necessary to bring the data in proper format
#' ##-Here, we use an example dataset found in in the `scaeData` package.
#'
#' ##-Find an alternative and recommended read in below as a second example
#'
#' example_data_5k <- scaeData::scaeDataGet(dataset="pbmc_5k")
#' lookup_name <- "pbmc_5k_lookup_table.csv"
#' lookup <- read.csv(system.file("extdata", lookup_name, package="scaeData"))
#'
#' barcode_loc <- file.path(example_data_5k$dir, example_data_5k$barcodes)
#' feature_loc <- file.path(example_data_5k$dir, example_data_5k$features)
#' matrix_loc  <- file.path(example_data_5k$dir, example_data_5k$matrix)
#'
#' feature_info <- utils::read.delim(feature_loc, header=FALSE)
#' cell_names   <- utils::read.csv(barcode_loc, sep="", header=FALSE)
#' mat <- t(Matrix::readMM(matrix_loc))
#'
#' ##-Prepare input data
#' colnames(feature_info) <- "Ensembl_ID"
#' sample_names <- "pbmc_5k"
#' sparse_mat <- as(mat, "CsparseMatrix")
#'
#' ##--colData
#' cell_info_list <- S4Vectors::DataFrame(Sample=rep(sample_names,
#'                                                  length(cell_names)),
#'                                       Barcode=cell_names$V1,
#'                                       row.names=NULL)
#' ##--rowData and count matrix
#' rownames(feature_info) <- feature_info[,1]
#' cnames <- cell_info_list$Barcode
#' colnames(sparse_mat) <- cnames
#'
#' scae <- SingleCellAlleleExperiment(assays=list(counts=sparse_mat),
#'                                    rowData=feature_info,
#'                                    colData=cell_info_list,
#'                                    lookup=lookup,
#'                                    verbose=TRUE)
#'
#' scae
#'
#' ##-OR, use the read in function `read_allele_counts()` !![RECOMMENDED]!!
#' ##-Find more examples in its documentation using `?read_allele_counts`
#'
#' # scae_2 <- read_allele_counts(example_data_5k$dir,
#' #                              sample_names="example_data",
#' #                              filter_mode="no",
#' #                              lookup_file=lookup,
#' #                              barcode_file=example_data_5k$barcodes,
#' #                              gene_file=example_data_5k$features,
#' #                              matrix_file=example_data_5k$matrix,
#' #                              verbose=TRUE)
#'
#' # scae_2
#'
#' @return A SingleCellAlleleExperiment object.
#' @export
SingleCellAlleleExperiment <- function(...,
                                       lookup,
                                       metadata=NULL,
                                       threshold=0,
                                       exp_type="ENS",
                                       log=TRUE,
                                       gene_symbols=FALSE,
                                       verbose=FALSE){
  sce <- SingleCellExperiment(...)

  ## input checks for optional package installation of optional functionalities
  check_valid_optional_package(log=log, gene_symbols=gene_symbols)


  sce_add_look <- ext_rd(sce, exp_type, gene_symbols, verbose=verbose)

  if (verbose){
    message("  Generating SCAE object: Extending rowData with new classifiers")
  }

  sce_filtered <- sce_add_look[, colSums(counts(sce_add_look)) > threshold]

  if (verbose){
    message("  Generating SCAE object: Filtering at ", threshold, " UMI counts.")
  }

  if(log){
    sce_filtered <- scuttle::computeLibraryFactors(sce_filtered)
    if (verbose){
      message("  Generating SCAE object: Compute Library Factors before adding new layers")
    }
  }

  scae <- alleles2genes(sce_filtered, lookup, exp_type, gene_symbols)

  if (verbose){
    message("  Generating SCAE object: Aggregating alleles corresponding to the same gene")
  }

  scae <- genes2functional(scae, lookup, exp_type, gene_symbols)

  if (verbose){
    message("  Generating SCAE object: Aggregating genes corresponding to the same functional groups")
  }

  if(log){
    normed_counts <- scuttle::normalizeCounts(scae, size_factors=sizeFactors(scae), transform="log")
    assays(scae)$logcounts  <- normed_counts
    logcounts(scae) <- DelayedArray::DelayedArray(logcounts(scae))
    if (verbose){
      message("  Generating SCAE object: Generate logcounts assay using Library Factors")
    }
  }

  counts(scae) <- DelayedArray::DelayedArray(counts(scae))
  metadata(scae)[["knee_info"]] <- metadata

  scae <- new("SingleCellAlleleExperiment", scae)

  return(scae)
}


##-----------------Internal functions used in the SCAE-Constructo-------------##
##-1--------------------------------ext_rd------------------------------------##

#' Extend rowData with new annotation columns
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param exp_type Internal character string parameter that determines in which
#' format the gene symbols in the input data are. Can be `c("ENS","noENS")`
#' @param gene_symbols A logical parameter to decide whether to compute additional
#' gene gene symbols in case the raw data only contains ENSEMBL gene identifiers.
#' @param verbose A logical parameter to decide if runtime-messages should be
#' shown during function execution. Use `FALSE` if no info runtime-messages
#' should be shown (default), and `TRUE` for showing runtime-messages.
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom SingleCellExperiment rowData
#' @return A SingleCellExperiment object
ext_rd <- function(sce, exp_type, gene_symbols, verbose=FALSE){

  allele_names_all <- find_allele_ids(sce)

  ## Group of genes for which extended informaton is stored
  rowData(sce[allele_names_all,])$NI_I <- "I"
  ## Allele level
  rowData(sce[allele_names_all,])$Quant_type <- "A"
  ## Group of genes for which classical (gene level) informaton is stored
  rn_in_alleles <- rownames(sce) %in% allele_names_all
  rowData(sce)[!(rn_in_alleles), ]$NI_I <- "NI"
  ## Gene level
  rowData(sce)[!(rn_in_alleles), ]$Quant_type <- "G"

  if (exp_type == "ENS" && gene_symbols){
    if (verbose){
      message("Using org.Hs to retrieve NCBI gene identifiers.")
    }

        gene_symbols <- get_ncbi_org(sce)
        rowData(sce)$Symbol <- gene_symbols
        rn_sce <- rownames(rowData(scae_subset_alleles(sce)))
        rowData(sce)[rn_sce,]$Symbol <- rn_sce
  }
  sce
}

## Code provided by Ahmad Al Ajami
#' Get NCBI genes using the org.HS.db package
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @importFrom methods as
#' @importFrom SingleCellExperiment rowData
#' @return A list of character strings for gene names.
get_ncbi_org <- function(sce){

      ensembl_ids <- rowData(sce)$Ensembl_ID
      ensembl_ids <- sub("\\..*", "", ensembl_ids)

      Hs_symbol <- org.Hs.eg.db::org.Hs.egSYMBOL
      Hs_ensembl <- org.Hs.eg.db::org.Hs.egENSEMBL
      mapped_Hs_genes_symbol <- AnnotationDbi::mappedkeys(Hs_symbol)
      mapped_Hs_genes_ensembl <- AnnotationDbi::mappedkeys(Hs_ensembl)
      Hs_symbol_df <- as.data.frame(Hs_symbol[mapped_Hs_genes_symbol])
      Hs_ensembl_df <- as.data.frame(Hs_ensembl[mapped_Hs_genes_ensembl])

      Hs_mapping <- merge(Hs_symbol_df, Hs_ensembl_df)

      indic <- match(ensembl_ids, Hs_mapping$ensembl_id)
      ncbi_symbols <- Hs_mapping$symbol[match(ensembl_ids, Hs_mapping$ensembl_id)]

      return(ncbi_symbols)
}


##-2-----------------------------allele2genes---------------------------------##

#' Identify rows containing allele information
#'
#' @description
#' Internal function used in `get_allelecounts()` to subsample the
#' quantification assay and only return the rows specifying
#' allele-quantification information.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @importFrom SingleCellExperiment counts
#' @return A SingleCellExperiment object
find_allele_ids <- function(sce){
  rn_counts_sce <- rownames(counts(sce))
  a <- grepl("*", rn_counts_sce, fixed=TRUE)
  if (sum(a) == 0){
    a <- grepl("HLA-", rn_counts_sce, fixed=TRUE)
  }

  allele_names_all <- rownames(counts(sce)[a,])
  allele_names_all
}

#' Get Subassay with allele gene names and raw allele quantification
#'
#' @description
#' Internal function used to build a subassay containing counts from raw alleles.
#' The rownames  of this subassay are already translated to the corresponding
#' immune gene identifier, which are extracted from the lookup table.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param lookup A data.frame object containing the lookup table.
#' @importFrom SingleCellExperiment counts
#' @return A SingleCellExperiment object
get_allelecounts <- function(sce, lookup){

  allele_ids_lookup <- find_allele_ids(sce)
  list_alid <- list()

  for (i in seq_along(allele_ids_lookup)){
    new_ids <- list(lookup[grepl(allele_ids_lookup[i], lookup$Allele, fixed=TRUE),]$Gene)
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
#' Internal function for the first assay extension used in the
#' `SingleCellAlleleExperiment()` constructor, computing the first of the two
#' new subassays that get appended to the quantification assay.
#' This subassay contains the allele gene identifiers instead of the
#' allele identifiers present in the raw data and sums up the expression counts
#' of alleles that have the same allele gene identifiers.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param lookup A data.frame object containing the lookup table.
#' @param exp_type Internal character string parameter that determines in which
#' format the gene symbols in the input data are. Can be `c("ENS","noENS")`
#' @param gene_symbols A logical parameter to decide whether to compute additional
#' gene gene symbols in case the raw data only contains ENSEMBL gene identifiers.
#' @importFrom SingleCellExperiment rowData colData SingleCellExperiment rbind
#' @importFrom SummarizedExperiment rowData<- colData<-
#' @importFrom Matrix colSums
#' @return A SingleCellExperiment object
alleles2genes <- function(sce, lookup, exp_type, gene_symbols){

  v_acounts <- get_allelecounts(sce, lookup)

  alleletogene_counts <- v_acounts[1][[1]]

  uniqs <- unique(rownames(alleletogene_counts))
  al_gene <- matrix(0, nrow=length(uniqs), ncol=ncol(alleletogene_counts))
  rownames(al_gene) <- uniqs

  for (i in seq_along(uniqs)){
    rn_a2g <- rownames(alleletogene_counts)
    uniq_sum <- colSums(alleletogene_counts[rn_a2g %in% uniqs[i], , drop=FALSE])
    al_gene[i,] <- uniq_sum
  }

  al_sce <- SingleCellExperiment(assays=list(counts=al_gene),
                                 colData=colData(sce))

  if (exp_type == "ENS"){
    rowData(al_sce)$Ensembl_ID <- rownames(al_gene)
  }
  if (gene_symbols){
    rowData(al_sce)$Symbol <- rownames(al_gene)
  }

  new_sce <- SingleCellExperiment::rbind(sce, al_sce)

  uniq_genes_sce <- rownames(new_sce) %in% uniqs
  rowData(new_sce[uniq_genes_sce])$NI_I <- "I"
  rowData(new_sce[uniq_genes_sce])$Quant_type <- "G"
  return(new_sce)
}


##-3------------------------------genes2func----------------------------------##

#' Building second new subassay for the SingleCellAlleleExperiment object
#'
#' @description
#' Internal function for the second assay extension used in the
#' `SingleCellAlleleExperiment()` constructor, computing the second of the two
#' new subassays that get appended to the quantification assay. This subassay
#' contains the functional allele classes and sums up the expression counts of
#' the allele genes that are in the same functional group.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param lookup A data.frame object containing the lookup table.
#' @param exp_type Internal character string parameter that determines in which
#' format the gene symbols in the input data are. Can be `c("ENS","noENS")`
#' @param gene_symbols A logical parameter to decide whether to compute additional
#' gene gene symbols in case the raw data only contains ENSEMBL gene identifiers.
#' @importFrom SingleCellExperiment colData counts SingleCellExperiment rbind
#' @importFrom SummarizedExperiment colData<- rowData<-
#' @importFrom Matrix colSums
#' @return A SingleCellExperiment object
genes2functional <- function(sce, lookup, exp_type, gene_symbols){

  ## find functional classes for each gene
  gene_names <- rownames(get_agenes(sce))

  list_func  <- list()
  for (i in seq_along(gene_names)){
    func_classes <- lookup$Function[lookup$Gene %in% gene_names[i]][1]
    list_func[[length(list_func) + 1]] <- func_classes
  }

  gene_func_names <- unlist(list_func)
  gene2func_counts <- counts(get_agenes(sce))

  rn_g2f <- rownames(gene2func_counts)
  rn_g2f <- gene_func_names
  uniqs <- unique(rn_g2f)
  gene_func <- matrix(0, nrow=length(uniqs), ncol=ncol(sce[1,]))
  rownames(gene_func) <- uniqs

  for (i in seq_along(uniqs)){
    gene_colsums <- colSums(gene2func_counts[rn_g2f %in% uniqs[i], , drop=FALSE])
    gene_func[i,] <- gene_colsums
  }

  func_sce <- SingleCellExperiment(assays=list(counts=gene_func),
                                   colData=colData(sce))
  if (exp_type == "ENS"){
    rowData(func_sce)$Ensembl_ID <- rownames(func_sce)
  }
  if (gene_symbols){
    rowData(func_sce)$Symbol <- rownames(func_sce)
  }

  final_scae <- SingleCellExperiment::rbind(sce, func_sce)
  uniq_func_sce <- rownames(final_scae) %in% uniqs

  ## Genes with extended quantification
  rowData(final_scae[uniq_func_sce])$NI_I <- "I"
  ## Functional level
  rowData(final_scae[uniq_func_sce])$Quant_type <- "F"
  final_scae
}

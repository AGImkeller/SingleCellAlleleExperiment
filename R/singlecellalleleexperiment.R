library(SingleCellExperiment)

#' SingleCellAlleleExperiment-Class
#'
#' @return definition for the scae class
.scae <- setClass("SingleCellAlleleExperiment", contains = "SingleCellExperiment")

SingleCellAlleleExperiment <- function(..., lookup){
  sce <- SingleCellExperiment(...)

  sce_add_look <- lookup_in(sce)
  scae <- alleles2genes(sce_add_look, lookup)
  scae <- genes2functional(scae, lookup)
  .scae(scae)
}

#-------------------------allele2genes + helpers--------------------------------
#####
find_allele_ids <- function(sce){
  allele_char = "*"
  all_genes = rownames(counts(sce))
  a <- grepl(allele_char, all_genes, fixed = TRUE)
  allele_names_all <- rownames(counts(sce)[a,])
  allele_names_all
}

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
  alleletogene_counts
}

#adds alleleids from the lookup to the ID-column of the allele_rows
add_allelids <- function(sce, lookup){
  wor_copy <- sce
  allele_ids_lookup <- find_allele_ids(wor_copy)

  list_alid_lookup <- list()
  for (i in 1:length(allele_ids_lookup)){
    new_aids <- list(lookup[grepl(allele_ids_lookup[i], lookup$Allele, fixed = TRUE),]$AlleleID)
    list_alid_lookup[[length(list_alid_lookup) + 1]] <- new_aids
  }
  list_alid_lookup
  allelids <- unlist(list_alid_lookup)

  rowData(wor_copy[allele_ids_lookup,])$ID <- allelids
  #print(rowData(wor_copy[allele_ids_lookup,]))
  wor_copy
}

#finds alleles, looks the corresponding gene_names in lookup_table and performs
#colSums for the alleles with the same gene_name
alleles2genes <- function(sce, lookup){
  w_copy <- sce
  #allele counts with allele_gene names as rownames
  alleletogene_counts <- get_allelecounts(w_copy, lookup)
  uniqs <- unique(rownames(alleletogene_counts))
  #build new matrix for the colSum
  al_gene <- matrix(0, nrow = length(uniqs), ncol = length(alleletogene_counts[1,]), )
  al_gene <- as(al_gene, "sparseMatrix")
  rownames(al_gene) <- uniqs

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
  #new_sce <- add_allelids(new_sce, lookup)
  #test to compare the counts
  test_counts <- sum(alleletogene_counts[rownames(alleletogene_counts) %in% uniqs[1],])
  test_colSum <- sum(al_gene[rownames(al_gene) %in% uniqs[1],])

  cat(paste("Counts of", uniqs[1] ,"before colSum: ",test_counts,"\n"))
  cat(paste("Counts of", uniqs[1], "after colSum: ",test_colSum, "\n"))
  cat(paste("---------------------------------------", "\n"))

  if (test_counts == test_colSum){
    return(new_sce)
  }else{
    print("Counts and ColSum do not match", "\n")
  }
}

#####

#--------------------------genes2func + helpers---------------------------------
#####
get_agene_names <- function(sce){
  func_work_copy <- sce
  allele_star <- rownames(counts(func_work_copy)[grepl("*", rownames(sce), fixed = TRUE),])
  last_allele <- rownames(counts(func_work_copy[allele_star,]))[length(allele_star)]
  index <- which(rownames(func_work_copy) == last_allele)
  hla_names <- rownames(counts(func_work_copy)[(index+1):nrow(counts(func_work_copy)),])
  hla_names
}

genes2functional <- function(sce, lookup){
  func_copy <- sce
  #alle hla_names zu gene_names benennen
  gene_names <- get_agene_names(func_copy)
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

  #für jede klasse kombinieren wir die gencounts von allen genen, dieser klasse
  for (i in 1:length(func_group_genes)){
    gene_colsums <- colSums(counts(func_copy[unlist(func_group_genes[i]),]))
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

  #testing the counts before and after colSums
  #just additional information during the progress
  func_one <- sum(counts(func_copy[func_group_genes[[1]],]))
  func_two <- sum(counts(func_copy[func_group_genes[[2]],]))

  cat(paste("Counts of", func_group_genes[1] ,"before colSum: ", func_one,"\n"))
  cat(paste("Counts of", func_group_genes[2] ,"before colSum: ", func_two,"\n"))

  cat(paste("Counts of", rownames(gene_func)[1], "after colSum: ", sum(gene_func[1,]), "\n"))
  cat(paste("Counts of", rownames(gene_func)[2], "after colSum: ", sum(gene_func[2,])))

  final_scae
}
#####
#-----------------------------common utils--------------------------------------
#####
lookup_in <- function(sce){
  new_sce <- sce
  allele_names_all <- find_allele_ids(new_sce)

  rowData(new_sce[allele_names_all,])$NI_I <- "I"
  rowData(new_sce[allele_names_all,])$Quant_type <- "A"

  rowData(new_sce)[!(rownames(new_sce) %in% allele_names_all), ]$NI_I <- "NI"
  rowData(new_sce)[!(rownames(new_sce) %in% allele_names_all), ]$Quant_type <- "G"
  new_sce
}

lookup_gene <- function(sce, lookup){
  new_sce_comb <- sce
  names_check <- get_allelecounts(new_sce_comb, lookup)
  uniqs_comb <- unique(rownames(names_check))

  rowData(new_sce_comb[rownames(new_sce_comb) %in% uniqs_comb])$NI_I <- "I"
  rowData(new_sce_comb[rownames(new_sce_comb) %in% uniqs_comb,])$Quant_type <- "G"
  new_sce_comb
}

get_alleles <- function(scae){
  alleles <- scae[rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "A"]
  alleles
}

get_agenes <- function(scae){
  agenes <- scae[rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "G"]
  agenes
}

get_func <- function(scae){
  func <- scae[rowData(scae)$NI_I == "I" & rowData(scae)$Quant_type == "F"]
  func
}

#####


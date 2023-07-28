################################################################################
############---function plots a knee plot ranking barcode counts---#############
################---majority of code taken from Ahmad Al Ajami---################
################################################################################


#------------------------------Knee plot---------------------------------------#

#####

#' Knee plot
#'
#' @description
#' Creates a Kneeplot ranking the barcodes according to their total UMI count. The plot is
#' used to determine a threshold for filtering barcodes in the preprocessing step.
#'
#'
#' @param path file path of the directory containing the input files as character string
#'
#' @import dplyr
#' @import tibble
#' @import ggplot2
#' @importFrom Matrix rowSums readMM
#' @importFrom utils read.csv read.delim
#'
#' @return returns a knee plot for determining a count threshold used for filtering out barcodes
#' @export
plotKnee <- function(path){

  barcodes_loc <- file.path(path, "cells_x_genes.barcodes.txt")
  features_loc <- file.path(path, "cells_x_genes.genes.txt")
  matrix_loc   <- file.path(path, "cells_x_genes.mtx")

  barcodes <- read.csv(barcodes_loc, sep = "", header = FALSE)
  features <- read.delim(features_loc, header = FALSE)
  matrix   <- readMM(matrix_loc)

  rownames(matrix) <- paste(barcodes$V1, "_D1", sep = "")
  colnames(matrix) <- features$V1
  colnames(matrix) <- gsub("\\|.*","", colnames(matrix))
  colnames(matrix) <- gsub("\\-",".", colnames(matrix))

  tot_counts <- Matrix::rowSums(matrix)

  df <- tibble(cells = rownames(matrix),
               total = tot_counts,
               rank = row_number(dplyr::desc(total))) %>%
    distinct() %>%
    dplyr::arrange(rank)

  total <- df$total

  ggplot(df, aes(total, rank)) +
    geom_path() +
    scale_x_log10() + scale_y_log10() + annotation_logticks() +
    labs(y = "Barcode rank", x = "Total UMI count")
}
#####

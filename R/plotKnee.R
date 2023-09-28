################################################################################
############---function plots a knee plot ranking barcode counts---#############
#############--Code for default knee plot is from Ahmad Al Ajami--##############
################################################################################

#------------------------------Knee plot---------------------------------------#

#####

#' Knee plot
#'
#' @description
#' Creates a knee plot ranking the barcodes according to their total UMI count. The plot is
#' used to determine a threshold for filtering barcodes in the preprocessing step. You can plot two different knee plots.
#' The `mode = "default` knee plot shows a default knee plot. The `mode = "advanced"` knee plot computes a knee- and inflection point
#' which can be used as a suggestion for filtering out low-quality cells.
#'
#'
#' @param path character string of the file path to the directory containing the expected input files
#' @param mode choose knee plot version. `mode = "default"` for the default knee plot. `mode = "advanced"` for the advanced knee plot including a knee- and inflection point. The standard value is the default knee plot.
#'
#' @importFrom dplyr desc arrange row_number distinct
#' @import tibble
#' @import ggplot2
#' @importFrom Matrix rowSums readMM
#' @importFrom utils read.csv read.delim
#' @importFrom S4Vectors metadata
#'
#' @return returns a knee plot for determining a count threshold used for filtering out barcodes
#' @export
plotKnee <- function(path, mode = NULL){

  # check validity of mode parameter
  if (is.null(mode)) {
    mode <- "default"
  } else {
    if (!mode %in% c("default", "advanced")) {
      stop("Invalid value for mode parameter. Allowed values are 'default' and 'advanced'.")
    }
  }

  # expected file names
  barcodes_loc <- file.path(path, "cells_x_genes.barcodes.txt")
  features_loc <- file.path(path, "cells_x_genes.genes.txt")
  matrix_loc   <- file.path(path, "cells_x_genes.mtx")

  barcodes <- read.csv(barcodes_loc, sep = "", header = FALSE)
  features <- read.delim(features_loc, header = FALSE)
  matrix   <- readMM(matrix_loc)

  # default knee plot
  if (mode == "default"){

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

    default <- ggplot(df, aes(x = rank, y = total)) +
                geom_path() +
                scale_x_log10() + scale_y_log10() + annotation_logticks() +
                labs(y = "Total UMI counts", x = "Barcode rank")

   return(default)
  }

  # advanced knee plot with knee- and inflection point
  if (mode == "advanced"){
    matrix <- t(matrix)

    br.out <- DropletUtils::barcodeRanks(matrix)
    names(br.out)
    fitteddf <- br.out$fitted

    data_df <- data.frame(rank = br.out$rank, total = br.out$total, fitteddf = br.out$fitted)

    gg <- ggplot(data_df, aes(x = rank, y = total)) +
          geom_point() +
          geom_line(aes(y = fitteddf), color = "red") +
          scale_x_log10() +
          scale_y_log10() +
          annotation_logticks() +
          labs(x = "Barcode rank", y = "Total UMI count") +
          geom_hline(yintercept = S4Vectors::metadata(br.out)$knee, color = "dodgerblue", linetype = "dashed") +
          geom_hline(yintercept = S4Vectors::metadata(br.out)$inflection, color = "forestgreen", linetype = "dashed") +
          annotate("text", x = 2, y = S4Vectors::metadata(br.out)$knee * 1.2, label = "knee", color = "dodgerblue") +
          annotate("text", x = 2.25, y = S4Vectors::metadata(br.out)$inflection * 1.2, label = "inflection", color = "forestgreen") +
          theme_bw()
    gg
  }
}
#####

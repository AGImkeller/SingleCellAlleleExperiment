################################################################################
############---function plots a knee plot ranking barcode counts---#############
################---majority of code taken from Ahmad Al Ajami---################
################################################################################


#------------------------------Knee plot---------------------------------------#

#####

#' Knee plot
#'
#' @description
#' Creates a knee plot ranking the barcodes according to their total UMI count. The plot is
#' used to determine a threshold for filtering barcodes in the preprocessing step.
#'
#'
#' @param matrix_file count matrix that was read in in the `read_from_sparse_allele()` function
#' @param genes_file list with gene identifiers that was read in in the `read_from_sparse_allele()` function
#' @param barcodes_file list with barcode identifiers that was read in in the `read_from_sparse_allele()` function
#'
#' @import dplyr
#' @import tibble
#' @import ggplot2
#' @importFrom Matrix rowSums readMM
#' @importFrom utils read.csv read.delim
#' @importFrom S4Vectors metadata
#'
#' @return returns a knee plot for determining a count threshold used for filtering out barcodes
#' @export
plotKnee <- function(matrix_file, genes_file, barcodes_file){

  barcodes <- barcodes_file
  features <- genes_file
  matrix <- matrix_file

  #advanced knee plot
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
  print(gg)

  inflection_p <- S4Vectors::metadata(br.out)$inflection
  return(inflection_p)

}
#####

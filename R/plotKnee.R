
#------------------------------Knee plot---------------------------------------#

#' Knee plot
#'
#' @description
#' Creates a knee plot ranking the barcodes according to their total UMI count. The plot is
#' used to determine a threshold for filtering barcodes in the preprocessing step.
#'
#'
#' @param matrix A sparse \code{\link{Matrix}} object containing the quantification data.
#' @param genes A data.frame object containing gene identifiers.
#' @param barcodes A data.frame object containing barcode identifiers.
#'
#' @import ggplot2
#' @importFrom Matrix rowSums readMM
#' @importFrom utils read.csv read.delim
#' @importFrom S4Vectors metadata
#'
#' @return A knee plot about the quantification data.
plot_knee <- function(matrix, genes, barcodes){

  barcodes <- barcodes
  features <- genes
  matrix <- matrix

  #advanced knee plot
  br_out <- DropletUtils::barcodeRanks(matrix)

  names(br_out)
  fitteddf <- br_out$fitted
  total <-  br_out$total

  data_df <- data.frame(rank=br_out$rank, total=br_out$total, fitteddf=br_out$fitted)

  gg <- ggplot(data_df, aes(x=rank, y=total)) +
        geom_point() +
        geom_line(aes(y=fitteddf), color="red") +
        scale_x_log10() +
        scale_y_log10() +
        annotation_logticks() +
        labs(x="Barcode rank", y="Total UMI count") +
        geom_hline(yintercept=S4Vectors::metadata(br_out)$knee, color="dodgerblue", linetype="dashed") +
        geom_hline(yintercept=S4Vectors::metadata(br_out)$inflection, color="forestgreen", linetype="dashed") +
        annotate("text", x=2, y=S4Vectors::metadata(br_out)$knee * 1.2, label="knee", color="dodgerblue") +
        annotate("text", x=2.25, y=S4Vectors::metadata(br_out)$inflection * 1.2, label="inflection", color="forestgreen") +
        theme_bw()
  suppressWarnings(print(gg))

  inflection_p <- S4Vectors::metadata(br_out)$inflection
  return(inflection_p)
}

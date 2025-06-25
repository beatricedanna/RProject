# OKAY
# To execute the function:
# peaks_vector <- TotalAccessibilityPerPeak(results)

#' TotalAccessibilityPerPeak
#'
#' @param atac_seq_dt A data.table with peak regions as rows and cells as columns. The first column should be 'feature_id'.
#' @return A named numeric vector: total accessibility per peak region.
#'
#' @export
TotalAccessibilityPerPeak <- function(atac_seq_dt){
  peaks_sums_vector <- rowSums(atac_seq_dt[, -1])
  names(peaks_sums_vector) <- atac_seq_dt$feature_id

  return(peaks_sums_vector)
}

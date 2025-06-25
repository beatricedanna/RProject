# OKAY
# To execute the function:
# peaks_vector <- ColSumPeaks(results)


ColSumPeaks <- function(atac_seq_dt){
  peaks_sums_vector <- rowSums(atac_seq_dt[, -1])
  names(peaks_sums_vector) <- atac_seq_dt$feature_id

  return(peaks_sums_vector)
}

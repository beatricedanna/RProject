# To execute the function:
# integration_results <- IntegrateAndSummarize(
# final_normalized_gex_dt,
# final_normalized_atac_dt,
# annotated_peaks_gr_final,
# protein_coding_gex
# )
utils::globalVariables(c(".SD", "nearest_gene_id"))

#' IntegrateAndSummarize
#'
#' @importFrom GenomicRanges seqnames start end mcols
#' @importFrom data.table .SD setnames := :=
#' @importFrom stats na.omit
#'
#' @param normalized_gex_dt A data.table containing normalized gene expression data.It must include a "feature_id" column with gene IDs.
#' @param normalized_atac_dt A data.table containing normalized ATAC-seq peak accessibility data.It must include a "feature_id" column with peak IDs.
#' @param annotated_peaks_gr_final A GRanges object containing peaks asscoiated with a gene. It must include a "nearest_gene_id" column specifying the closest gene each peak has associated with.
#' @param finalized_expression_gr A GRanges object of protein-coding gene expression data including gene IDs in the metadata gene_id column.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{merged_data}{A data.table merging normalized gene expression and aggregated peak accessibility by gene ID.}
#'   \item{unmerged_peaks_summary}{A data.table summarizing total, merged, and unmerged peaks counts.}
#'   \item{unassociated_genes_summary}{A data.table summarizing total genes, genes with associated peaks, and genes without associated peaks.}
#'   \item{unassociated_gene_ids}{A vector of protein-coding gene IDs (from finalized_expression_gr) with no associated peaks.}
#' }
#' @export
IntegrateAndSummarize <- function(normalized_gex_dt, normalized_atac_dt, annotated_peaks_gr_final, finalized_expression_gr) {

  # Assign "chr:start-end" style names to each peak in annotated_peaks_gr_final
  names(annotated_peaks_gr_final) <- paste0(
    seqnames(annotated_peaks_gr_final), ":",
    start(annotated_peaks_gr_final), "-",
    end(annotated_peaks_gr_final)
  )
  # Create peak to gene map (peak_gene_map)
  peak_gene_map <- data.table(
    feature_id = names(annotated_peaks_gr_final), # To be sure that ID column is called "feature_id" as in normalized_atac_dt
    nearest_gene_id = mcols(annotated_peaks_gr_final)$nearest_gene_id
  ) # Create a data.table with 2 columns: "feature_id" (from annotated_peaks_gr_final) and "nearest_gene_id" (from annotated_peaks_gr_final)

  # Merge of normalized_atac_dt and peak_gene_map by "feature_id" (contained in both dt)
  atac_annotated_dt <- merge(
    normalized_atac_dt,
    peak_gene_map,
    by = "feature_id",
    all.x = TRUE # Mantain all peaks (insert NA if no association)
  )

  atac_annotated_dt <- na.omit(atac_annotated_dt, cols = "nearest_gene_id") # Delete peaks with NA
  cell_cols <- grep("^[A-Z]", names(atac_annotated_dt), value = TRUE) # Save in cell_cols all columns starting with a capital letter
  aggregated_atac_dt <- atac_annotated_dt[, lapply(.SD, sum), by = nearest_gene_id, .SDcols = cell_cols] # Add values for peaks associated with the same nearest_gene_id (by = nearest_gene_id groups them and lapply(.SD, sum) add them), (.SDcols = cell_cols specifies that the sum should be done only on cell columns)

  # Rename the column in aggregated_atac_dt and normalized_gex_dt to merge them
  setnames(aggregated_atac_dt, "nearest_gene_id", "gene_id")
  setnames(normalized_gex_dt, "feature_id", "gene_id")

  # Merge aggregated_atac_dt and normalized_gex_dt to obtain a new data.table containing for each "gene_id" the nomalized gene expression columns and the accessibility columns for that gene
  merged_data <- merge(normalized_gex_dt, aggregated_atac_dt, by = "gene_id")

  merged_gene_ids <- merged_data$gene_id # Extract gene_id from merged_data and save them in merged_gene_ids

  unmerged_peak_ids <- peak_gene_map[!nearest_gene_id %in% merged_gene_ids, feature_id] # Create a vector of non associated peaks
  # Create a summary data.table with number of total peaks, number of associated peaks and number of non associated peaks
  unmerged_peaks_summary <- data.table(
    Category = c("Total Peaks", "Merged Peaks", "Unmerged Peaks"),
    Count = c(length(annotated_peaks_gr_final), length(annotated_peaks_gr_final) - length(unmerged_peak_ids), length(unmerged_peak_ids))
  )

  all_expression_gene_ids <- finalized_expression_gr$gene_id # Create a vector of all gene ids from finalized_expression_gr
  unassociated_gene_ids <- setdiff(all_expression_gene_ids, merged_gene_ids) # Creates a vector of genes non asscoiated with a peak (not present in merged_gene_ids)
  # Create a summary data.table with number of total genes, number of genes with associated peaks and number of genes without associated peaks
  unassociated_genes_summary <- data.table(
    Category = c("Total Genes", "Genes with Peaks", "Genes without Peaks"),
    Count = c(length(all_expression_gene_ids), length(all_expression_gene_ids) - length(unassociated_gene_ids), length(unassociated_gene_ids))
  )

  return(list(
    merged_data = merged_data,
    unmerged_peaks_summary = unmerged_peaks_summary,
    unassociated_genes_summary = unassociated_genes_summary,
    unassociated_gene_ids = unassociated_gene_ids
  ))
}

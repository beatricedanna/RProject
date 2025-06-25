# Before running the function:

# Load the libraries:
# library(ggplot2)

# To execute the function:
# summary_plots <- GenerateSummaryPlots(
#  annotated_peaks_gr_final,
#  protein_coding_gex,
#  integration_results$unassociated_gene_ids
# )

# To print the graphs:
# peak_intensity_plot <- summary_plots$peak_intensity_plot
# print(peak_intensity_plot)
# unassociated_gex_plot <- summary_plots$unassociated_gex_plot
# print(unassociated_gex_plot)

utils::globalVariables(c("log10_total_accessibility", "log10_total_expression"))

#' GenerateSummaryPlots
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_bw theme element_text
#' @importFrom gtools mixedsort
#'
#' @param annotated_peaks_gr A GenomicRanges object containing annotated ATAC-seq peaks,with metadata column "total_accessibility" representing peak intensities.
#' @param finalized_expression_gr A GenomicRanges object containing gene expression data, with metadata column "total_expression".
#' @param unassociated_gene_ids A character vector of gene IDs that are not associated with any ATAC peak.
#'
#' @return A named list containing two ggplot2 plot objects:
#' \item{peak_intensity_plot}{Boxplot of peak intensities per chromosome.}
#' \item{unassociated_gex_plot}{Boxplot of gene expression for genes without associated peaks, per chromosome.}
#'
#' @export
GenerateSummaryPlots <- function(annotated_peaks_gr, finalized_expression_gr, unassociated_gene_ids) {

  cat("Plot of peak intensity distribution chromosome by chromosome...\n")

  peaks_plot_df <- as.data.frame(annotated_peaks_gr) # Convert annotated_peaks_gr (GRanges obj) into a data frame
  peaks_plot_df$log10_total_accessibility <- log10(peaks_plot_df$total_accessibility + 1) # Add a column named "log10_total_accessibility" that calculates log10 of total accessibility per each peak

  chrom_order <- gtools::mixedsort(unique(as.character(peaks_plot_df$seqnames))) # Order "naturally" the chromosomes (from "seqnames" column of peaks_plot_df) and convert them into characters (removing also duplicates)
  peaks_plot_df$seqnames <- factor(peaks_plot_df$seqnames, levels = chrom_order) # Transform "seqnames" column in an ordered column (based on chrom_order)

  plot_peak_intensity <- ggplot(peaks_plot_df, aes(x = seqnames, y = log10_total_accessibility, fill = seqnames)) +
    geom_boxplot(outlier.shape = NA) + # Draws the boxplot for each chromosome, hiding outliers
    labs(
      title = "Peak Intensity Distribution per chromosome",
      x = "Chromosome",
      y = "Log10(Total Accessibility + 1)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

  cat("Plot of gene expression distribution of genes without associated peaks...\n")

  valid_unassociated_ids <- intersect(unassociated_gene_ids, finalized_expression_gr$gene_id) # Find intersection between unassociated genes and genes present in finalized_expression_gr
  unassociated_genes_gr <- finalized_expression_gr[finalized_expression_gr$gene_id %in% valid_unassociated_ids, ] #Create a subset using only valid ids, preventing addition of NA

  gex_plot_df <- as.data.frame(unassociated_genes_gr)
  gex_plot_df$log10_total_expression <- log10(gex_plot_df$total_expression + 1)
  gex_plot_df$seqnames <- factor(gex_plot_df$seqnames, levels = chrom_order) # Chromosome order
  gex_plot_df <- gex_plot_df[!is.na(gex_plot_df$seqnames), ] # Delete all NA if existing

  plot_unassociated_gex <- ggplot(gex_plot_df, aes(x = seqnames, y = log10_total_expression, fill = seqnames)) +
    geom_boxplot(outlier.shape = NA) +
    labs(
      title = "Gene expression distribution of genes without associated peaks",
      x = "Chromosome",
      y = "Log10(Total Expression + 1)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") # No legend, there is X axis

  return(list(
    peak_intensity_plot = plot_peak_intensity,
    unassociated_gex_plot = plot_unassociated_gex
  ))
}

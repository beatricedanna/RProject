# Before running the function:

# Assign merged_data:
# merged_data <- integration_results$merged_data

# To execute the function:
# plot_data <- PrepareDataForScatterPlot(merged_data, protein_coding_gex)

utils::globalVariables(c(".", "gene_id", "nearest_gene_id"))

#' PrepareDataForScatterPlot
#'
#' @importFrom data.table as.data.table setnames data.table
#' @importFrom stats na.omit
#'
#' @param merged_dt A data.table containing merged gene expression and ATAC-seq data. Columns ending in `.x` are expression, `.y` are ATAC.
#' @param genes_gr A GRanges object containing gene metadata, including gene IDs and chromosomal locations.
#'
#' @return A data.table with columns:
#' \describe{
#'   \item{gene_id}{Gene identifier}
#'   \item{avg_gex_signal}{Average expression signal across all cells}
#'   \item{avg_atac_signal}{Average ATAC-seq accessibility signal across all cells}
#'   \item{chromosome}{Chromosome where the gene is located}
#' }
#'
#' @export
PrepareDataForScatterPlot <- function(merged_dt, genes_gr) {

  # Extract columns from merged_dt after merge: gene expression columns end in '.x', ATAC columns end in '.y'
  gex_cols <- names(merged_dt)[grep("\\.x$", names(merged_dt))]
  atac_cols <- names(merged_dt)[grep("\\.y$", names(merged_dt))]

  cat("Calculating average expression and ATAC signal per gene...\n")

  # Create a new summary data.table (plot_data) with 3 columns
  plot_data <- merged_dt[, .(
    gene_id = gene_id, # Column 1: gene_id
    avg_gex_signal = rowMeans(.SD[, gex_cols, with = FALSE]), # Column 2: Average (all cells) gene expression for each gene
    avg_atac_signal = rowMeans(.SD[, atac_cols, with = FALSE]) # Column 3: Average (all cells) accessibility for each gene
  )]

  cat("Adding chromosome information...\n")

  # Create a data table matching each gene to its chromosome
  gene_location_map <- as.data.table(genes_gr)[, .(gene_id, seqnames)] # Extract "gene_id" and "seqnames" (chromosome name) columns from the gene GRanges object
  setnames(gene_location_map, "seqnames", "chromosome") # Rename "seqnames" column as "chromosome"

  # Merge the plot data with the location map by "gene_id"
  plot_data_final <- merge(plot_data, gene_location_map, by = "gene_id")


  return(plot_data_final)
}

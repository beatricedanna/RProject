# OKAY
# To execute the function:
# gene_vector <- TotalExpressionPerGene(results)

#' TotalExpressionPerGene
#'
#' @param gene_expression_dt A data.table with gene ids (and gene expression) as rows and cells as columns. The first column should be 'feature_id'. Previously assigned from ExtractData as gene_expression_dt <- results$gene_expression_dt .
#' @return A named numeric vector: total accessibility per peak region.
#'
#' @export
TotalExpressionPerGene <- function(gene_expression_dt) {
  gene_sums_vector <- rowSums(gene_expression_dt[, -1]) # Creates a vector made
  names(gene_sums_vector) <- gene_expression_dt$feature_id

  return(gene_sums_vector)
}

# OKAY
# To execute the function:
# results <- ExtractData(combined_dt)

# Recall gene_expression_dt and atac_seq_dt:
# gene_expression_dt <- results$gene_expression
# atac_seq_dt <- results$atac_seq

# Check:
# print(gene_expression_dt[1:5, 1:5])
# print(atac_seq_dt[1:5, 1:5])

# Declare 'feature_id' as a global variable to avoid NOTES during R CMD check
# This is necessary because 'feature_id' is used as a column name inside data.table expressions
utils::globalVariables("feature_id")

#' ExtractData
#'
#' @importFrom data.table %like%
#' @param combined_dt A combined data table containing feature_id, cell_barcodes and values of both peaks and gene expression.
#'
#' @return A list with two data.tables:
#' \describe{
#'   \item{gene_expression_dt}{Contains rows with `feature_id` starting with "ENSG", representing gene expression data}
#'   \item{atac_seq_dt}{Contains rows with `feature_id` starting with "chr", representing chromatin accessibility data}
#' }
#' @export
ExtractData <- function (combined_dt){
  # Subset creation
  gene_expression_dt <- combined_dt[feature_id %like% "^ENSG"] # Assign to a subset (gene_expression_dt) all values of feature_id that start with ENSG
  atac_seq_dt <- combined_dt[feature_id %like% "^chr"]  # Assign to a subset (atac_seq_dt) all values of feature_id that start with chr

  cat("Gene expression data:", nrow(gene_expression_dt), "genes x", ncol(gene_expression_dt)-1, "cells\n") #ncol-1 because one column is for the ids
  cat("ATAC-seq data:", nrow(atac_seq_dt), "peaks x", ncol(atac_seq_dt)-1, "cells\n") #ncol-1 because one column is for the ids

  return(list(
    gene_expression_dt,
    atac_seq_dt
  ))
}

# To execute the function:
# final_normalized_gex_dt <- NormalizeCPM(gene_expression_dt)
# final_normalized_atac_dt <- NormalizeCPM(atac_seq_dt)

utils::globalVariables("..numeric_col_names")

#' NormalizeCPM
#'
#' @param data_table A data.table containing counts data. It must include a column named "feature_id".
#' @param chunk_size Integer, optional (default = 1000). Number of columns to process at a time.
#'
#' @return A data.table with the same dimensions as data_table, where counts are normalized to CPM and log2-transformed. The "feature_id" column is mantained as the first column.
#' @export
NormalizeCPM <- function(data_table, chunk_size = 1000) {

  id_col_name <- "feature_id" # Create a variable with "feature_id", which is a column name, for easy usage

  feature_ids <- data_table[[id_col_name]] # Extract the column from data_table and assign it to feature_ids variable

  numeric_col_names <- setdiff(names(data_table), id_col_name) # Extract all numeric columns and assign them to numeric_col_names variable

  numeric_dt <- data_table[, ..numeric_col_names] # Create a new version of the data_table containing only numeric columns (to solve 'x must be numeric' error previously receieved)

  cat("Sum of the columns...\n")
  col_sums <- colSums(numeric_dt)

 # Chunking in order to not receive dimension errors
  normalized_chunks_list <- list() # Create an empty list used to store normalized data
  num_cells <- ncol(numeric_dt) # Assign to num_cells the number of column of numeric_dt (it represents the number of cells)

  cat("Processing in chunks of", chunk_size, "cells...\n")
  for (i in seq(1, num_cells, by = chunk_size)) {
    start_col <- i
    end_col <- min(i + chunk_size - 1, num_cells)

    cat("Processing columns from", start_col, "to", end_col, "\n")

    chunk_matrix <- as.matrix(numeric_dt[, start_col:end_col, with = FALSE]) # Extract from the data table columns from start_col to end_col and convert them into a matrix
    current_col_sums <- col_sums[start_col:end_col] # Saves in current_col_sums the sum of the columns processed in this chunk

    matrix_frac <- sweep(chunk_matrix, 2, current_col_sums, FUN = "/") # Divides each column in chunk_matrix by the sum of columns present in current_col_sums
    matrix_norm <- log2((matrix_frac * 10^6) + 1) # Multiply each fraction per 1 million (CPM) and calculates log2

    normalized_chunks_list[[length(normalized_chunks_list) + 1]] <- as.data.table(matrix_norm) # Convert the matrix in a data table and add the result to normalized_chunks_list
  }

  cat("Re-assemble...\n")
  final_normalized_dt <- do.call(cbind, normalized_chunks_list) # Put together all data.table from normalized_chunks_list and join them by column
  final_normalized_dt <- cbind(feature_id = feature_ids, final_normalized_dt) # Add "feature_id" column as the first column

  return(final_normalized_dt)
}

#' SparseToDataTable
#' @importFrom data.table as.data.table
#' @importFrom utils head
#' @param mtx A sparse matrix (e.g., class dgTMatrix or dgCMatrix).
#'
#' @return A data.table with columns: feature_id, cell_barcode, value.
#'
#' @export
SparseToDataTable <- function(mtx){
  full_mat <- as.matrix(mtx)
  data_table <- as.data.table(full_mat,keep.rownames = "row.names")
  print(head(data_table))
  return(data_table)
}

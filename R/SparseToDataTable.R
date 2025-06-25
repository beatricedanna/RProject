# OKAY
# Before running the function:

# Load the libraries
# library(data.table)
# library(Matrix)

# Set working directory
# setwd("~/Desktop/project")

# Load the files
# matrix_file   <- "matrix.mtx.gz"
# features_file <- "features.tsv.gz"
# barcodes_file <- "barcodes.tsv.gz"

# Read the files:
# mtx <- readMM(file = matrix_file)
# features <- fread(features_file, header = FALSE)
# barcodes <- fread(barcodes_file, header = FALSE)

# Check of the dimensions:
# cat("Sparse matrix dimension:", dim(mtx)[1], "features x", dim(mtx)[2], "barcodes\n")

# To execute the function:
# combined_dt <- SparseToDataTable(mtx, features, barcodes)

#' SparseToDataTable
#'
#' @importFrom data.table as.data.table
#' @importFrom utils head
#' @param mtx A sparse matrix (e.g., class dgTMatrix or dgCMatrix).
#' @param features Gene ids or chromosome ids.
#' @param barcodes Cell barcodes for identification.
#'
#' @return A data.table with columns: feature_id, cell_barcode, value.
#'
#' @export
SparseToDataTable <- function (mtx, features, barcodes){
  if (nrow(mtx) == nrow(barcodes)) {
    cat("Matrix found with the wrong orientation -> Transposition in progres...\n")
    mtx <- t(mtx)
  }
  # Assign rownames and colnames
  rownames(mtx) <- features$V1
  colnames(mtx) <- barcodes$V1

  # Conversion in dense matrix
  full_mat <- as.matrix(mtx)

  # Conversion in data.table with row names assigned as feature_id
  combined_dt <- as.data.table(full_mat, keep.rownames = "feature_id")

  return(combined_dt)
}

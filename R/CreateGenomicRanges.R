# OKAY
# Before running the function:

# Load the libraries:
# library(GenomicRanges)
# library(ensembldb)
# library(EnsDb.Hsapiens.v86)

# To execute the function:
# granges_objects <- CreateGenomicRanges(gene_vector, peaks_vector)

# Assign genes_gr and peaks_gr:
# genes_gr <- granges_objects$genes_gr
# peaks_gr <- granges_objects$peaks_gr

# Check:
# class(genes_gr)
# class(peaks_gr)
# as.data.frame(genes_gr[1:5])
# as.data.frame(peaks_gr[1:5])

#' CreateGenomicRanges
#'
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom S4Vectors mcols mcols<-
#'
#' @param gene_vector A named numeric vector: total expression per each gene.
#' @param peaks_vector A named numeric vector: total accessibility per peak region.
#'
#' @return A list containing two GRanges objects:
#' \describe{
#'   \item{genes_gr}{GRanges of genes with total expression as metadata}
#'   \item{peaks_gr}{GRanges of peak regions with total accessibility as metadata}
#' }
#' @export
CreateGenomicRanges <- function(gene_vector, peaks_vector) {

  #CREATING GENE GR

  # Clean gene ids
  names(gene_vector) <- sub("\\..*$", "", names(gene_vector))

  # Load Ensembl Database
  ensdb <- EnsDb.Hsapiens.v86

  # Extract all genes data (contains gene_id column) from EnsDb
  all_genes_gr <- genes(ensdb)

  cat("Finding common genes...\n")
  common_gene_ids <- intersect(names(gene_vector), all_genes_gr$gene_id)
  cat(" Found ", length(common_gene_ids), "genes in  common.\n")

  if (length(common_gene_ids) == 0) {
    stop("No genes in common. ERROR")
  }

  gene_vector_common <- gene_vector[common_gene_ids] # Create a subset of gene_vector keeping only genes in common
  genes_gr_common <- all_genes_gr[all_genes_gr$gene_id %in% common_gene_ids, ] # Create a GRanges object (genes_gr_common) containing only genes present in common_gene_ids
  genes_gr_ordered <- genes_gr_common[match(names(gene_vector_common), genes_gr_common$gene_id), ] # Create an ordered GRanges object following the same order of gene_vector_common, adding NA for each gene not found in common
  mcols(genes_gr_ordered)$total_expression <- gene_vector_common #Add a metadata column of total expression data (originally contained in gene_vector_common)

 # CREATE PEAK GR

  peaks_gr <- GRanges(names(peaks_vector)) #Create a GRanges object with seqnames = chr and ranges = (start,end)
  mcols(peaks_gr)$total_accessibility <- peaks_vector  #Add a metadata column of total accessibility data (originally contained in peaks_vector)

  cat("GenomicRanges Objects created successfully.\n\n")

  return(list(
    genes_gr = genes_gr_ordered,
    peaks_gr = peaks_gr
  ))
}

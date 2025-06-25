# Before running the function:

# Install BiocManager and others:
# install.packages("BiocManager")
# BiocManager::install("GenomeInfoDb")

# Load the libraries:
# library(rtracklayer)
# library(GenomicRanges)
# library(GenomeInfoDb)

# Load the files:
# gtf_file <- system.file("extdata", "Homo_sapiens.GRCh38.114.gtf.gz", package = "RProject")
# full_gtf_gr <- rtracklayer::import(gtf_file)

# To execute the function:
# annotated_peaks_gr_final <- AnnotatePeaks(peaks_gr, full_gtf_gr)

#' AnnotatePeaks
#'
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqlevelsStyle seqlevelsStyle<-
#' @importFrom S4Vectors mcols
#'
#' @param peaks_gr A GRanges object representing peak regions.
#' @param full_gtf_gr A GRanges object imported from a GTF file containing gene annotations.
#'
#' @return A GRanges object identical to peaks_gr but filtered for having only peaks with a nearest protein-coding gene
#' @export
AnnotatePeaks <- function(peaks_gr, full_gtf_gr) {

  original_peak_names <- names(peaks_gr) # Save original chromosome ids

  cat("Filter for coding genes only...\n")
  gene_annotation_gr <- full_gtf_gr[full_gtf_gr$type == "gene" & full_gtf_gr$gene_biotype == "protein_coding"] # Isolation of protein coding genes only containing chromosomal position of these genes (available with seqnames(gene_annotation_gr))

  # Standardization of chromosome names in both peaks_gr and gene_annotation_gr using UCSC style ("chrN")
  seqlevelsStyle(peaks_gr) <- "UCSC"
  seqlevelsStyle(gene_annotation_gr) <- "UCSC"

  cat("Intersection between chromosomes of our peaks data and chromosomes of protein coding genes... \n")
  common_seqlevels <- intersect(seqlevels(peaks_gr), seqlevels(gene_annotation_gr)) # Keep only chromosomes present in both our peak data and in protein coding genes GRanges obj

  # Filter the GRanges obj keeping only genes present in common_seqlevels (removes all genes that do not have a peak in protein coding genes (chromosomes))
  # pruning.mode="coarse" removes all chromosomal regions associated with chromosomes not present anymore in peaks_gr and gene_annotation_gr
  seqlevels(peaks_gr, pruning.mode="coarse") <- common_seqlevels
  seqlevels(gene_annotation_gr, pruning.mode="coarse") <- common_seqlevels

  cat("Finding the nearest gene for each peak...\n")
  nearest_index <- nearest(peaks_gr, gene_annotation_gr) # For each peak in peaks_gr gives the position of the nearest gene, returning a numeric vector of the same lenght of peaks_gr containing the position in gene_annotation_gr of the nearest gene

  valid_indices <- !is.na(nearest_index) # Create a vector returning FALSE for each NA
  peaks_gr <- peaks_gr[valid_indices, ] # Removes all peaks without an associated gene
  nearest_index <- nearest_index[valid_indices] # Updates nearest_index mantaining only valid indices
  final_peak_names <- original_peak_names[valid_indices] # Updates the original list of peak ids keeping only peaks asscoiated with a protein coding gene

  cat("Metadata addition...\n")
  nearest_gene_info <- mcols(gene_annotation_gr[nearest_index]) # Data frame containing metadata from gene_annotation_gr of the genes closest to each peak in peaks_gr
  mcols(peaks_gr)$nearest_gene_id <- nearest_gene_info$gene_id # Addition of nearest_gene_id column to peaks_gr (from nearest_gene_info)
  mcols(peaks_gr)$nearest_gene_name <- nearest_gene_info$gene_name # Addition of nearest_gene_name column to peaks_gr (from nearest_gene_info)

  names(peaks_gr) <- final_peak_names # Substitution of the original peak names with the list of peak names associated with a gene

  return(peaks_gr)
}

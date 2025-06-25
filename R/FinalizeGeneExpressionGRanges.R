# To execute the function:
# protein_coding_gex <- FinalizeGeneExpressionGRanges(genes_gr, full_gtf_gr)

#' FinalizeGeneExpressionGRanges
#'
#' @importFrom S4Vectors mcols
#'
#' @param genes_gr A GRanges object containing gene expression data with `gene_id` metadata.
#' @param full_gtf_gr A GRanges object imported from a full GTF annotation file.
#'
#' @return A GRanges object subset of protein-coding genes, with an added metadata column `gene_name`.
#' @export
FinalizeGeneExpressionGRanges <- function(genes_gr, full_gtf_gr) {

  # Create a reference data frame from the GTF, we need only gene features
  gtf_gene_features <- full_gtf_gr[full_gtf_gr$type == "gene"]

  # Create a simple data.frame for easier lookup. We remove duplicates.
  gene_annotation_df <- as.data.frame(mcols(gtf_gene_features)) # Convert metadata extracted from gtf_gene_features (which is Bioconductor DataFrame) in a data.frame
  gene_annotation_df <- gene_annotation_df[, c("gene_id", "gene_name", "gene_biotype")] # Select only "gene_id", "gene_name" and "gene_biotype" columns
  gene_annotation_df <- unique(gene_annotation_df) # Removes duplicates

  cat("Filtering for protein-coding genes...\n")

  protein_coding_ids <- gene_annotation_df$gene_id[gene_annotation_df$gene_biotype == "protein_coding"] # Get the list of protein-coding gene ids

  # Subset the genes_gr object
  # The gene IDs are already in the 'gene_id' metadata column from a previous step
  protein_coding_expression_gr <- genes_gr[genes_gr$gene_id %in% protein_coding_ids] # Create a subset of genes_gr containing only gene_id of protein_coding_genes (contained in protein_coding_ids)

  cat ( "Total of ", length(protein_coding_expression_gr), "protein-coding genes.\n")

  cat("3. Adding gene IDs...\n")

  matched_indices <- match(
    mcols(protein_coding_expression_gr)$gene_id,
    gene_annotation_df$gene_id
  ) # Create numeric vector matched_indices of the same lenght of protein_coding_expression_gr containing indices giving the position of the protein coding gene in gene_annotation_df$gene_id data.frame

  mcols(protein_coding_expression_gr)$gene_name <- gene_annotation_df$gene_name[matched_indices] # Add "gene_name" column from gene_annotation_df$gene_name based on the indices of matched_indices

  return(protein_coding_expression_gr)
}

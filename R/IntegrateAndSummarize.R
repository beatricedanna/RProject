# Funzione IntegrateAndSummarize - Versione che garantisce la coerenza degli ID

IntegrateAndSummarize <- function(gex_dt, atac_dt, peaks_gr, genes_gr) {

  if (!require("data.table", quietly = TRUE)) stop("Package 'data.table' non trovato.")
  if (!require("GenomicRanges", quietly = TRUE)) stop("Package 'GenomicRanges' non trovato.")

  cat("--- Inizio Integrazione (Garantendo Coerenza ID) ---\n")

  # --- INIZIO CORREZIONE FONDAMENTALE ---

  # 1. Pulisci gli ID Ensembl in TUTTI i data.table di input
  #    Questo assicura che tutti i merge e i confronti successivi funzionino.
  cat("1. Pulizia degli ID Ensembl (rimozione versione)...\n")
  gex_dt$feature_id <- sub("\\..*$", "", gex_dt$feature_id)
  # Nota: atac_dt non ha ID Ensembl, quindi non serve pulirlo.

  # Anche l'oggetto GRanges dei geni finali deve avere ID puliti (anche se dovrebbe già essere così)
  genes_gr$gene_id <- sub("\\..*$", "", genes_gr$gene_id)

  # --- FINE CORREZIONE ---

  # 2. Creazione della mappa Picco-Gene
  cat("2. Creazione della mappa Picco-Gene...\n")
  map_dt <- data.table(
    feature_id = names(peaks_gr),
    nearest_gene_id = mcols(peaks_gr)$nearest_gene_id
  )

  # 3. Merge e Aggregazione (ora lavorano tutti con ID puliti)
  cat("3. Unione e Aggregazione...\n")
  atac_annotated_dt <- merge(atac_dt, map_dt, by = "feature_id", all.x = TRUE)
  atac_annotated_dt <- na.omit(atac_annotated_dt, cols = "nearest_gene_id")
  cell_cols <- grep("^[A-Z]", names(atac_annotated_dt), value = TRUE)
  aggregated_atac_dt <- atac_annotated_dt[, lapply(.SD, sum), by = nearest_gene_id, .SDcols = cell_cols]

  # 4. Merge Finale
  cat("4. Unione finale...\n")
  setnames(aggregated_atac_dt, "nearest_gene_id", "gene_id")
  setnames(gex_dt, "feature_id", "gene_id")
  merged_data <- merge(gex_dt, aggregated_atac_dt, by = "gene_id")

  # 5. Riepiloghi
  cat("5. Creazione delle tabelle di riepilogo...\n")
  merged_gene_ids <- merged_data$gene_id
  unmerged_peak_ids <- map_dt[!nearest_gene_id %in% merged_gene_ids, feature_id]
  unmerged_peaks_summary <- data.table(
    Category = c("Total Peaks", "Merged Peaks", "Unmerged Peaks"),
    Count = c(length(peaks_gr), length(peaks_gr) - length(unmerged_peak_ids), length(unmerged_peak_ids))
  )

  all_expression_gene_ids <- genes_gr$gene_id
  unassociated_gene_ids <- setdiff(all_expression_gene_ids, merged_gene_ids)
  unassociated_genes_summary <- data.table(
    Category = c("Total Genes", "Genes with Peaks", "Genes without Peaks"),
    Count = c(length(all_expression_gene_ids), length(all_expression_gene_ids) - length(unassociated_gene_ids), length(unassociated_gene_ids))
  )

  cat("--- Fine Integrazione e Riepilogo ---\n\n")

  return(list(
    merged_data = merged_data,
    unmerged_peaks_summary = unmerged_peaks_summary,
    unassociated_genes_summary = unassociated_genes_summary,
    unassociated_gene_ids = unassociated_gene_ids
  ))
}

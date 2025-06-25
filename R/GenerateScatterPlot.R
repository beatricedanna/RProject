# Before running the function:

# Load the libraries:
# library(gtools)

# To execute the function:
# scatter_plots <- GenerateScatterPlots(plot_data)

# To print the graphs:
# print(scatter_plots$total_plot)
# print(scatter_plots$faceted_plot)

if (getRversion() >= "2.15.1")  utils::globalVariables(c(
  "avg_gex_signal", "avg_atac_signal", "gene_id", "chromosome"
))

#' GenerateScatterPlots
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth facet_wrap labs theme_bw theme element_rect
#' @importFrom gtools mixedsort
#'
#' @param plot_data A data.table or data.frame containing columns: "avg_gex_signal", "avg_atac_signal", "chromosome".
#' @return A list with two ggplot objects:
#'  \describe{
#'    \item{total_plot}{Overall scatter plot of gene expression vs. accessibility.}
#'    \item{faceted_plot}{Scatter plot faceted by chromosome.}
#' }
#'
#' @export
GenerateScatterPlots <- function(plot_data) {

  cat("Generating the total scatter plot...\n")

  total_plot <- ggplot(plot_data, aes(x = avg_gex_signal, y = avg_atac_signal)) +
    geom_point(alpha = 0.1, size = 0.5) + # Use transparency to reveal density and smaller points for better readabilty
    geom_smooth(method = "lm", color = "red", se = FALSE) + # Add a linear (lm) regression line (se = FALSE -> no confidence interval)
    labs(
      title = "Overall Relationship between Gene Expression and Chromatin Accessibility",
      x = "Average log2(CPM+1) Expression",
      y = "Average Aggregated log2(CPM+1) ATAC Signal"
    ) +
    theme_bw()

  cat("Generating the faceted scatter plot by chromosome...\n")

    chrom_order <- gtools::mixedsort(unique(as.character(plot_data$chromosome))) # Order "naturally" the chromosomes (from "chromosome" column of plot_data) and convert them into characters (removing also duplicates)
    plot_data$chromosome <- factor(plot_data$chromosome, levels = chrom_order) # Transform "chromosome" column in an ordered column (based on chrom_order)

  faceted_plot <- ggplot(plot_data, aes(x = avg_gex_signal, y = avg_atac_signal)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_smooth(method = "lm", color = "blue", se = FALSE, linewidth = 0.5) +
    facet_wrap(~ chromosome, ncol = 6, scales = "free") + # Create small graphs for each chromosome value, ncol = 6 -> displays graphs into 6 columns, scales = "free" -> each graph has independent axis
    labs(
      title = "Gene Expression vs. Chromatin Accessibility by Chromosome",
      x = "Average log2(CPM+1) Expression",
      y = "Average Aggregated ATAC Signal"
    ) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "grey90")) # Style facet labels

  return(list(
    total_plot = total_plot,
    faceted_plot = faceted_plot
  ))
}

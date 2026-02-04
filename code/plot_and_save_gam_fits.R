# Load required libraries
library(tidygam)
library(ggplot2)
library(dplyr)

# Function to plot GAM fitted values along with actual data points

#' Generate and Save GAM Fitted Plots for Multiple Genes
#'
#' @param gam_results Output from run_gam_analysis() containing models
#' @param genes_list Character vector of gene names to plot
#' @param actual_data Data frame with actual gene expression values
#' @param gene_col Name of gene column (default: "gene_name")
#' @param group_col Name of group column (default: "gam_group")
#' @param time_col Name of time column (default: "timePoint_num")
#' @param exp_col Name of expression column (default: "gene_exp")
#' @param output_dir Directory to save plots (default: "gam_plots")
#' @param length_out Number of points for fitted line (default: 20)
#' @param plot_width Plot width in inches (default: 8)
#' @param plot_height Plot height in inches (default: 6)
#' @param time_breaks Numeric vector of time points for x-axis
#' @param time_labels Character vector of labels for time points
#' @param color_palette Color palette (default: MetBrewer::Klimt)
#' @param show_ribbon Show confidence interval ribbon (default: TRUE)
#' @param show_jitter Show jittered actual data points (default: TRUE)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return Data frame of all fitted values and saves plots to output_dir

plot_and_save_gam_fits <- function(gam_results,
                                   genes_list,
                                   actual_data,
                                   gene_col = "gene_name",
                                   group_col = "gam_group",
                                   time_col = "timePoint_num",
                                   exp_col = "gene_exp",
                                   output_dir = "gam_plots",
                                   length_out = 20,
                                   plot_width = 8,
                                   plot_height = 6,
                                   time_breaks = c(0, 0.167, 1, 2),
                                   time_labels = c("0", "4 hr", "24 hr", "48 hr"),
                                   color_palette = NULL,
                                   show_ribbon = TRUE,
                                   show_jitter = TRUE,
                                   verbose = TRUE) {
  

  
  # ============================================
  # Validation
  # ============================================
  
  if (verbose) cat("Validating inputs...\n")
  
  # Check if models exist in gam_results
  if (!"models" %in% names(gam_results)) {
    stop("gam_results must contain 'models' component. Did you set save_models appropriately?")
  }
  
  # Check if genes exist in models
  missing_genes <- setdiff(genes_list, names(gam_results$models))
  if (length(missing_genes) > 0) {
    warning(paste0("The following genes not found in models: ", 
                   paste(missing_genes, collapse = ", ")))
    genes_list <- intersect(genes_list, names(gam_results$models))
  }
  
  if (length(genes_list) == 0) {
    stop("No valid genes to plot")
  }
  
  # Check required columns in actual_data
  required_cols <- c(gene_col, group_col, time_col, exp_col)
  missing_cols <- setdiff(required_cols, names(actual_data))
  if (length(missing_cols) > 0) {
    stop(paste0("Missing columns in actual_data: ", 
                paste(missing_cols, collapse = ", ")))
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat(paste0("✓ Created output directory: ", output_dir, "\n"))
  }
  
  # Set default color palette if not provided
  if (is.null(color_palette)) {
    color_palette <- as.vector(paletteer::paletteer_d("MetBrewer::Klimt"))
  }
  
  if (verbose) {
    cat(paste0("✓ Validation complete\n"))
    cat(paste0("✓ Processing ", length(genes_list), " genes\n"))
    cat(paste0("✓ Output directory: ", output_dir, "\n\n"))
  }
  
  # ============================================
  # Generate fitted values for all genes
  # ============================================
  
  if (verbose) cat("Generating fitted values...\n")
  
  gene_gam_fit_values <- NULL
  
  for (k in seq_along(genes_list)) {
    gene_name <- genes_list[k]
    
    if (verbose && k %% 10 == 0) {
      cat(paste0("  Processing gene ", k, "/", length(genes_list), ": ", gene_name, "\n"))
    }
    
    tryCatch({
      # Extract model for this gene
      gene_model <- gam_results$models[[gene_name]]$overall$full
      
      if (is.null(gene_model)) {
        warning(paste0("Model not found for gene: ", gene_name))
        next
      }
      
      # Predict fitted values
      gene_fit_values_k <- tidygam::predict_gam(gene_model, length_out = length_out) %>%
        mutate(!!sym(gene_col) := gene_name) %>%
        dplyr::rename(fitted_value = gene_exp)
      
      # Combine with all fitted values
      gene_gam_fit_values <- bind_rows(gene_gam_fit_values, gene_fit_values_k)
      
    }, error = function(e) {
      warning(paste0("Error processing gene ", gene_name, ": ", e$message))
    })
  }
  
  if (is.null(gene_gam_fit_values)) {
    stop("No fitted values were generated successfully")
  }
  
  if (verbose) {
    cat(paste0("✓ Generated fitted values for ", 
               length(unique(gene_gam_fit_values[[gene_col]])), 
               " genes\n\n"))
  }
  
  # ============================================
  # Create and save plots for each gene
  # ============================================
  
  if (verbose) cat("Creating and saving plots...\n")
  
  n_plots_saved <- 0
  all_genes_list <- unique(gene_gam_fit_values[[gene_col]])
  for (current_gene in all_genes_list) {
    
    tryCatch({
      # Filter data for this gene
      gene_fitted <- gene_gam_fit_values %>%
        filter(!!sym(gene_col) == current_gene)
      
      gene_actual <- actual_data %>%
        filter(!!sym(gene_col) == current_gene)
      
      # Create plot
      p <- ggplot(gene_fitted) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 0, size = 14),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.title = element_blank(),
          legend.text = element_text(size = 12)
        ) +
        scale_x_continuous(breaks = time_breaks, labels = time_labels) +
        scale_color_manual(values = color_palette) +
        scale_fill_manual(values = color_palette) +
        labs(
          title = current_gene,
          x = "Time point",
          y = "Gene expression"
        )
      
      # Add ribbon if requested
      if (show_ribbon) {
        p <- p + geom_ribbon(
          aes(x = !!sym(time_col), 
              y = fitted_value,
              ymin = lower_ci,
              ymax = upper_ci,
              fill = !!sym(group_col)),
          alpha = 0.3
        )
      }
      
      # Add fitted line
      p <- p + geom_line(
        aes(x = !!sym(time_col),
            y = fitted_value,
            color = !!sym(group_col)),
        linewidth = 2
      )
      
      # Add jittered points if requested
      if (show_jitter) {
        p <- p + geom_jitter(
          data = gene_actual,
          aes(x = !!sym(time_col),
              y = !!sym(exp_col),
              color = !!sym(group_col)),
          alpha = 0.5,
          width = 0.05,
          height = 0
        )
      }
      
      # Remove fill legend (redundant with color)
      p <- p + guides(fill = "none")
      
      # Save plot
      filename <- file.path(output_dir, paste0(current_gene, ".pdf"))
      ggsave(
        filename = filename,
        plot = p,
        width = plot_width,
        height = plot_height,
        dpi = 300
      )
      
      n_plots_saved <- n_plots_saved + 1
      
      if (verbose && n_plots_saved %% 10 == 0) {
        cat(paste0("  Saved ", n_plots_saved, " plots\n"))
      }
      
    }, error = function(e) {
      warning(paste0("Error creating plot for gene ", current_gene, ": ", e$message))
    })
  }
  
  # ============================================
  # Summary
  # ============================================
  
  if (verbose) {
    cat("\n============================================\n")
    cat("Plotting Complete!\n")
    cat("============================================\n")
    cat(paste0("Plots saved: ", n_plots_saved, "/", length(genes_list), "\n"))
    cat(paste0("Output directory: ", output_dir, "\n"))
    cat("============================================\n\n")
  }
  
  # Return fitted values data frame
  return(invisible(gene_gam_fit_values))
}



# ============================================
# Example Usage
# ============================================

# Generate and save plots for all genes
# plot_and_save_gam_fits(
#   gam_results = gam_results,
#   genes_list = signf_genes$gene_name,
#   actual_data = exp_all_genes_final,
#   output_dir = file.path(mastCellsDir, "IL33/gam_plots"),  # Your custom folder
#   time_breaks = c(0, 0.167, 1, 2),
#   time_labels = c("0", "4 hr", "24 hr", "48 hr"),
#   plot_width = 8,
#   plot_height = 4,
#   show_ribbon = FALSE
# )



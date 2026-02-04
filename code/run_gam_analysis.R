library(dplyr)
library(purrr)
library(mgcv)
library(future)
library(furrr)
library(progressr)



#' Run GAM Analysis Across All Genes with Overall and Pairwise Comparisons
#'
#' @param data Data frame containing gene expression data (prepared with baseline copies)
#' @param formula_full Formula for full GAM model
#' @param formula_null Formula for null GAM model (nested in full)
#' @param gene_col Name of gene column (default: "gene_name")
#' @param group_col Name of group column (default: "gam_group")
#' @param num_cores Number of cores for parallel processing (default: 1)
#' @param ci_level Confidence level for CIs (default: 0.95)
#' @param min_obs_per_group Minimum observations per group for pairwise comparisons (default: 5)
#' @param save_models Character, which models to save: "none", "overall_full" (default), "unordered_full", "all"
#' @param verbose Print detailed messages (default: TRUE)
#'
#' @return List with: results (data frame), warnings (data frame), summary (list), models (list of models if requested)
run_gam_analysis <- function(data,
                             formula_full,
                             formula_null,
                             gene_col = "gene_name",
                             group_col = "gam_group",
                             num_cores = 1,
                             ci_level = 0.95,
                             min_obs_per_group = 5,
                             save_models = "overall_full",
                             verbose = TRUE) {
  
  start_time <- Sys.time()
  
  # ============================================
  # Pre-flight Validation
  # ============================================
  
  if (verbose) cat("Running pre-flight checks...\n")
  
  # Validate save_models parameter
  valid_save_options <- c("none", "overall_full", "unordered_full", "all")
  if (!save_models %in% valid_save_options) {
    stop(paste0("save_models must be one of: ", 
                paste(valid_save_options, collapse = ", ")))
  }
  
  # Check required columns exist
  if (!gene_col %in% names(data)) {
    stop(paste0("Gene column '", gene_col, "' not found in data"))
  }
  if (!group_col %in% names(data)) {
    stop(paste0("Group column '", group_col, "' not found in data"))
  }
  
  # Check group column is factor
  if (!is.factor(data[[group_col]])) {
    stop(paste0("Column '", group_col, "' must be a factor"))
  }
  
  # Extract variables from formulas
  extract_vars <- function(formula) {
    all.vars(formula)
  }
  
  full_vars <- extract_vars(formula_full)
  null_vars <- extract_vars(formula_null)
  
  # Check all formula variables exist in data
  missing_full <- setdiff(full_vars, names(data))
  if (length(missing_full) > 0) {
    stop(paste0("Variables in formula_full not found in data: ", 
                paste(missing_full, collapse = ", ")))
  }
  
  missing_null <- setdiff(null_vars, names(data))
  if (length(missing_null) > 0) {
    stop(paste0("Variables in formula_null not found in data: ", 
                paste(missing_null, collapse = ", ")))
  }
  
  # Get group levels
  group_levels <- levels(data[[group_col]])
  n_groups <- length(group_levels)
  
  if (verbose) {
    cat(paste0("✓ Data validation passed\n"))
    cat(paste0("✓ Found ", n_groups, " groups: ", 
               paste(group_levels, collapse = ", "), "\n"))
  }
  
  # Generate all pairwise combinations
  pairwise_combos <- combn(group_levels, 2, simplify = FALSE)
  n_pairs <- length(pairwise_combos)
  
  if (verbose) {
    cat(paste0("✓ Will perform ", n_pairs, " pairwise comparisons\n"))
    cat(paste0("✓ Note: Each comparison fits 2 models (unordered for fixed, ordered for smooth)\n"))
    if (save_models != "none") {
      cat(paste0("✓ Models to save: ", save_models, "\n"))
    }
  }
  
  # Get unique genes
  all_genes <- unique(data[[gene_col]])
  n_genes <- length(all_genes)
  
  if (verbose) {
    cat(paste0("✓ Processing ", n_genes, " genes\n"))
    cat(paste0("✓ Using ", num_cores, " cores\n\n"))
  }
  
  # ============================================
  # Set up parallel processing
  # ============================================
  
  plan(multisession, workers = num_cores)
  
  # ============================================
  # Process single gene function
  # ============================================
  
  process_single_gene <- function(gene, p) {
    
    # Update progress
    p()
    
    # Initialize result row
    result <- data.frame(gene_name = gene)
    warnings_list <- list()
    
    # Initialize model storage based on save_models option
    gene_models <- NULL
    if (save_models != "none") {
      gene_models <- list(gene_name = gene)
      if (save_models %in% c("overall_full", "all")) {
        gene_models$overall <- list()
      }
      if (save_models %in% c("unordered_full", "all")) {
        gene_models$pairwise <- list()
      }
    }
    
    # Subset data for this gene
    gene_data <- data %>% filter(!!sym(gene_col) == gene)
    
    # Check for zero variance
    if (sd(gene_data$gene_exp, na.rm = TRUE) == 0) {
      warnings_list[[length(warnings_list) + 1]] <- data.frame(
        gene_name = gene,
        warning_type = "zero_variance",
        message = "Gene has zero variance, skipping"
      )
      result$overall_fixed_p <- NA
      result$overall_smooth_p <- NA
      return(list(result = result, 
                  warnings = do.call(rbind, warnings_list),
                  models = gene_models))
    }
    
    # ============================================
    # Overall Tests (unordered factors)
    # ============================================
    
    overall_success <- FALSE
    
    tryCatch({
      # Fit full model
      gmt_full <- gam(formula_full, data = gene_data)
      
      # Fit null model
      gmt_null <- gam(formula_null, data = gene_data)
      
      # Save overall models (if save_models == "overall_full" or "all")
      if (save_models %in% c("overall_full", "all")) {
        gene_models$overall$full <- gmt_full
        gene_models$overall$null <- gmt_null
      }
      
      # Overall fixed effect p-value (parametric term)
      gmt_summary <- summary(gmt_full)
      if (!is.null(gmt_summary$p.table) && nrow(gmt_summary$p.table) > 0) {
        # Get p-value for group effect (typically row with group_col name)
        group_rows <- grep(group_col, rownames(gmt_summary$p.table))
        if (length(group_rows) > 0) {
          # Use F-test or chi-square test for overall group effect
          overall_anova <- anova(gmt_full)
          if (!is.null(overall_anova$pTerms.table)) {
            result$overall_fixed_p <- overall_anova$pTerms.pv[1]
          } else {
            result$overall_fixed_p <- NA
          }
        } else {
          result$overall_fixed_p <- NA
        }
      } else {
        result$overall_fixed_p <- NA
      }
      
      # Overall smooth effect p-value (LRT)
      lrt_test <- anova(gmt_full, gmt_null, test = "Chisq")
      result$overall_smooth_p <- lrt_test$`Pr(>Chi)`[2]
      
      overall_success <- TRUE
      
    }, error = function(e) {
      warnings_list[[length(warnings_list) + 1]] <<- data.frame(
        gene_name = gene,
        warning_type = "overall_model_failed",
        message = as.character(e$message)
      )
      result$overall_fixed_p <<- NA
      result$overall_smooth_p <<- NA
    })
    
    # ============================================
    # Pairwise Comparisons (TWO MODELS)
    # ============================================
    
    for (i in seq_along(pairwise_combos)) {
      pair <- pairwise_combos[[i]]
      group1 <- pair[1]
      group2 <- pair[2]
      pair_name <- paste0(group1, "_vs_", group2)
      
      # Initialize pairwise columns
      result[[paste0(pair_name, "_fixed_p")]] <- NA
      result[[paste0(pair_name, "_fixed_estimate")]] <- NA
      result[[paste0(pair_name, "_fixed_se")]] <- NA
      result[[paste0(pair_name, "_fixed_ci_lower")]] <- NA
      result[[paste0(pair_name, "_fixed_ci_upper")]] <- NA
      result[[paste0(pair_name, "_smooth_p")]] <- NA
      
      # Subset to two groups
      pair_data_base <- gene_data %>%
        filter(!!sym(group_col) %in% c(group1, group2))
      
      # Check minimum observations
      group_counts <- pair_data_base %>%
        group_by(!!sym(group_col)) %>%
        summarise(n = n(), .groups = "drop")
      
      if (any(group_counts$n < min_obs_per_group)) {
        warnings_list[[length(warnings_list) + 1]] <- data.frame(
          gene_name = gene,
          warning_type = paste0("pairwise_", pair_name),
          message = paste0("Insufficient observations (min: ", 
                           min(group_counts$n), ")")
        )
        next
      }
      
      # ============================================
      # MODEL 1: Unordered factor (for fixed effects)
      # ============================================
      
      tryCatch({
        # Create unordered factor data
        pair_data_unordered <- pair_data_base %>%
          mutate(!!sym(group_col) := factor(!!sym(group_col), 
                                            levels = c(group1, group2),
                                            ordered = FALSE))
        
        # Fit models with unordered factor
        gmt_pair_full_unord <- gam(formula_full, data = pair_data_unordered)
        
        # Save pairwise unordered model
        if (save_models == "unordered_full") {
          # Store only unordered full model
          if (is.null(gene_models$pairwise)) {
            gene_models$pairwise <- list()
          }
          gene_models$pairwise[[pair_name]] <- gmt_pair_full_unord
        } else if (save_models == "all") {
          # Store all models with full details
          if (is.null(gene_models$pairwise[[pair_name]])) {
            gene_models$pairwise[[pair_name]] <- list()
          }
          gene_models$pairwise[[pair_name]]$unordered <- gmt_pair_full_unord
          gene_models$pairwise[[pair_name]]$data_unordered <- pair_data_unordered
        }
        
        # Extract fixed effect results
        pair_summary <- summary(gmt_pair_full_unord)
        
        if (!is.null(pair_summary$p.table) && nrow(pair_summary$p.table) > 1) {
          # Second row is the group difference (first group is reference)
          coef_idx <- 2
          
          result[[paste0(pair_name, "_fixed_estimate")]] <- 
            pair_summary$p.table[coef_idx, "Estimate"]
          
          result[[paste0(pair_name, "_fixed_se")]] <- 
            pair_summary$p.table[coef_idx, "Std. Error"]
          
          result[[paste0(pair_name, "_fixed_p")]] <- 
            pair_summary$p.table[coef_idx, "Pr(>|t|)"]
          
          # Calculate CI
          z_crit <- qnorm((1 + ci_level) / 2)
          est <- pair_summary$p.table[coef_idx, "Estimate"]
          se <- pair_summary$p.table[coef_idx, "Std. Error"]
          
          result[[paste0(pair_name, "_fixed_ci_lower")]] <- est - z_crit * se
          result[[paste0(pair_name, "_fixed_ci_upper")]] <- est + z_crit * se
        }
        
      }, error = function(e) {
        warnings_list[[length(warnings_list) + 1]] <<- data.frame(
          gene_name = gene,
          warning_type = paste0("pairwise_", pair_name, "_unordered"),
          message = as.character(e$message)
        )
      })
      
      # ============================================
      # MODEL 2: Ordered factor (for smooth effects)
      # ============================================
      
      tryCatch({
        # Create ordered factor data
        pair_data_ordered <- pair_data_base %>%
          mutate(!!sym(group_col) := factor(!!sym(group_col), 
                                            levels = c(group1, group2),
                                            ordered = TRUE))
        
        # Fit model with ordered factor
        gmt_pair_full_ord <- gam(formula_full, data = pair_data_ordered)
        
        # Save pairwise ordered model (only if save_models == "all")
        if (save_models == "all") {
          gene_models$pairwise[[pair_name]]$ordered <- gmt_pair_full_ord
          gene_models$pairwise[[pair_name]]$data_ordered <- pair_data_ordered
        }
        
        # Extract smooth effect p-value from summary
        pair_summary_ord <- summary(gmt_pair_full_ord)
        
        if (!is.null(pair_summary_ord$s.table) && nrow(pair_summary_ord$s.table) > 0) {
          # Look for the by-group smooth term (contains group_col name)
          smooth_rows <- grep(group_col, rownames(pair_summary_ord$s.table))
          if (length(smooth_rows) > 0) {
            # Get p-value from the by-group smooth term
            result[[paste0(pair_name, "_smooth_p")]] <- 
              pair_summary_ord$s.table[smooth_rows[1], "p-value"]
          }
        }
        
      }, error = function(e) {
        warnings_list[[length(warnings_list) + 1]] <<- data.frame(
          gene_name = gene,
          warning_type = paste0("pairwise_", pair_name, "_ordered"),
          message = as.character(e$message)
        )
      })
    }
    
    # Combine warnings
    warnings_df <- if (length(warnings_list) > 0) {
      do.call(rbind, warnings_list)
    } else {
      NULL
    }
    
    return(list(result = result, 
                warnings = warnings_df,
                models = gene_models))
  }
  
  # ============================================
  # Run analysis with progress bar
  # ============================================
  
  if (verbose) cat("Starting analysis...\n\n")
  
  with_progress({
    p <- progressor(steps = n_genes)
    
    results_list <- future_map(all_genes, 
                               ~process_single_gene(.x, p),
                               .options = furrr_options(seed = TRUE))
  })
  
  # ============================================
  # Combine results
  # ============================================
  
  if (verbose) cat("\nCombining results...\n")
  
  results_df <- map_dfr(results_list, "result")
  warnings_df <- map(results_list, "warnings") %>%
    discard(is.null) %>%
    bind_rows()
  
  # Extract models if saved
  models_list <- NULL
  if (save_models != "none") {
    models_list <- map(results_list, "models") %>%
      set_names(all_genes)
  }
  
  # ============================================
  # Summary statistics
  # ============================================
  
  end_time <- Sys.time()
  time_elapsed <- difftime(end_time, start_time, units = "mins")
  
  n_successful <- sum(!is.na(results_df$overall_fixed_p))
  n_failed <- sum(is.na(results_df$overall_fixed_p))
  n_warnings <- nrow(warnings_df)
  
  summary_stats <- list(
    n_genes_total = n_genes,
    n_genes_successful = n_successful,
    n_genes_failed = n_failed,
    n_warnings = n_warnings,
    n_groups = n_groups,
    n_pairwise_comparisons = n_pairs,
    n_models_per_gene = 2 + (n_pairs * 2),  # overall(2) + pairwise(unord_full + ord_full)*n_pairs
    models_saved = save_models,
    time_elapsed_mins = as.numeric(time_elapsed),
    ci_level = ci_level,
    min_obs_per_group = min_obs_per_group
  )
  
  if (verbose) {
    cat("\n============================================\n")
    cat("Analysis Complete!\n")
    cat("============================================\n")
    cat(paste0("Total genes: ", n_genes, "\n"))
    cat(paste0("Successfully processed: ", n_successful, "\n"))
    cat(paste0("Failed: ", n_failed, "\n"))
    cat(paste0("Warnings generated: ", n_warnings, "\n"))
    cat(paste0("Models per gene: ", summary_stats$n_models_per_gene, "\n"))
    cat(paste0("Models saved: ", save_models, "\n"))
    cat(paste0("Time elapsed: ", round(time_elapsed, 2), " minutes\n"))
    cat("============================================\n\n")
  }
  
  # Close parallel backend
  plan(sequential)
  
  # Return results
  output <- list(
    results = results_df,
    warnings = warnings_df,
    summary = summary_stats
  )
  
  if (save_models != "none") {
    output$models <- models_list
  }
  
  return(output)
}






# ============================================
# Example Usage
# ============================================

# Define formulas
# formula_full <- formula("gene_exp ~ gam_group + s(timePoint_num, k=3) + 
#                          s(timePoint_num, k=3, by=gam_group) + 
#                          s(donor_id, bs='re')")
# 
# formula_null <- formula("gene_exp ~ gam_group + s(timePoint_num, k=3) + 
#                          s(donor_id, bs='re')")
# 
# # Run analysis
# gam_results <- run_gam_analysis(
#   data = exp_all_genes_final,
#   formula_full = formula_full,
#   formula_null = formula_null,
#   gene_col = "gene_name",
#   group_col = "gam_group",
#   num_cores = 4,
#   ci_level = 0.95,
#   min_obs_per_group = 5,
#   verbose = TRUE
# )
# 
# # Access results
# results_df <- gam_results$results
# warnings_df <- gam_results$warnings
# summary_stats <- gam_results$summary
# 


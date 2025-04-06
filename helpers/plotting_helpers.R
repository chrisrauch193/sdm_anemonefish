# helpers/plotting_helpers.R

library(ggplot2)
library(corrplot)
library(car) # Needed for car::vif in plot_vif_results usage context


# Define the technical variable names and their corresponding display names
# This lookup table should contain ALL potential variables you might use across scenarios.
display_names_lookup <- c(
  "par_baseline_depthsurf_mean"  = "PAR Surface Mean",
  "sws_baseline_depthsurf_mean"  = "Wind Stress Surf Mean", # Shortened slightly
  "thetao_baseline_depthmax_mean"= "Temp MaxDepth Mean",
  "thetao_baseline_depthmax_range"= "Temp MaxDepth Range",
  "so_baseline_depthmax_mean"    = "Salinity MaxDepth Mean",
  "no3_baseline_depthmax_mean"   = "Nitrate MaxDepth Mean",
  "no3_baseline_depthmax_range"  = "Nitrate MaxDepth Range",
  "chl_baseline_depthmax_mean"   = "Chlorophyll MaxDepth Mean",
  "o2_baseline_depthmax_range"   = "Oxygen MaxDepth Range",
  "ph_baseline_depthmax_mean"    = "PH MaxDepth Mean",
  
  "par_ssp119_depthsurf_dec50_mean"  = "PAR Surface Mean",
  "sws_ssp119_depthsurf_dec50_mean"  = "Wind Stress Surf Mean", # Shortened slightly
  "thetao_ssp119_depthmax_dec50_mean"= "Temp MaxDepth Mean",
  "thetao_ssp119_depthmax_dec50_range"= "Temp MaxDepth Range",
  "so_ssp119_depthmax_dec50_mean"    = "Salinity MaxDepth Mean",
  "no3_ssp119_depthmax_dec50_mean"   = "Nitrate MaxDepth Mean",
  "no3_ssp119_depthmax_dec50_range"  = "Nitrate MaxDepth Range",
  "chl_ssp119_depthmax_dec50_mean"   = "Chlorophyll MaxDepth Mean",
  "o2_ssp119_depthmax_dec50_range"   = "Oxygen MaxDepth Range",
  "ph_ssp119_depthmax_dec50_mean"    = "PH MaxDepth Mean",
  
  "par_ssp119_depthsurf_dec100_mean"  = "PAR Surface Mean",
  "sws_ssp119_depthsurf_dec100_mean"  = "Wind Stress Surf Mean", # Shortened slightly
  "thetao_ssp119_depthmax_dec100_mean"= "Temp MaxDepth Mean",
  "thetao_ssp119_depthmax_dec100_range"= "Temp MaxDepth Range",
  "so_ssp119_depthmax_dec100_mean"    = "Salinity MaxDepth Mean",
  "no3_ssp119_depthmax_dec100_mean"   = "Nitrate MaxDepth Mean",
  "no3_ssp119_depthmax_dec100_range"  = "Nitrate MaxDepth Range",
  "chl_ssp119_depthmax_dec100_mean"   = "Chlorophyll MaxDepth Mean",
  "o2_ssp119_depthmax_dec100_range"   = "Oxygen MaxDepth Range",
  "ph_ssp119_depthmax_dec100_mean"    = "PH MaxDepth Mean",
  
  "par_ssp585_depthsurf_dec50_mean"  = "PAR Surface Mean",
  "sws_ssp585_depthsurf_dec50_mean"  = "Wind Stress Surf Mean", # Shortened slightly
  "thetao_ssp585_depthmax_dec50_mean"= "Temp MaxDepth Mean",
  "thetao_ssp585_depthmax_dec50_range"= "Temp MaxDepth Range",
  "so_ssp585_depthmax_dec50_mean"    = "Salinity MaxDepth Mean",
  "no3_ssp585_depthmax_dec50_mean"   = "Nitrate MaxDepth Mean",
  "no3_ssp585_depthmax_dec50_range"  = "Nitrate MaxDepth Range",
  "chl_ssp585_depthmax_dec50_mean"   = "Chlorophyll MaxDepth Mean",
  "o2_ssp585_depthmax_dec50_range"   = "Oxygen MaxDepth Range",
  "ph_ssp585_depthmax_dec50_mean"    = "PH MaxDepth Mean",
  
  "par_ssp585_depthsurf_dec100_mean"  = "PAR Surface Mean",
  "sws_ssp585_depthsurf_dec100_mean"  = "Wind Stress Surf Mean", # Shortened slightly
  "thetao_ssp585_depthmax_dec100_mean"= "Temp MaxDepth Mean",
  "thetao_ssp585_depthmax_dec100_range"= "Temp MaxDepth Range",
  "so_ssp585_depthmax_dec100_mean"    = "Salinity MaxDepth Mean",
  "no3_ssp585_depthmax_dec100_mean"   = "Nitrate MaxDepth Mean",
  "no3_ssp585_depthmax_dec100_range"  = "Nitrate MaxDepth Range",
  "chl_ssp585_depthmax_dec100_mean"   = "Chlorophyll MaxDepth Mean",
  "o2_ssp585_depthmax_dec100_range"   = "Oxygen MaxDepth Range",
  "ph_ssp585_depthmax_dec100_mean"    = "PH MaxDepth Mean",
  
  "bathymetry_mean"              = "Bathymetry Mean",
  "distcoast"                    = "Distance to Coast",
  "rugosity"                     = "Rugosity"
  # Add other potential variables here if they exist in your env data
)


plot_vif_results <- function(vif_result, save_path = NULL) {
  # Ensure vif_result has names (technical variable names)
  if (is.null(names(vif_result))) {
    stop("vif_result must have names (technical variable names).")
  }
  
  # Check if vif_result comes directly from car::vif (single vector)
  # or potentially from usdm::vif (data frame - handle common cases)
  if (is.data.frame(vif_result) && "VIF" %in% colnames(vif_result) && "Variables" %in% colnames(vif_result)) {
    df <- data.frame(Variable = vif_result$Variables, VIF = as.numeric(vif_result$VIF))
  } else if (is.numeric(vif_result) && !is.null(names(vif_result))) {
    df <- data.frame(Variable = names(vif_result), VIF = as.numeric(vif_result)) # Ensure VIF is numeric
  } else {
    stop("vif_result format not recognized. Expected named vector (like from car::vif) or data.frame with 'Variables' and 'VIF' columns (like from usdm::vifstep/vifcor).")
  }
  
  
  # Get display names, use technical name if lookup fails (NA)
  original_vars <- df$Variable
  display_labels <- display_names_lookup[original_vars]
  display_labels[is.na(display_labels)] <- original_vars[is.na(display_labels)] # Fallback
  
  # Replace technical names with display names in the data frame
  # Use factor to maintain the original order in the plot
  df$Variable <- factor(display_labels, levels = display_labels[order(match(original_vars, names(vif_result)))]) # Preserve original order
  
  # Define the color palette based on the number of variables
  num_vars <- nrow(df)
  colors <- colorRampPalette(c("#1F3F8C", "#A9192A"))(num_vars) # Dark blue to dark red
  
  # Ensure VIF values are non-negative
  df$VIF <- pmax(0, df$VIF)
  
  # ggplot now uses the display names stored in df$Variable
  p <- ggplot(df, aes(x = Variable, y = VIF)) +
    geom_bar(stat = "identity", aes(fill = Variable), show.legend = FALSE) +
    # Assign colors based on the factor levels (which are the display names)
    scale_fill_manual(values = setNames(colors, levels(df$Variable))) +
    labs(title = "VIF Analysis Results",
         x     = "Environment Variables",
         y     = "VIF Values") +
    scale_y_continuous(limits = c(0, ceiling(max(df$VIF, 10, na.rm = TRUE))), # Ensure limit is at least 10, handle NAs
                       breaks = seq(0, ceiling(max(df$VIF, 10, na.rm = TRUE)), by = 2)) + # Adjust breaks
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") + # Add VIF threshold line
    theme_minimal() +
    # Adjust text size if labels are long
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) # Slightly smaller text
  
  if (!is.null(save_path)) {
    # Adjust saved plot dimensions if needed for longer labels
    ggsave(save_path, plot = p, width = 9, height = 7) # Slightly wider/taller
    cat("VIF plot saved to", save_path, "\n")
  }
  
  return(p)
}


plot_correlation_results <- function(env_extract, save_path = NULL) {
  
  # Check if env_extract has columns
  if (ncol(env_extract) == 0) {
    warning("env_extract has no columns. Cannot generate correlation plot.")
    return(NULL) # Or handle appropriately
  }
  if (nrow(env_extract) < 2) {
    warning("env_extract has less than 2 rows. Cannot calculate correlations.")
    return(NULL)
  }
  # Remove columns with zero variance if any exist before correlation
  variances <- apply(env_extract, 2, var, na.rm = TRUE)
  zero_var_cols <- variances == 0 | is.na(variances) # Also check for NA variance
  if (any(zero_var_cols)) {
    removed_cols <- names(variances[zero_var_cols])
    warning("Removing columns with zero or NA variance before correlation: ",
            paste(removed_cols, collapse=", "))
    env_extract <- env_extract[, !zero_var_cols, drop = FALSE]
    
    # Update display_names_lookup if columns were removed
    # valid_vars <- colnames(env_extract)
    # display_names_lookup <- display_names_lookup[names(display_names_lookup) %in% valid_vars]
    
    if (ncol(env_extract) < 2) {
      warning("Not enough columns left after removing zero/NA variance columns. Cannot generate correlation plot.")
      return(NULL)
    }
  }
  
  
  # Pearson correlation analysis
  env.cor <- round(cor(env_extract, method = "pearson", use = "pairwise.complete.obs"), 3)
  
  # Alternative p-value calculation (robust to NAs from pairwise deletion)
  n_vars <- ncol(env_extract)
  env.p <- matrix(NA, nrow = n_vars, ncol = n_vars)
  colnames(env.p) <- rownames(env.p) <- colnames(env_extract)
  for (i in 1:n_vars) {
    for (j in i:n_vars) {
      if (i == j) {
        env.p[i, j] <- 0
      } else {
        # Only run cor.test if there are enough complete pairs
        complete_pairs <- sum(complete.cases(env_extract[, c(i, j)]))
        if(complete_pairs > 2) {
          test_result <- tryCatch(cor.test(env_extract[[i]], env_extract[[j]], method = "pearson"), error = function(e) NULL)
          if (!is.null(test_result)) {
            env.p[i, j] <- env.p[j, i] <- test_result$p.value
          }
        } else {
          # Not enough pairs to calculate p-value, leave as NA
          env.p[i, j] <- env.p[j, i] <- NA
        }
      }
    }
  }
  env.p <- round(env.p, 3)
  
  
  # --- Map technical names to display names for matrix dimensions ---
  technical_names <- colnames(env.cor) # Or rownames
  if (is.null(technical_names)) {
    warning("Correlation matrix has no column names. Cannot apply display names.")
    return(NULL)
  }
  # Ensure display_names_lookup only contains names present in technical_names
  display_names_lookup_subset <- display_names_lookup[names(display_names_lookup) %in% technical_names]
  # Get display labels, fallback to technical names if not found
  display_labels <- display_names_lookup_subset[technical_names]
  display_labels[is.na(display_labels)] <- technical_names[is.na(display_labels)]
  
  # Apply the display names to the matrices
  colnames(env.cor) <- rownames(env.cor) <- display_labels
  if (!is.null(env.p)) {
    colnames(env.p)   <- rownames(env.p)   <- display_labels
  }
  # --- End of Mapping ---
  
  # Define the plotting commands within a function for cleaner saving/displaying
  plot_commands <- function() {
    # Adjust cex (character expansion) values if labels are long
    corrplot(corr = env.cor, type = "upper", tl.pos = "tp",
             tl.col = "black", tl.cex = 0.7, cl.cex = 0.8,
             # Only add p.mat if it exists and has the right dimensions
             p.mat = if (!is.null(env.p) && all(dim(env.p) == dim(env.cor))) env.p else NULL,
             insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 1, pch.col = "grey", # Adjust significance levels/display
             order = "original") # Using display names now
    corrplot(
      corr = env.cor, type = "lower", add = TRUE, method = "number",
      tl.pos = "n", # No text labels on lower half
      col = "black", # Color of numbers
      diag = FALSE, cl.pos = "n", # No color legend on lower half
      number.cex = 0.7, number.font = 1, order = "original", # Reduced number.cex
      # Add correlation coefficients significance stars to lower plot numbers
      addCoef.col = "black", # Color of the coefficients
      # Only add p.mat if it exists and has the right dimensions
      p.mat = if (!is.null(env.p) && all(dim(env.p) == dim(env.cor))) env.p else NULL,
      sig.level = c(.001, .01, .05), insig = "label_sig", pch.col="transparent", pch.cex=0.01 # Make symbols transparent
    )
  }
  
  # Saving the plot
  if (!is.null(save_path)) {
    png(filename = save_path, width = 9, height = 9, units = "in", res = 300) # Increased size for labels
    op <- par(mar = c(1, 1, 3, 1)) # Adjust margins: bottom, left, top, right
    plot_commands() # Execute plotting commands to the PNG device
    title(main = paste("Correlation Matrix:", basename(dirname(save_path))), line = 1.5, adj = 0.5) # Add title
    par(op) # Restore original par settings
    dev.off()
    cat("Correlation plot saved to", save_path, "\n")
  }
  
  # Plotting to the current device (optional, can be slow)
  # op <- par(mfrow = c(1, 1), mar = c(1, 1, 3, 1)) # Set margins for screen display
  # plot_commands() # Plot to current device (e.g., RStudio Plots pane)
  # title(main = paste("Correlation Matrix:", basename(dirname(save_path))), line = 1.5, adj=0.5)
  # # recorded_plot <- recordPlot() # Capture it if needed
  # par(op) # Restore original par settings
  
  # Return the correlation matrix itself
  return(env.cor)
}
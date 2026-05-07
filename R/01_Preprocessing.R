# ==============================================================================
# Module 1: Data Preparation, Cleaning, and Preprocessing (Advanced Flow Cytometry Edition)
# Filename: R/01_Preprocessing.R
# ==============================================================================

#' Automatically clean FCS data (QC: debris removal + doublet removal)
Auto_Clean_FCS <- function(fcs) {
  # 1. Remove debris (based on the 0.5% lower quantile cutoff of FSC-A)
  if ("FSC-A" %in% colnames(fcs)) {
    fsc_vals <- flowCore::exprs(fcs)[, "FSC-A"]
    threshold <- quantile(fsc_vals, probs = 0.005, na.rm = TRUE)
    fcs <- fcs[fsc_vals > threshold, ]
  }
  
  # 2. Remove doublets (Singlet Gating)
  cols <- colnames(fcs)
  if ("FSC-A" %in% cols && "FSC-H" %in% cols) {
    tryCatch({
      sg <- flowStats::singletGate(fcs, area = "FSC-A", height = "FSC-H", prediction_level = 0.99)
      fcs <- flowCore::Subset(fcs, sg)
    }, error = function(e) { warning("      Automatic doublet removal failed; skipping this step.") })
  }
  return(fcs)
}

#' Read and process FCS files at a professional level
#' @param do_QC Whether to perform doublet removal and dead-cell/debris filtering
#' @param do_compensation Whether to automatically detect and apply the internal fluorescence compensation matrix in the FCS file
#' @param transform_method Choose the data transformation method: "arcsinh" (common flow-cytometry standard), "logicle" (dynamic biexponential transform), or "none"
#' @param cofactor Cofactor for ArcSinh transformation; typically 150 for flow cytometry and 5 for CyTOF
Prepare_Raw_Data <- function(work_dir, group_keywords, n_cells_per_sample = 5000, 
                             do_QC = TRUE, 
                             do_compensation = TRUE, 
                             transform_method = c("arcsinh", "logicle", "none"),
                             cofactor = 150) {
  
  transform_method <- match.arg(transform_method)
  files <- list.files(work_dir, pattern = ".fcs$", full.names = TRUE)
  if(length(files) == 0) stop("  No FCS files found")
  
  combined_data <- data.frame()
  message(">>>   Processing ", length(files), " FCS files...")
  message(sprintf("        Preprocessing settings: QC=%s | compensation=%s | transform=%s", do_QC, do_compensation, transform_method))
  message("---------------------------------------------------------------")
  
  for (file in files) {
    # Read raw data (disable built-in transformations and preserve the original matrix)
    fcs <- flowCore::read.FCS(file, transformation = FALSE, truncate_max_range = FALSE)
    n_original <- nrow(fcs)
    
    # ==========================================
    #    Step 1: Quality Control (QC)
    # ==========================================
    if (do_QC) {
      fcs <- Auto_Clean_FCS(fcs)
      n_clean <- nrow(fcs)
      percent_removed <- round((n_original - n_clean) / n_original * 100, 1)
      message(sprintf("     %s\n      Original: %d -> QC Cleaned: %d |    Debris removed: %s%%",
                      basename(file), n_original, n_clean, percent_removed))
    } else {
      message(sprintf("     %s |   Loaded: %d cells", basename(file), n_original))
    }
    
    # ==========================================
    #   Step 2: Compensation (fluorescence compensation)
    # ==========================================
    if (do_compensation) {
      # Automatically detect the spillover matrix keyword written by the instrument
      spill_matrix <- flowCore::keyword(fcs)$SPILL
      if (is.null(spill_matrix)) spill_matrix <- flowCore::keyword(fcs)$spillover
      
      if (!is.null(spill_matrix)) {
        fcs <- flowCore::compensate(fcs, spill_matrix)
      } else {
        warning("         No embedded compensation matrix (SPILL/spillover) detected; skipping compensation.")
      }
    }
    
    # ==========================================
    #   Step 3: Data Transformation (nonlinear transformation)
    # ==========================================
    # Automatically identify fluorescence channels to transform (excluding physical/time channels such as FSC, SSC, and Time)
    fluor_cols <- colnames(fcs)[
      !grepl("Time|TIME|FSC|SSC", colnames(fcs), ignore.case = TRUE)
    ]    
    if (transform_method == "arcsinh") {
      #   Fix: the correct function name in flowCore is arcsinhTransform
      trans_list <- flowCore::transformList(fluor_cols, flowCore::arcsinhTransform(a = 0, b = 1/cofactor, c = 0))
      fcs <- flowCore::transform(fcs, trans_list)
      
    } else if (transform_method == "logicle") {
      # Logicle is a dynamic biexponential transform that handles negative tails more effectively
      tryCatch({
        lgcl_trans <- flowCore::estimateLogicle(fcs, channels = fluor_cols)
        fcs <- flowCore::transform(fcs, lgcl_trans)
      }, error = function(e) {
        warning("         Logicle dynamic estimation failed (possibly due to extreme outliers); automatically falling back to ArcSinh transformation.")
        #   Fix: also corrected the spelling in the fallback mechanism
        trans_list <- flowCore::transformList(fluor_cols, flowCore::arcsinhTransform(a = 0, b = 1/cofactor, c = 0))
        fcs <- flowCore::transform(fcs, trans_list)
      })
    }
    
    # ==========================================
    #    Step 4: Extract matrix and downsample
    # ==========================================
    exprs_data <- as.data.frame(flowCore::exprs(fcs))
    
    # Extract valid columns (by convention, keep fluorescence area channels ending in "-A")
    use_cols <- grep("-A$", colnames(exprs_data), value = TRUE)
    if(length(use_cols) == 0) use_cols <- colnames(exprs_data) # Handle instruments that do not use the -A suffix
    use_cols <- use_cols[!grepl("Time|TIME", use_cols)]
    
    data_clean <- exprs_data[, use_cols, drop = FALSE]
    
    # Random downsampling
    if (nrow(data_clean) > n_cells_per_sample) {
      data_sub <- data_clean[sample(nrow(data_clean), n_cells_per_sample), ]
    } else {
      data_sub <- data_clean
    }
    
    # Add group metadata labels
    data_sub$SampleID <- basename(file)
    data_sub$Group <- "Other"
    for (key in group_keywords) {
      if (grepl(key, basename(file), ignore.case = TRUE)) { data_sub$Group <- key; break }
    }
    
    combined_data <- dplyr::bind_rows(combined_data, data_sub)
  }
  
  row.names(combined_data) <- paste0("Cell_", 1:nrow(combined_data))
  
  message("---------------------------------------------------------------")
  message(">>>   Data loading and standardized preprocessing completed! (Real biological differences were preserved)")
  return(combined_data)
}

#' Apply channel-to-protein name mapping
Apply_Panel_Mapping <- function(df, panel_list) {
  message(">>>   Applying panel channel mapping...")
  current_cols <- colnames(df)
  found_count <- 0
  for (channel in names(panel_list)) {
    if (channel %in% current_cols) {
      colnames(df)[current_cols == channel] <- panel_list[[channel]]
      found_count <- found_count + 1
    }
  }
  message(paste0(">>>   Successfully replaced ", found_count, " channel names!"))
  return(df)
}


#' Automatically filter the large cell population using the exported threshold template (multidimensional bounding-box trimming)
#' @param df Large raw data frame to be filtered
#' @param target_node Tree node name to use as the global gate (for example, "CD45" or "Live_CD45")
#' @param template_path Path to the template file (default: "Gating_Thresholds_Template.csv")
#' @return Filtered clean data frame
Apply_Gating_Template <- function(df, target_node, template_path = "Gating_Thresholds_Template.csv", margin = 0) {

  if (!file.exists(template_path)) {
    warning(paste("   Template file not found:", template_path, ", template-based filtering will be skipped."))
    return(df)
  }

  gating_tpl <- read.csv(template_path)
  gate_rules <- gating_tpl[gating_tpl$Tree_Node == target_node, ]

  if (nrow(gate_rules) == 0) {
    warning(sprintf("   Node named '%s' not found in the template; skipping filtering. Please check that the spelling matches the gating step.", target_node))
    return(df)
  }

  message(sprintf("\n>>>    Applying template rule [%s] to stringently filter the large cell population...", target_node))
  if (margin != 0) message(sprintf("   margin = %.4f applied to tighten gate boundaries", margin))
  original_n <- nrow(df)

  # Traverse the threshold rules for this node across all channels and apply multidimensional bounding-box trimming
  for (i in 1:nrow(gate_rules)) {
    chn <- gate_rules$Channel[i]
    min_val <- gate_rules$Min_Threshold[i] + margin
    max_val <- gate_rules$Max_Threshold[i] - margin

    if (min_val >= max_val) {
      warning(sprintf("   margin=%.4f is too large for channel %s (min=%.3f >= max=%.3f); skipping this channel.", margin, chn, min_val, max_val))
      next
    }

    if (chn %in% colnames(df)) {
      df <- df[df[[chn]] >= min_val & df[[chn]] <= max_val, ]
    }
  }

  retention_rate <- (nrow(df) / original_n) * 100
  message(sprintf(">>>   Filtering completed! Debris successfully removed; clean cells remaining: %d (retention: %.1f%%)", nrow(df), retention_rate))

  return(df)
}

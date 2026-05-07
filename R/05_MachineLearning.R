# ==============================================================================
# Module 5: Machine Learning Engine   Automated Gating Projection (Industrial-Grade ML Projection)
# Filename: R/05_MachineLearning.R
# New: panel-compatibility checks, cv_folds=1 bug fix, and lightweight Project_New_Sample version
# ==============================================================================

# ==============================================================================
# ADTnorm-inspired landmark alignment for flow/CITE-like protein batches
# ==============================================================================

.EFML_default_ignore_cols <- function() {
  c("tSNE1", "tSNE2", "UMAP1", "UMAP2", "Cluster", "SampleID", "Group",
    "Time", "TIME", "CellType", "Pseudotime", "Cluster_Name", "Gate_Status",
    "Batch", "id")
}

.EFML_quantile_landmarks <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) < 10 || length(unique(x)) < 2) return(numeric(0))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8))
}

.EFML_density_landmarks <- function(x,
                                    density_adjust = 1,
                                    grid_n = 512,
                                    min_peak_height = 0.05,
                                    min_peak_distance_frac = 0.03,
                                    max_internal_landmarks = 7) {
  x <- x[is.finite(x)]
  if (length(x) < 50 || length(unique(x)) < 5) return(numeric(0))
  
  den <- tryCatch(
    stats::density(x, n = grid_n, adjust = density_adjust, na.rm = TRUE),
    error = function(e) NULL
  )
  if (is.null(den) || length(den$x) < 5) return(numeric(0))
  
  y <- den$y
  xx <- den$x
  n <- length(y)
  peak_idx <- which(y[2:(n - 1)] >= y[1:(n - 2)] & y[2:(n - 1)] > y[3:n]) + 1
  if (length(peak_idx) == 0) return(numeric(0))
  
  peak_idx <- peak_idx[y[peak_idx] >= max(y, na.rm = TRUE) * min_peak_height]
  if (length(peak_idx) == 0) return(numeric(0))
  
  min_dist <- diff(range(xx, na.rm = TRUE)) * min_peak_distance_frac
  ranked <- peak_idx[order(y[peak_idx], decreasing = TRUE)]
  kept <- integer(0)
  for (idx in ranked) {
    if (length(kept) == 0 || all(abs(xx[idx] - xx[kept]) >= min_dist)) {
      kept <- c(kept, idx)
    }
  }
  peak_idx <- sort(kept)
  
  valley_x <- numeric(0)
  if (length(peak_idx) >= 2) {
    for (i in seq_len(length(peak_idx) - 1)) {
      span <- peak_idx[i]:peak_idx[i + 1]
      valley_x <- c(valley_x, xx[span[which.min(y[span])]])
    }
  } else {
    right_tail <- which(xx > xx[peak_idx])
    if (length(right_tail) > 5) {
      valley_x <- xx[right_tail[which.min(y[right_tail])]]
    }
  }
  
  internal <- sort(unique(c(xx[peak_idx], valley_x)))
  if (length(internal) > max_internal_landmarks) {
    keep <- unique(round(seq(1, length(internal), length.out = max_internal_landmarks)))
    internal <- internal[keep]
  }
  internal
}

.EFML_prepare_mapping <- function(source_landmarks, target_landmarks) {
  if (length(source_landmarks) != length(target_landmarks)) return(NULL)
  keep <- is.finite(source_landmarks) & is.finite(target_landmarks)
  source_landmarks <- source_landmarks[keep]
  target_landmarks <- target_landmarks[keep]
  if (length(source_landmarks) < 2) return(NULL)
  
  ord <- order(source_landmarks)
  source_landmarks <- source_landmarks[ord]
  target_landmarks <- target_landmarks[ord]
  
  keep_unique <- !duplicated(source_landmarks)
  source_landmarks <- source_landmarks[keep_unique]
  target_landmarks <- target_landmarks[keep_unique]
  if (length(source_landmarks) < 2) return(NULL)
  
  list(source = source_landmarks, target = target_landmarks)
}

.EFML_warp_values <- function(x, mapping) {
  out <- x
  finite_idx <- is.finite(x)
  out[finite_idx] <- stats::approx(
    x = mapping$source,
    y = mapping$target,
    xout = x[finite_idx],
    rule = 2,
    ties = "ordered"
  )$y
  out
}

.EFML_sample_for_density <- function(x, max_n = 50000) {
  x <- x[is.finite(x)]
  if (length(x) > max_n) x <- sample(x, max_n)
  x
}

.EFML_add_density_lines <- function(df, feature, batch_col, batches, max_cells_per_batch = 50000) {
  pal <- grDevices::rainbow(length(batches))
  has_plot <- FALSE
  
  for (i in seq_along(batches)) {
    idx <- as.character(df[[batch_col]]) == batches[i]
    idx[is.na(idx)] <- FALSE
    vals <- .EFML_sample_for_density(df[idx, feature], max_cells_per_batch)
    if (length(vals) < 10 || length(unique(vals)) < 2) next
    den <- tryCatch(stats::density(vals, na.rm = TRUE), error = function(e) NULL)
    if (is.null(den)) next
    
    if (!has_plot) {
      plot(den, col = pal[i], lwd = 2, main = feature, xlab = "Expression", ylab = "Density")
      has_plot <- TRUE
    } else {
      lines(den, col = pal[i], lwd = 2)
    }
  }
  
  if (has_plot) {
    legend("topright", legend = batches, col = pal, lwd = 2, cex = 0.7, bty = "n")
  } else {
    plot.new()
    title(main = paste(feature, "(insufficient variation)"))
  }
}

.EFML_feature_medians <- function(df, features, fallback_df = NULL) {
  vals <- vapply(features, function(feature) {
    x <- df[[feature]]
    med <- suppressWarnings(stats::median(x[is.finite(x)], na.rm = TRUE))
    if (!is.finite(med) && !is.null(fallback_df) && feature %in% colnames(fallback_df)) {
      y <- fallback_df[[feature]]
      med <- suppressWarnings(stats::median(y[is.finite(y)], na.rm = TRUE))
    }
    if (!is.finite(med)) med <- 0
    med
  }, numeric(1))
  vals
}

.EFML_impute_nonfinite_features <- function(df, features, impute_values, context = "data") {
  report <- vector("list", length(features))
  for (i in seq_along(features)) {
    feature <- features[i]
    bad_idx <- !is.finite(df[[feature]])
    n_bad <- sum(bad_idx)
    if (n_bad > 0) df[[feature]][bad_idx] <- impute_values[[feature]]
    report[[i]] <- data.frame(
      Feature = feature,
      Context = context,
      Imputed_Cells = n_bad,
      Impute_Value = impute_values[[feature]],
      stringsAsFactors = FALSE
    )
  }
  report <- do.call(rbind, report)
  list(data = df, report = report[report$Imputed_Cells > 0, , drop = FALSE])
}

.EFML_safe_scale_features <- function(df, features, context = "data") {
  mat <- as.data.frame(df[, features, drop = FALSE])
  centers <- vapply(mat, function(x) {
    val <- mean(x[is.finite(x)], na.rm = TRUE)
    if (!is.finite(val)) 0 else val
  }, numeric(1))
  scales <- vapply(mat, function(x) {
    val <- stats::sd(x[is.finite(x)], na.rm = TRUE)
    if (!is.finite(val) || val == 0) 1 else val
  }, numeric(1))
  
  scaled <- sweep(as.matrix(mat), 2, centers, "-")
  scaled <- sweep(scaled, 2, scales, "/")
  nonfinite_n <- sum(!is.finite(scaled))
  if (nonfinite_n > 0) {
    scaled[!is.finite(scaled)] <- 0
    message(sprintf("         Safe scaling replaced %d non-finite values in %s with 0.", nonfinite_n, context))
  }
  as.data.frame(scaled)
}

#' ADTnorm-inspired batch-effect removal for protein marker matrices
#'
#' This lightweight implementation follows the core idea described in ADTnorm:
#' detect marker-wise density landmarks and align them across batches using a
#' monotone one-to-one transformation. If peak/valley landmarks are ambiguous,
#' the function falls back to robust quantile landmarks.
#'
#' @param data Cell-level data frame containing marker columns and a batch column
#' @param batch_col Column that identifies experimental batch
#' @param reference_batch Batch used as the expression-scale reference
#' @param features Marker columns to correct. If NULL, numeric non-metadata columns are used
#' @param quantile_probs Quantiles used when density landmarks cannot be matched
#' @param density_adjust Bandwidth adjustment passed to density()
#' @param diagnostic_pdf Optional PDF file with before/after density diagnostics
#' @return A corrected data.frame with ADTnorm_Lite_Info stored as an attribute
Remove_Batch_Effects_ADTnorm_Lite <- function(data,
                                              batch_col = "Batch",
                                              reference_batch = NULL,
                                              features = NULL,
                                              quantile_probs = c(0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999),
                                              density_adjust = 1,
                                              min_peak_height = 0.05,
                                              min_peak_distance_frac = 0.03,
                                              max_internal_landmarks = 7,
                                              diagnostic_pdf = NULL,
                                              diagnostic_csv = "ADTnorm_Lite_Batch_Correction_Diagnostics.csv",
                                              max_cells_per_density = 50000) {
  if (!batch_col %in% colnames(data)) {
    stop(sprintf("Batch column '%s' was not found in data.", batch_col))
  }
  
  batches <- unique(as.character(data[[batch_col]]))
  batches <- batches[!is.na(batches)]
  if (length(batches) < 2) {
    message(">>> ADTnorm-lite skipped: fewer than two batches were detected.")
    return(data)
  }
  
  if (is.null(reference_batch)) reference_batch <- batches[1]
  if (!reference_batch %in% batches) {
    stop(sprintf("reference_batch '%s' was not found in %s.", reference_batch, batch_col))
  }
  
  if (is.null(features)) {
    numeric_cols <- colnames(data)[sapply(data, is.numeric)]
    features <- setdiff(numeric_cols, .EFML_default_ignore_cols())
  }
  features <- intersect(features, colnames(data))
  features <- features[sapply(data[, features, drop = FALSE], is.numeric)]
  if (length(features) == 0) stop("No numeric marker features were available for ADTnorm-lite correction.")
  
  message("\n>>> [ADTnorm-lite] Removing marker-wise batch effects by landmark alignment...")
  message(sprintf("     - Batch column: %s", batch_col))
  message(sprintf("     - Reference batch: %s", reference_batch))
  message(sprintf("     - Number of corrected markers: %d", length(features)))
  
  original_data <- data
  corrected <- data
  diag_rows <- list()
  
  if (!is.null(diagnostic_pdf)) {
    grDevices::pdf(diagnostic_pdf, width = 11, height = 5.5)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  for (feature in features) {
    ref_vals <- original_data[as.character(original_data[[batch_col]]) == reference_batch, feature]
    ref_quant <- .EFML_quantile_landmarks(ref_vals, quantile_probs)
    ref_density <- .EFML_density_landmarks(
      ref_vals,
      density_adjust = density_adjust,
      min_peak_height = min_peak_height,
      min_peak_distance_frac = min_peak_distance_frac,
      max_internal_landmarks = max_internal_landmarks
    )
    
    if (length(ref_quant) < 2) {
      diag_rows[[length(diag_rows) + 1]] <- data.frame(
        Feature = feature, Batch = reference_batch, Method = "skipped_reference_low_variation",
        Ref_Landmarks = length(ref_density), Batch_Landmarks = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    for (batch in setdiff(batches, reference_batch)) {
      idx <- as.character(original_data[[batch_col]]) == batch
      batch_vals <- original_data[idx, feature]
      batch_quant <- .EFML_quantile_landmarks(batch_vals, quantile_probs)
      batch_density <- .EFML_density_landmarks(
        batch_vals,
        density_adjust = density_adjust,
        min_peak_height = min_peak_height,
        min_peak_distance_frac = min_peak_distance_frac,
        max_internal_landmarks = max_internal_landmarks
      )
      
      method <- "quantile_fallback"
      mapping <- NULL
      
      if (length(ref_density) >= 2 && length(ref_density) == length(batch_density)) {
        ref_anchor <- sort(c(ref_quant[c(1, length(ref_quant))], ref_density))
        batch_anchor <- sort(c(batch_quant[c(1, length(batch_quant))], batch_density))
        mapping <- .EFML_prepare_mapping(batch_anchor, ref_anchor)
        if (!is.null(mapping)) method <- "density_peak_valley"
      }
      
      if (is.null(mapping) && length(batch_quant) == length(ref_quant)) {
        mapping <- .EFML_prepare_mapping(batch_quant, ref_quant)
      }
      
      if (is.null(mapping)) {
        method <- "skipped_low_variation"
      } else {
        corrected[idx, feature] <- .EFML_warp_values(batch_vals, mapping)
      }
      
      diag_rows[[length(diag_rows) + 1]] <- data.frame(
        Feature = feature,
        Batch = batch,
        Method = method,
        Ref_Landmarks = length(ref_density),
        Batch_Landmarks = length(batch_density),
        stringsAsFactors = FALSE
      )
    }
    
    if (!is.null(diagnostic_pdf)) {
      old_par <- graphics::par(no.readonly = TRUE)
      graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
      .EFML_add_density_lines(original_data, feature, batch_col, batches, max_cells_per_density)
      graphics::title(main = paste0(feature, " - before correction"))
      .EFML_add_density_lines(corrected, feature, batch_col, batches, max_cells_per_density)
      graphics::title(main = paste0(feature, " - after correction"))
      graphics::par(old_par)
    }
  }
  
  diag_df <- if (length(diag_rows) > 0) do.call(rbind, diag_rows) else data.frame()
  attr(corrected, "ADTnorm_Lite_Info") <- diag_df
  
  if (!is.null(diagnostic_csv) && nrow(diag_df) > 0) {
    utils::write.csv(diag_df, diagnostic_csv, row.names = FALSE)
    message(sprintf("     - Diagnostic table saved to: %s", diagnostic_csv))
  }
  if (!is.null(diagnostic_pdf)) {
    message(sprintf("     - Diagnostic density PDF saved to: %s", diagnostic_pdf))
  }
  
  message(">>> ADTnorm-lite batch correction completed.")
  corrected
}

#' Automatic projection engine for ultra-large single-cell atlases
#' @param ref_obj Reference object (containing annotated Data and UMAP_Model)
#' @param massive_data Large raw dataset to be predicted
#' @param target_label Label column to learn and project (for example, "Cluster_Name" or "CellType")
#' @param cv_folds Number of cross-validation folds (default 5; set to 1 to skip CV)
#' @param rejection_threshold Rejection threshold (0~1). Cells with prediction probabilities below this threshold are labeled "Unassigned"
#' @param use_class_weights Whether to use class weights to balance rare cell types (default TRUE)
#' @param rare_cell_boost Boost factor for rare cell types (default 2.0, higher = more protection for rare types)
#' @param batch_correction "none" or "adtnorm_lite"; use "adtnorm_lite" to align marker distributions across batches before projection
#' @param batch_col Batch column used for ADTnorm-lite correction
#' @param reference_batch Batch used as the landmark-alignment reference
Project_Atlas_To_Massive_Data <- function(ref_obj,
                                          massive_data,
                                          target_label = "Cluster_Name",
                                          cv_folds = 5,
                                          rejection_threshold = 0.5,
                                          use_class_weights = TRUE,
                                          rare_cell_boost = 2.0,
                                          batch_correction = c("none", "adtnorm_lite"),
                                          batch_col = NULL,
                                          reference_batch = NULL,
                                          batch_diagnostic_pdf = NULL) {

  if (!requireNamespace("ranger", quietly = TRUE)) stop("  Missing the 'ranger' random-forest package")
  if (!requireNamespace("uwot",   quietly = TRUE)) stop("  Missing the 'uwot' package; UMAP projection cannot be performed")
  batch_correction <- match.arg(batch_correction)

  ref_data <- ref_obj$Data
  if (!(target_label %in% colnames(ref_data))) stop("  Target label column not found in the reference atlas: ", target_label)

  message(sprintf("\n======================================================="))
  message(sprintf(">>> [Ultra-large-scale ML projection] Starting industrial-grade projection engine"))
  message(sprintf("     - Target number of large-scale cells: %d", nrow(massive_data)))
  message(sprintf("     - Target label to predict: %s", target_label))
  message(sprintf("     - OOD rejection threshold (OOD Rejection): %.2f", rejection_threshold))
  if (use_class_weights) {
    message(sprintf("     - Class weighting: ENABLED (rare_cell_boost=%.1f)", rare_cell_boost))
  } else {
    message(sprintf("     - Class weighting: DISABLED"))
  }
  message(sprintf("=======================================================\n"))

  ignore_cols <- .EFML_default_ignore_cols()
  numeric_cols <- colnames(ref_data)[sapply(ref_data, is.numeric)]
  features     <- setdiff(numeric_cols, ignore_cols)

  #   Change 1: panel-compatibility checks to prevent silent failure
  missing_features <- setdiff(features, colnames(massive_data))
  if (length(missing_features) > 0) {
    stop(sprintf("  Panel mismatch: massive_data is missing the following feature columns:\n   %s\nPlease check that the channel mapping is consistent.",
                 paste(missing_features, collapse = ", ")))
  }
  
  batch_info <- NULL
  if (batch_correction == "adtnorm_lite") {
    if (is.null(batch_col)) stop("batch_col must be provided when batch_correction = 'adtnorm_lite'.")
    massive_data <- Remove_Batch_Effects_ADTnorm_Lite(
      data = massive_data,
      batch_col = batch_col,
      reference_batch = reference_batch,
      features = features,
      diagnostic_pdf = batch_diagnostic_pdf
    )
    batch_info <- attr(massive_data, "ADTnorm_Lite_Info")
  }
  
  feature_impute_values <- .EFML_feature_medians(ref_data, features, fallback_df = massive_data)
  massive_impute <- .EFML_impute_nonfinite_features(
    massive_data,
    features,
    feature_impute_values,
    context = "massive_data"
  )
  massive_data <- massive_impute$data
  imputation_info <- massive_impute$report
  if (nrow(imputation_info) > 0) {
    message("   [Feature imputation] Non-finite marker values were found and replaced with reference medians:")
    for (i in seq_len(nrow(imputation_info))) {
      message(sprintf("         %s: %d cells -> %.4f",
                      imputation_info$Feature[i],
                      imputation_info$Imputed_Cells[i],
                      imputation_info$Impute_Value[i]))
    }
  }

  train_df <- ref_data[!is.na(ref_data[[target_label]]) & ref_data[[target_label]] != "Ungated",
                       c(features, target_label)]
  train_impute <- .EFML_impute_nonfinite_features(
    train_df,
    features,
    feature_impute_values,
    context = "reference_training_data"
  )
  train_df <- train_impute$data
  if (nrow(train_impute$report) > 0) {
    imputation_info <- rbind(imputation_info, train_impute$report)
  }
  train_df[[target_label]] <- as.factor(train_df[[target_label]])

  #   NEW: Calculate class weights to protect rare cell types
  class_weights <- NULL
  if (use_class_weights) {
    class_counts <- table(train_df[[target_label]])
    class_props <- class_counts / sum(class_counts)

    # Calculate inverse frequency weights with boost for rare types
    # Formula: weight = (1 / proportion) ^ rare_cell_boost
    # This gives exponentially higher weights to rarer classes
    class_weights <- (1 / class_props) ^ rare_cell_boost

    # Normalize weights so mean = 1
    class_weights <- class_weights / mean(class_weights)

    # Assign weights to each sample
    sample_weights <- class_weights[as.character(train_df[[target_label]])]

    message("   [Class weight summary]")
    weight_summary <- data.frame(
      CellType = as.character(names(class_counts)),
      Count = as.numeric(class_counts),
      Proportion = round(as.numeric(class_props) * 100, 2),
      Weight = round(as.numeric(class_weights), 2),
      stringsAsFactors = FALSE
    )
    weight_summary <- weight_summary[order(weight_summary$Count), ]

    # Show top 5 rarest types with highest weights
    message("         Top 5 rarest types (highest protection):")
    for (i in 1:min(5, nrow(weight_summary))) {
      message(sprintf("           %s: n=%d (%.2f%%), weight=%.2fx",
                      weight_summary$CellType[i],
                      weight_summary$Count[i],
                      weight_summary$Proportion[i],
                      weight_summary$Weight[i]))
    }
  }
  
  #   Change 2: fix the bug where cv_accuracies was not initialized when cv_folds=1
  cv_accuracies <- numeric(0)  # Initialize as an empty vector so mean() returns NaN rather than an error for empty input

  if (cv_folds > 1) {
    message(sprintf("   [1/4] Performing %d-fold cross-validation model evaluation...", cv_folds))
    set.seed(42)
    folds         <- sample(rep(1:cv_folds, length.out = nrow(train_df)))
    cv_accuracies <- numeric(cv_folds)

    for (i in 1:cv_folds) {
      fold_train <- train_df[folds != i, ]
      fold_test  <- train_df[folds == i, ]

      # Apply class weights in CV if enabled
      if (use_class_weights) {
        fold_weights <- sample_weights[folds != i]
        fold_model <- ranger::ranger(
          dependent.variable.name = target_label, data = fold_train,
          num.trees = 50, classification = TRUE,
          case.weights = fold_weights,
          num.threads = parallel::detectCores() - 1
        )
      } else {
        fold_model <- ranger::ranger(
          dependent.variable.name = target_label, data = fold_train,
          num.trees = 50, classification = TRUE,
          num.threads = parallel::detectCores() - 1
        )
      }

      fold_pred       <- predict(fold_model, data = fold_test[, features])$predictions
      cv_accuracies[i] <- sum(fold_pred == fold_test[[target_label]]) / nrow(fold_test)
    }
    message(sprintf("         CV mean validation accuracy: %.2f%% (+-%.2f%%)",
                    mean(cv_accuracies) * 100, sd(cv_accuracies) * 100))
  } else {
    message("   [1/4] cv_folds=1 Skipping cross-validation.")
  }

  message("   [2/4] Training the final decision boundary on the full reference atlas...")
  set.seed(42)

  # Train final model with or without class weights
  if (use_class_weights) {
    rf_model <- ranger::ranger(
      dependent.variable.name = target_label, data = train_df,
      num.trees = 100, probability = TRUE, importance = "impurity",
      case.weights = sample_weights,
      num.threads = parallel::detectCores() - 1
    )
  } else {
    rf_model <- ranger::ranger(
      dependent.variable.name = target_label, data = train_df,
      num.trees = 100, probability = TRUE, importance = "impurity",
      num.threads = parallel::detectCores() - 1
    )
  }
  
  importance_df <- data.frame(
    Marker     = names(rf_model$variable.importance),
    Importance = as.numeric(rf_model$variable.importance)
  )
  importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]
  top_5 <- head(importance_df, 5)
  message("         [Model decision basis] Top 5 top driving features (Feature Importance):")
  for (k in 1:nrow(top_5)) {
    message(sprintf("             %d. %s (Score: %.1f)", k, top_5$Marker[k], top_5$Importance[k]))
  }
  
  message("   [3/4] Projecting labels onto the large-scale cells and performing OOD rejection...")
  pred <- predict(rf_model, data = massive_data[, features])
  
  prob_matrix       <- pred$predictions
  max_probs         <- apply(prob_matrix, 1, max)
  predicted_classes <- colnames(prob_matrix)[max.col(prob_matrix, ties.method = "first")]
  
  rejected_idx <- max_probs < rejection_threshold
  predicted_classes[rejected_idx] <- "Unassigned"
  
  n_rejected <- sum(rejected_idx)
  message(sprintf("         Rejection completed: %d cells with confidence < %.2f were assigned to 'Unassigned' (%.1f%%)",
                  n_rejected, rejection_threshold, (n_rejected / nrow(massive_data)) * 100))
  
  massive_data[[target_label]]                          <- predicted_classes
  massive_data[[paste0(target_label, "_Confidence")]]   <- round(max_probs, 3)
  
  message("   [4/4] Embedding large-scale cells into the established UMAP manifold...")
  if (!is.null(ref_obj$Info$UMAP_Model)) {
    umap_input <- massive_data[, features]
    
    if (isTRUE(ref_obj$Info$Scaled)) {
      message("         Detected that the reference model was scaled; automatically aligning the Z-score feature space...")
      umap_input <- .EFML_safe_scale_features(umap_input, features, context = "large-scale UMAP input")
    } else {
      umap_nonfinite <- .EFML_impute_nonfinite_features(
        umap_input,
        features,
        feature_impute_values,
        context = "large-scale UMAP input"
      )
      umap_input <- umap_nonfinite$data
    }
    if (any(!is.finite(as.matrix(umap_input[, features, drop = FALSE])))) {
      stop("UMAP input still contains non-finite values after imputation and safe scaling.")
    }
    
    umap_proj        <- uwot::umap_transform(X = as.matrix(umap_input),
                                             model = ref_obj$Info$UMAP_Model, verbose = FALSE)
    massive_data$UMAP1 <- umap_proj[, 1]
    massive_data$UMAP2 <- umap_proj[, 2]
  } else {
    warning("   Reference atlas lacks the underlying UMAP model; skipping coordinate projection.")
  }
  
  #   Change 2 continued: when cv_accuracies is empty, mean() returns NaN; replace it with NA for clearer semantics
  cv_mean_acc <- if (length(cv_accuracies) > 0) mean(cv_accuracies) else NA_real_
  
  message(">>> Ultra-large-scale projection completed!")
  return(list(
    Data = massive_data,
    Info = list(
      Method       = "Massive_Atlas_Projection_with_Rejection",
      Model_CV_Acc = cv_mean_acc,
      Features     = features,
      RF_Model     = rf_model,
      Batch_Correction = list(
        Method = batch_correction,
        Batch_Column = batch_col,
        Reference_Batch = reference_batch,
        Diagnostics = batch_info
      ),
      Feature_Imputation = imputation_info,
      Feature_Impute_Values = feature_impute_values
    )
  ))
}


#'   New: lightweight single-sample fast projection (no CV, suitable for routine rapid prediction)
#' @param ref_obj Trained reference object (must contain Info$RF_Model and Info$Features)
#' @param new_data New sample data frame to predict
#' @param target_label Predicted label column name
#' @param rejection_threshold Rejection threshold
#' @details This function reuses the RF model trained by Project_Atlas_To_Massive_Data,
#'   skips retraining and CV, and is suitable for fast batch projection of multiple new samples.
#'   Before use, make sure ref_obj is the object returned by Project_Atlas_To_Massive_Data
#'   (with RF_Model and Features stored in Info).
Project_New_Sample <- function(ref_obj,
                               new_data,
                               target_label = "Cluster_Name",
                               rejection_threshold = 0.5) {
  
  if (!requireNamespace("ranger", quietly = TRUE)) stop("  Missing the 'ranger' random-forest package")
  if (!requireNamespace("uwot",   quietly = TRUE)) stop("  Missing the 'uwot' package; UMAP projection cannot be performed")
  
  rf_model <- ref_obj$Info$RF_Model
  features <- ref_obj$Info$Features
  feature_impute_values <- ref_obj$Info$Feature_Impute_Values
  
  if (is.null(rf_model)) stop("  ref_obj$Info$RF_Model is empty; please run Project_Atlas_To_Massive_Data() first to generate the model.")
  if (is.null(features)) stop("  ref_obj$Info$Features is empty; unable to determine the feature columns.")
  
  # Panel compatibility check
  missing_features <- setdiff(features, colnames(new_data))
  if (length(missing_features) > 0) {
    stop(sprintf("  Panel mismatch: new_data is missing the following feature columns:\n   %s",
                 paste(missing_features, collapse = ", ")))
  }
  if (is.null(feature_impute_values)) {
    feature_impute_values <- .EFML_feature_medians(new_data, features)
  }
  new_impute <- .EFML_impute_nonfinite_features(
    new_data,
    features,
    feature_impute_values,
    context = "new_sample_data"
  )
  new_data <- new_impute$data
  if (nrow(new_impute$report) > 0) {
    message("     Non-finite values in the new sample were replaced before projection.")
  }
  
  message(sprintf("\n>>> [Lightweight projection] Rapidly projecting %d cells...", nrow(new_data)))
  
  pred        <- predict(rf_model, data = new_data[, features])
  prob_matrix <- pred$predictions
  max_probs   <- apply(prob_matrix, 1, max)
  pred_labels <- colnames(prob_matrix)[max.col(prob_matrix, ties.method = "first")]
  
  rejected_idx <- max_probs < rejection_threshold
  pred_labels[rejected_idx] <- "Unassigned"
  
  n_rejected <- sum(rejected_idx)
  message(sprintf("     Rejected: %d cells (%.1f%%) were assigned to 'Unassigned'",
                  n_rejected, n_rejected / nrow(new_data) * 100))
  
  new_data[[target_label]]                        <- pred_labels
  new_data[[paste0(target_label, "_Confidence")]] <- round(max_probs, 3)
  
  # UMAP projection
  if (!is.null(ref_obj$Info$UMAP_Model)) {
    umap_input <- new_data[, features]
    if (isTRUE(ref_obj$Info$Scaled)) {
      umap_input <- .EFML_safe_scale_features(umap_input, features, context = "new-sample UMAP input")
    } else {
      umap_nonfinite <- .EFML_impute_nonfinite_features(
        umap_input,
        features,
        feature_impute_values,
        context = "new-sample UMAP input"
      )
      umap_input <- umap_nonfinite$data
    }
    if (any(!is.finite(as.matrix(umap_input[, features, drop = FALSE])))) {
      stop("New-sample UMAP input still contains non-finite values after imputation and safe scaling.")
    }
    umap_proj      <- uwot::umap_transform(X = as.matrix(umap_input),
                                           model = ref_obj$Info$UMAP_Model, verbose = FALSE)
    new_data$UMAP1 <- umap_proj[, 1]
    new_data$UMAP2 <- umap_proj[, 2]
  }
  
  message(">>> Lightweight projection completed!")
  return(new_data)
}

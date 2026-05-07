# ==============================================================================
# Module 7: Industrial-grade statistics engine, trajectory inference, and automated HTML reporting
# Filename: R/07_Stats_and_Export.R
# ==============================================================================

# ==============================================================================
# 1. Differential cell-abundance analysis
# ==============================================================================

#' Calculate cell-subpopulation proportions and perform statistical testing
#' @param method "wilcoxon", "t.test", "kruskal", "anova", "lmm"
#' @param transform "none", "asin_sqrt", "clr"
#' @param random_effect Random-effect column name (for LMM only)
#' @param sample_col Sample column name
Perform_Abundance_Stats <- function(
    res_obj,
    target_label  = "Cluster_Name",
    group_col     = "Group",
    method        = c("wilcoxon", "t.test", "kruskal", "anova", "lmm"),
    transform     = c("none", "asin_sqrt", "clr"),
    random_effect = NULL,
    sample_col    = "SampleID") {
  
  df        <- res_obj$Data
  method    <- match.arg(method)
  transform <- match.arg(transform)
  
  message(sprintf("\n>>> [Statistics engine] Running %s abundance differential analysis...", target_label))
  message(sprintf("     - Test method: %s | Data transformation: %s", toupper(method), toupper(transform)))
  
  if (is.null(df) || !is.data.frame(df)) stop("res_obj$Data does not exist or is not a data.frame.")
  
  required_cols <- c(sample_col, target_label, group_col)
  missing_cols  <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("res_obj$Data is missing required columns: %s\nCurrently available columns: %s",
                 paste(missing_cols, collapse = ", "),
                 paste(colnames(df), collapse = ", ")))
  }
  
  if (!is.null(random_effect) && !random_effect %in% colnames(df)) {
    stop(sprintf("random_effect '%s' does not exist in res_obj$Data.", random_effect))
  }
  
  count_df <- as.data.frame(table(
    SampleID = df[[sample_col]],
    Cluster  = df[[target_label]],
    useNA    = "no"
  ))
  count_df <- count_df[count_df$Freq > 0, , drop = FALSE]
  
  prop_df <- count_df |>
    dplyr::group_by(SampleID) |>
    dplyr::mutate(Proportion = Freq / sum(Freq), Total_Cells = sum(Freq)) |>
    dplyr::ungroup()
  
  meta_cols <- unique(c(sample_col, group_col, random_effect))
  meta_cols <- meta_cols[!is.na(meta_cols) & nzchar(meta_cols)]
  meta_df   <- df[, meta_cols, drop = FALSE] |> dplyr::distinct()
  colnames(meta_df)[colnames(meta_df) == sample_col] <- "SampleID"
  
  dup_check <- meta_df |>
    dplyr::distinct(SampleID, .data[[group_col]]) |>
    dplyr::count(SampleID) |>
    dplyr::filter(n > 1)
  if (nrow(dup_check) > 0) stop("Detected multiple group labels for the same SampleID; please clean the metadata first.")
  
  stat_df <- merge(prop_df, meta_df, by = "SampleID", all.x = TRUE)
  if (any(is.na(stat_df[[group_col]]))) stop(sprintf("Some samples were not matched to group information for '%s'.", group_col))
  
  if (transform == "asin_sqrt") {
    stat_df$Value <- asin(sqrt(stat_df$Proportion))
  } else if (transform == "clr") {
    positive_props <- stat_df$Proportion[stat_df$Proportion > 0]
    if (length(positive_props) == 0) stop("All proportions are 0; CLR transformation cannot be performed.")
    pseudo_count <- min(positive_props) / 2
    stat_df <- stat_df |>
      dplyr::group_by(SampleID) |>
      dplyr::mutate(
        Prop_adj = ifelse(Proportion == 0, pseudo_count, Proportion),
        Value    = log(Prop_adj) - mean(log(Prop_adj))
      ) |>
      dplyr::ungroup()
  } else {
    stat_df$Value <- stat_df$Proportion
  }
  
  clusters <- unique(stat_df$Cluster)
  res_list <- list()
  
  for (cl in clusters) {
    sub_df   <- stat_df[stat_df$Cluster == cl, , drop = FALSE]
    n_groups <- length(unique(stats::na.omit(sub_df[[group_col]])))
    if (n_groups < 2) next
    
    #   Change: assign the tryCatch return value to a local variable to fully fix the closure-scope issue
    cl_result <- tryCatch({
      if (method == "wilcoxon" && n_groups == 2) {
        test_res <- wilcox.test(stats::as.formula(paste("Value ~", group_col)), data = sub_df)
        data.frame(Cluster = cl, p_value = test_res$p.value, Method = "Wilcoxon")
      } else if (method == "t.test" && n_groups == 2) {
        test_res <- t.test(stats::as.formula(paste("Value ~", group_col)), data = sub_df)
        data.frame(Cluster = cl, p_value = test_res$p.value, Method = "T-Test")
      } else if (method == "kruskal" && n_groups >= 2) {
        test_res <- kruskal.test(stats::as.formula(paste("Value ~", group_col)), data = sub_df)
        data.frame(Cluster = cl, p_value = test_res$p.value, Method = "Kruskal-Wallis")
      } else if (method == "anova" && n_groups >= 2) {
        fit      <- aov(stats::as.formula(paste("Value ~", group_col)), data = sub_df)
        test_res <- broom::tidy(fit)
        data.frame(Cluster = cl, p_value = test_res$p.value[1], Method = "ANOVA")
      } else if (method == "lmm") {
        if (is.null(random_effect)) stop("LMM requires the random_effect parameter")
        if (!requireNamespace("lmerTest", quietly = TRUE)) {
          stop("Package 'lmerTest' is required for method='lmm'. Install it with install.packages('lmerTest').")
        }
        if (!requireNamespace("broom.mixed", quietly = TRUE)) {
          stop("Package 'broom.mixed' is required for method='lmm'. Install it with install.packages('broom.mixed').")
        }
        form     <- stats::as.formula(paste("Value ~", group_col, "+ (1|", random_effect, ")"))
        fit      <- getExportedValue("lmerTest", "lmer")(form, data = sub_df)
        test_res <- getExportedValue("broom.mixed", "tidy")(fit)
        p_val    <- test_res$p.value[grepl(paste0("^", group_col), test_res$term)]
        p_val    <- if (length(p_val) > 0) p_val[1] else NA_real_
        data.frame(Cluster = cl, p_value = p_val, Method = "Linear Mixed Model")
      } else {
        NULL
      }
    }, error = function(e) {
      data.frame(Cluster = cl, p_value = NA_real_,
                 Method = paste0("Error/NA: ", conditionMessage(e)))
    })
    
    if (!is.null(cl_result)) res_list[[cl]] <- cl_result
  }
  
  final_res <- do.call(rbind, res_list)
  if (!is.null(final_res) && nrow(final_res) > 0) {
    final_res$FDR          <- p.adjust(final_res$p_value, method = "fdr")
    final_res$Significance <- cut(final_res$FDR,
                                  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                  labels = c("***", "**", "*", "ns"))
    final_res <- final_res[order(final_res$p_value), , drop = FALSE]
    rownames(final_res) <- NULL
  } else {
    final_res <- data.frame(Cluster = character(), p_value = numeric(),
                            Method = character(), FDR = numeric(), Significance = character())
  }
  
  mean_prop <- stat_df |>
    dplyr::group_by(Cluster, .data[[group_col]]) |>
    dplyr::summarise(Mean_Prop = mean(Proportion, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = .data[[group_col]], values_from = Mean_Prop,
                       names_prefix = "Mean_Prop_")
  
  final_res <- merge(final_res, mean_prop, by = "Cluster", all.x = TRUE)
  
  write.csv(final_res, paste0("08_Differential_Abundance_Stats_", method, ".csv"), row.names = FALSE)
  message(">>> Statistical calculations completed and exported to CSV.")
  
  return(list(RawData = stat_df, Stats = final_res))
}

# ==============================================================================
#   New: Marker differential expression analysis
# ==============================================================================

#' Compare marker-expression differences between two groups within a specified cluster
#' @param cluster_name Cluster name to analyze (NULL means analyze all clusters)
#' @param group_col Grouping column name
#' @param markers List of markers to test (NULL means auto-detect all markers)
#' @param method Test method: "wilcoxon" or "t.test"
#' @return Data frame containing log2FC, p_value, and FDR
Perform_Marker_DE <- function(res_obj,
                              cluster_name  = NULL,
                              group_col     = "Group",
                              cluster_col   = "Cluster_Name",
                              markers       = NULL,
                              method        = c("wilcoxon", "t.test")) {
  method <- match.arg(method)
  df     <- res_obj$Data
  
  if (!group_col   %in% colnames(df)) stop(sprintf("Grouping column '%s' not found.", group_col))
  if (!cluster_col %in% colnames(df)) stop(sprintf("Cluster column '%s' not found.", cluster_col))
  
  groups <- unique(df[[group_col]])
  if (length(groups) != 2) stop("Perform_Marker_DE currently supports only two-group comparisons.")
  
  ignore_cols <- c("tSNE1", "tSNE2", "UMAP1", "UMAP2", "Cluster", "Cluster_Name",
                   "SampleID", "Group", "Time", "TIME", "CellType", "Pseudotime",
                   "Gate_Status", "id")
  if (is.null(markers)) {
    markers <- setdiff(colnames(df)[sapply(df, is.numeric)], ignore_cols)
  }
  
  clusters_to_run <- if (!is.null(cluster_name)) cluster_name else unique(df[[cluster_col]])
  
  all_results <- list()
  
  for (cl in clusters_to_run) {
    sub_df <- df[df[[cluster_col]] == cl, , drop = FALSE]
    if (nrow(sub_df) < 10) next
    
    grp_a <- sub_df[sub_df[[group_col]] == groups[1], ]
    grp_b <- sub_df[sub_df[[group_col]] == groups[2], ]
    if (nrow(grp_a) < 3 || nrow(grp_b) < 3) next
    
    res_rows <- lapply(markers, function(m) {
      if (!m %in% colnames(sub_df)) return(NULL)
      vals_a <- grp_a[[m]]
      vals_b <- grp_b[[m]]
      
      p_val <- tryCatch({
        if (method == "wilcoxon") {
          wilcox.test(vals_a, vals_b)$p.value
        } else {
          t.test(vals_a, vals_b)$p.value
        }
      }, error = function(e) NA_real_)
      
      mean_a  <- mean(vals_a, na.rm = TRUE)
      mean_b  <- mean(vals_b, na.rm = TRUE)
      log2fc  <- log2((mean_b + 1e-6) / (mean_a + 1e-6))
      
      data.frame(Cluster = cl, Marker = m,
                 Mean_GroupA = round(mean_a, 4),
                 Mean_GroupB = round(mean_b, 4),
                 Log2FC = round(log2fc, 4),
                 p_value = p_val,
                 stringsAsFactors = FALSE)
    })
    
    cl_df <- do.call(rbind, Filter(Negate(is.null), res_rows))
    if (!is.null(cl_df) && nrow(cl_df) > 0) {
      cl_df$FDR <- p.adjust(cl_df$p_value, method = "fdr")
      all_results[[cl]] <- cl_df
    }
  }
  
  if (length(all_results) == 0) {
    message(">>> No eligible clusters found, or the cell count is insufficient.")
    return(data.frame())
  }
  
  final_de <- do.call(rbind, all_results)
  colnames(final_de)[colnames(final_de) == "Mean_GroupA"] <- paste0("Mean_", groups[1])
  colnames(final_de)[colnames(final_de) == "Mean_GroupB"] <- paste0("Mean_", groups[2])
  final_de <- final_de[order(final_de$FDR), , drop = FALSE]
  rownames(final_de) <- NULL
  
  out_file <- paste0("DE_Markers_", method, ".csv")
  write.csv(final_de, out_file, row.names = FALSE)
  message(sprintf(">>> Marker DE analysis completed; results exported to: %s", out_file))
  
  return(final_de)
}

# ==============================================================================
# 2. Internal helper functions
# ==============================================================================

.safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

.pick_embedding_columns <- function(df, use_dim = "UMAP") {
  use_dim    <- toupper(use_dim)
  candidates <- if (use_dim == "UMAP") {
    list(c("UMAP1", "UMAP2"), c("UMAP_1", "UMAP_2"))
  } else {
    list(c("tSNE1", "tSNE2"), c("TSNE1", "TSNE2"), c("tSNE_1", "tSNE_2"))
  }
  for (cand in candidates) {
    if (all(cand %in% colnames(df))) return(cand)
  }
  stop(sprintf("Unable to find %s coordinate columns. Current columns include: %s", use_dim, paste(colnames(df), collapse = ", ")))
}

.pick_cluster_column <- function(df) {
  for (nm in c("Cluster_Name", "Cluster")) {
    if (nm %in% colnames(df)) return(nm)
  }
  stop("Cluster column not found; please provide at least 'Cluster_Name' or 'Cluster'.")
}

.summarize_scalar_info <- function(x, max_len = 120) {
  if (is.null(x) || length(x) == 0)    return(NA_character_)
  if (is.data.frame(x))                  return(sprintf("data.frame [%d x %d]", nrow(x), ncol(x)))
  if (is.matrix(x))                      return(sprintf("matrix [%d x %d]", nrow(x), ncol(x)))
  if (is.list(x) && !is.atomic(x))       return(sprintf("list [%d elements]", length(x)))
  if (length(x) > 8)                     return(sprintf("%s [%d values]", class(x)[1], length(x)))
  out <- paste(as.character(x), collapse = ", ")
  if (nchar(out) > max_len) out <- paste0(substr(out, 1, max_len), " ...")
  out
}

.save_png_plot <- function(plot_obj, file, width = 2200, height = 1600, res = 220) {
  grDevices::png(file, width = width, height = height, res = res)
  try(print(plot_obj), silent = TRUE)
  grDevices::dev.off()
  invisible(file)
}

.detect_marker_columns <- function(df) {
  meta_exact <- c("tSNE1", "tSNE2", "TSNE1", "TSNE2", "tSNE_1", "tSNE_2",
                  "UMAP1", "UMAP2", "UMAP_1", "UMAP_2",
                  "Cluster", "Cluster_Name", "CellType", "Group", "SampleID",
                  "Time", "TIME", "Pseudotime", "Gate_Status", "id")
  num_cols <- colnames(df)[vapply(df, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, meta_exact)
  num_cols[!grepl("^(Prob_|Probability_|RF_|Lineage_|Trajectory_)", num_cols)]
}

.best_label_col <- function(df) {
  for (nm in c("Cluster_Name", "CellType", "Cluster")) {
    if (nm %in% colnames(df)) return(nm)
  }
  NULL
}

# ==============================================================================
# 3. Trajectory inference (Slingshot) - pure ggplot2 plotting, replacing aes_string + graphics::lines
# ==============================================================================

#' Cell developmental trajectory inference (Slingshot)
#' @param start_by_col Optional; specify a cell-annotation column to define the start point, such as "CellType"
#' @param start_by_value Used together with start_by_col, e.g. "Naive CD4 T"
#' @param start_mode "majority" or "centroid"
Infer_Trajectory <- function(res_obj,
                             start_cluster  = NULL,
                             start_by_col   = NULL,
                             start_by_value = NULL,
                             start_mode     = c("majority", "centroid"),
                             use_dim        = "UMAP",
                             out_dir        = "09_Trajectory_Outputs",
                             prefix         = "Trajectory",
                             point_size     = 0.5,
                             seed           = 123) {
  
  if (!requireNamespace("slingshot",            quietly = TRUE)) stop("  Please install the 'slingshot' package first")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) stop("  Please install the 'SingleCellExperiment' package first")
  
  start_mode <- match.arg(start_mode)
  set.seed(seed)
  message("\n>>> Running Slingshot developmental trajectory inference...")
  
  df <- res_obj$Data
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) stop("res_obj$Data is empty.")
  
  embed_cols  <- .pick_embedding_columns(df, use_dim = use_dim)
  cluster_col <- .pick_cluster_column(df)
  coords      <- as.matrix(df[, embed_cols, drop = FALSE])
  rownames(coords) <- seq_len(nrow(df))
  
  if (any(!is.finite(coords))) stop(sprintf("%s coordinates contain NA/Inf.", toupper(use_dim)))
  
  cluster_labels <- as.character(df[[cluster_col]])
  if (all(is.na(cluster_labels))) stop("Cluster labels are all NA.")
  
  if (!is.null(start_cluster) && !start_cluster %in% unique(cluster_labels)) {
    warning(sprintf("start_cluster='%s' is not present in the cluster labels; automatic mapping will be attempted.", start_cluster))
    start_cluster <- NULL
  }
  
  resolved_start_cluster <- start_cluster
  resolved_start_n_cells <- NA_integer_
  
  if (is.null(resolved_start_cluster) && !is.null(start_by_col) && !is.null(start_by_value)) {
    if (!start_by_col %in% colnames(df)) {
      warning(sprintf("start_by_col='%s' does not exist; ignoring the custom start point.", start_by_col))
    } else {
      idx <- which(as.character(df[[start_by_col]]) %in% as.character(start_by_value))
      resolved_start_n_cells <- length(idx)
      if (length(idx) == 0) {
        warning(sprintf("Value '%s' was not found in column '%s'; Slingshot will infer the start point automatically.",
                        start_by_col, paste(start_by_value, collapse = ", ")))
      } else if (start_by_col == cluster_col) {
        resolved_start_cluster <- as.character(start_by_value[1])
      } else if (start_mode == "majority") {
        tb <- sort(table(cluster_labels[idx]), decreasing = TRUE)
        tb <- tb[names(tb) != ""]
        if (length(tb) > 0) resolved_start_cluster <- names(tb)[1]
      } else {
        target_center <- colMeans(coords[idx, , drop = FALSE], na.rm = TRUE)
        cl_levels     <- unique(cluster_labels)
        centers       <- sapply(cl_levels, function(cl) {
          colMeans(coords[cluster_labels == cl, , drop = FALSE], na.rm = TRUE)
        })
        if (is.null(dim(centers))) {
          centers <- matrix(centers, nrow = ncol(coords),
                            dimnames = list(colnames(coords), cl_levels))
        }
        dists <- apply(centers, 2, function(z) sqrt(sum((z - target_center)^2)))
        resolved_start_cluster <- names(which.min(dists))
      }
    }
  }
  
  if (!is.null(resolved_start_cluster) && !resolved_start_cluster %in% unique(cluster_labels)) {
    warning(sprintf("Resolved start cluster='%s' is not present in the cluster labels; Slingshot will infer it automatically.",
                    resolved_start_cluster))
    resolved_start_cluster <- NULL
  }
  
  dummy_expr <- matrix(0, nrow = 1, ncol = nrow(df),
                       dimnames = list("dummy_feature", seq_len(nrow(df))))
  col_df <- S4Vectors::DataFrame(cluster = factor(cluster_labels, levels = unique(cluster_labels)))
  if ("Group"    %in% colnames(df)) col_df$Group    <- df$Group
  if ("SampleID" %in% colnames(df)) col_df$SampleID <- df$SampleID
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = dummy_expr),
    colData = col_df
  )
  SingleCellExperiment::reducedDims(sce) <- S4Vectors::SimpleList(dimred = coords)
  
  sce <- slingshot::slingshot(sce, clusterLabels = "cluster",
                              reducedDim = "dimred",
                              start.clus = resolved_start_cluster)
  
  pst <- as.matrix(slingshot::slingPseudotime(sce))
  if (nrow(pst) != nrow(df)) stop("Slingshot returned pseudotime values inconsistent with the number of input cells.")
  
  df$Pseudotime <- suppressWarnings(rowMeans(pst, na.rm = TRUE))
  df$Pseudotime[!is.finite(df$Pseudotime)] <- NA_real_
  
  curve_weights <- slingshot::slingCurveWeights(sce)
  if (!is.null(curve_weights)) {
    curve_weights <- as.matrix(curve_weights)
    if (nrow(curve_weights) == nrow(df)) {
      df$Trajectory_Lineage <- apply(curve_weights, 1, function(x) {
        if (all(is.na(x)) || all(x <= 0)) return(NA_character_)
        paste0("Lineage_", which.max(x))
      })
    }
  }
  
  res_obj$Data       <- df
  res_obj$Trajectory <- list(
    method          = "slingshot",
    use_dim         = toupper(use_dim),
    embedding_cols  = embed_cols,
    cluster_col     = cluster_col,
    start_cluster   = resolved_start_cluster,
    start_specified = list(start_cluster  = start_cluster,
                           start_by_col   = start_by_col,
                           start_by_value = start_by_value,
                           start_mode     = start_mode),
    start_cells_n   = resolved_start_n_cells,
    n_lineages      = length(slingshot::slingCurves(sce))
  )
  
  .safe_dir_create(out_dir)
  base_name <- file.path(out_dir, prefix)
  plot_df   <- df
  
  # Extract Slingshot curve coordinates for pure ggplot2 geom_path (replacing graphics::lines)
  curve_df <- tryCatch({
    curves    <- slingshot::slingCurves(sce)
    curve_dfs <- lapply(seq_along(curves), function(i) {
      pts <- curves[[i]]$s[curves[[i]]$ord, ]
      data.frame(x = pts[, 1], y = pts[, 2], Lineage = paste0("Lineage_", i))
    })
    do.call(rbind, curve_dfs)
  }, error = function(e) NULL)
  
  subtitle_start <- ifelse(is.null(resolved_start_cluster), "auto", resolved_start_cluster)
  
  # Plot 1: colored by pseudotime
  png1 <- paste0(base_name, "_Pseudotime.png")
  g1 <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = .data[[embed_cols[1]]],
                                     y = .data[[embed_cols[2]]],
                                     color = Pseudotime)) +
    ggplot2::geom_point(size = point_size, alpha = 0.75, na.rm = TRUE) +
    ggplot2::scale_color_viridis_c(option = "plasma", na.value = "grey80") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(title    = sprintf("Slingshot Trajectory on %s", toupper(use_dim)),
                  subtitle = sprintf("Colored by pseudotime | start: %s", subtitle_start),
                  x = embed_cols[1], y = embed_cols[2])
  if (!is.null(curve_df)) {
    g1 <- g1 + ggplot2::geom_path(
      data = curve_df,
      ggplot2::aes(x = x, y = y, group = Lineage),
      color = "black", linewidth = 1.2, inherit.aes = FALSE
    )
  }
  .save_png_plot(g1, png1, width = 2200, height = 1800)
  
  # Plot 2: colored by cluster
  png2 <- paste0(base_name, "_Cluster.png")
  g2 <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = .data[[embed_cols[1]]],
                                     y = .data[[embed_cols[2]]],
                                     color = .data[[cluster_col]])) +
    ggplot2::geom_point(size = point_size, alpha = 0.75, na.rm = TRUE) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(title    = sprintf("Slingshot Trajectory on %s", toupper(use_dim)),
                  subtitle = sprintf("Colored by %s", cluster_col),
                  x = embed_cols[1], y = embed_cols[2], color = cluster_col)
  if (!is.null(curve_df)) {
    g2 <- g2 + ggplot2::geom_path(
      data = curve_df,
      ggplot2::aes(x = x, y = y, group = Lineage),
      color = "black", linewidth = 1.2, inherit.aes = FALSE
    )
  }
  .save_png_plot(g2, png2, width = 2200, height = 1800)
  
  # Plot 3: histogram of pseudotime distribution
  png3 <- paste0(base_name, "_Pseudotime_Distribution.png")
  g3 <- ggplot2::ggplot(
    plot_df[is.finite(plot_df$Pseudotime), , drop = FALSE],
    ggplot2::aes(x = Pseudotime)
  ) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.8, fill = "#3498DB", color = "white") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(title = "Pseudotime Distribution", x = "Pseudotime", y = "Cell count")
  .save_png_plot(g3, png3, width = 1800, height = 1400)
  
  res_obj$Trajectory$files <- list(
    pseudotime_png   = png1,
    cluster_png      = png2,
    distribution_png = png3
  )
  
  message(sprintf(">>> Trajectory inference completed! Start cluster: %s.", subtitle_start))
  return(list(res_obj = res_obj, sce = sce, plot_files = res_obj$Trajectory$files))
}


# ==============================================================================
# 4. FCS export
# ==============================================================================

#' Export clean FCS files (importable into FlowJo)
Export_To_FCS <- function(res_obj, out_dir = "Annotated_FCS_Outputs", sample_col = "SampleID") {
  if (!requireNamespace("flowCore", quietly = TRUE)) stop("  Please install the 'flowCore' package first")
  if (!requireNamespace("Biobase",  quietly = TRUE)) stop("  Please install the 'Biobase' package first")
  
  .safe_dir_create(out_dir)
  df <- res_obj$Data
  if (!sample_col %in% colnames(df)) stop(sprintf("Sample column '%s' not found.", sample_col))
  samples <- unique(df[[sample_col]])
  
  message("\n>>> Exporting data to FlowJo-compatible FCS format...")
  for (samp in samples) {
    sub_df <- df[df[[sample_col]] == samp, , drop = FALSE]
    num_df <- sub_df[, sapply(sub_df, is.numeric), drop = FALSE]
    
    if (ncol(num_df) == 0) {
      warning(sprintf("Sample %s has no numeric columns and was skipped.", samp))
      next
    }
    
    if ("Cluster"      %in% colnames(sub_df)) num_df$Cluster_ID  <- as.numeric(as.factor(sub_df$Cluster))
    if ("Cluster_Name" %in% colnames(sub_df)) num_df$CellType_ID <- as.numeric(as.factor(sub_df$Cluster_Name))
    
    meta          <- data.frame(name = colnames(num_df), desc = colnames(num_df), stringsAsFactors = FALSE)
    meta$range    <- apply(num_df, 2, function(x) max(x, na.rm = TRUE))
    meta$minRange <- apply(num_df, 2, function(x) min(x, na.rm = TRUE))
    meta$maxRange <- meta$range
    
    ff <- flowCore::flowFrame(as.matrix(num_df), parameters = Biobase::AnnotatedDataFrame(meta))
    flowCore::write.FCS(ff, file.path(out_dir, paste0("Annotated_", samp, ".fcs")))
  }
  message(">>> All FCS files were successfully exported to: ", out_dir)
}

# ==============================================================================
# 5. Final HTML report
# ==============================================================================

#' Final streamlined HTML report
#' @param template "manuscript" manuscript supplement; "review" internal review; "full" full version
Generate_Final_HTML_Report <- function(
    res_obj,
    stats_res    = NULL,
    project_name = "EasyFlow2 Final Report",
    out_file     = "EasyFlow2_Final_Report.html",
    sample_col   = "SampleID",
    group_col    = "Group",
    template     = c("manuscript", "review", "full"),
    include_plots = NULL,
    top_markers_n = 10,
    top_stats_n   = 15,
    true_label    = "CellType",
    pred_label    = "Cluster_Name") {
  
  template <- match.arg(template)
  df       <- res_obj$Data
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) stop("res_obj$Data is empty.")
  if (!sample_col %in% colnames(df)) stop(sprintf("Missing sample column '%s'.", sample_col))
  
  has_group  <- group_col %in% colnames(df)
  label_col  <- .best_label_col(df)
  embed_cols <- tryCatch(.pick_embedding_columns(df, "UMAP"), error = function(e) NULL)
  embed_label <- if (!is.null(embed_cols)) "UMAP" else NULL
  if (is.null(embed_cols)) {
    embed_cols  <- tryCatch(.pick_embedding_columns(df, "TSNE"), error = function(e) NULL)
    if (!is.null(embed_cols)) embed_label <- "tSNE"
  }
  marker_cols <- .detect_marker_columns(df)
  
  default_plots <- switch(
    template,
    manuscript = c("sample_counts", "embedding", "composition", "marker_dotplot", "trajectory", "abundance"),
    review     = c("sample_counts", "embedding", "composition", "marker_dotplot", "sample_pca", "trajectory", "abundance", "confusion"),
    full       = c("sample_counts", "embedding", "composition", "marker_dotplot", "sample_pca", "trajectory", "abundance", "confusion", "metadata")
  )
  if (is.null(include_plots)) include_plots <- default_plots
  include_plots <- unique(include_plots)
  
  asset_dir  <- .safe_dir_create(file.path(getwd(), "09_Final_Report_Assets"))
  plot_files <- list()
  
  # 1. Sample counts
  counts_df <- as.data.frame(table(df[[sample_col]],
                                   if (has_group) df[[group_col]] else factor("All"),
                                   useNA = "no"))
  colnames(counts_df) <- c("Sample", "Group", "Cell_Count")
  counts_df <- counts_df[counts_df$Cell_Count > 0, , drop = FALSE]
  if ("sample_counts" %in% include_plots && nrow(counts_df) > 0) {
    f <- file.path(asset_dir, "01_sample_counts.png")
    p <- ggplot2::ggplot(counts_df, ggplot2::aes(x = Sample, y = Cell_Count, fill = Group)) +
      ggplot2::geom_col() +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = "Cells retained per sample", x = "Sample", y = "Cell count")
    .save_png_plot(p, f, height = 1400)
    plot_files$sample_counts <- normalizePath(f, winslash = "/", mustWork = FALSE)
  }
  
  # 2. Embedding
  if ("embedding" %in% include_plots && !is.null(embed_cols) && !is.null(label_col)) {
    f <- file.path(asset_dir, "02_embedding.png")
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[embed_cols[1]]],
                                          y = .data[[embed_cols[2]]],
                                          color = .data[[label_col]])) +
      ggplot2::geom_point(size = 0.25, alpha = 0.7, na.rm = TRUE) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::labs(title    = sprintf("%s overview", embed_label),
                    subtitle = sprintf("Colored by %s", label_col),
                    x = embed_cols[1], y = embed_cols[2], color = label_col)
    .save_png_plot(p, f, height = 1800)
    plot_files$embedding <- normalizePath(f, winslash = "/", mustWork = FALSE)
  }
  
  # 3. Composition
  composition_df <- data.frame()
  if ("composition" %in% include_plots && has_group && !is.null(label_col)) {
    sample_comp <- df |>
      dplyr::count(.data[[sample_col]], .data[[group_col]], .data[[label_col]], name = "n_cells") |>
      dplyr::group_by(.data[[sample_col]]) |>
      dplyr::mutate(prop = n_cells / sum(n_cells)) |>
      dplyr::ungroup()
    
    composition_df <- sample_comp |>
      dplyr::group_by(.data[[group_col]], .data[[label_col]]) |>
      dplyr::summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop")
    
    f <- file.path(asset_dir, "03_composition.png")
    p <- ggplot2::ggplot(composition_df,
                         ggplot2::aes(x = .data[[group_col]], y = mean_prop,
                                      fill = .data[[label_col]])) +
      ggplot2::geom_col(position = "fill") +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::labs(title = "Mean sample-level composition by group",
                    x = group_col, y = "Relative abundance", fill = label_col)
    .save_png_plot(p, f, height = 1600)
    plot_files$composition <- normalizePath(f, winslash = "/", mustWork = FALSE)
  }
  
  # 4. Marker dotplot
  selected_markers <- character()
  if ("marker_dotplot" %in% include_plots && !is.null(label_col) && length(marker_cols) > 1) {
    marker_var       <- vapply(marker_cols, function(m) stats::var(df[[m]], na.rm = TRUE), numeric(1))
    selected_markers <- names(utils::head(sort(marker_var, decreasing = TRUE), top_markers_n))
    
    scaled_df <- df
    for (m in selected_markers) scaled_df[[m]] <- as.numeric(scale(scaled_df[[m]]))
    
    dotplot_df <- scaled_df |>
      dplyr::select(dplyr::all_of(c(label_col, selected_markers))) |>
      tidyr::pivot_longer(cols = dplyr::all_of(selected_markers),
                          names_to = "Marker", values_to = "Expression") |>
      dplyr::group_by(.data[[label_col]], Marker) |>
      dplyr::summarise(Mean_Exp = mean(Expression, na.rm = TRUE),
                       Pct_Exp  = mean(Expression > 0, na.rm = TRUE) * 100,
                       .groups  = "drop")
    
    f <- file.path(asset_dir, "04_marker_dotplot.png")
    p <- ggplot2::ggplot(dotplot_df,
                         ggplot2::aes(x = Marker, y = .data[[label_col]])) +
      ggplot2::geom_point(ggplot2::aes(size = Pct_Exp, color = Mean_Exp)) +
      ggplot2::scale_size_continuous(range = c(0, 6), limits = c(0, 100), name = "% expressing") +
      ggplot2::scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3",
                                     midpoint = 0, name = "Mean Z-score") +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = "Top variable marker dotplot", x = "Marker", y = label_col)
    .save_png_plot(p, f, height = 1700)
    plot_files$marker_dotplot <- normalizePath(f, winslash = "/", mustWork = FALSE)
  }
  
  # 5. Sample PCA
  if ("sample_pca" %in% include_plots && has_group && length(marker_cols) >= 2) {
    sample_mat <- df |>
      dplyr::group_by(.data[[sample_col]], .data[[group_col]]) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(marker_cols),
                                     ~stats::median(.x, na.rm = TRUE)), .groups = "drop")
    if (nrow(sample_mat) >= 3) {
      pca <- stats::prcomp(sample_mat[, marker_cols, drop = FALSE], center = TRUE, scale. = TRUE)
      sample_pca_df <- data.frame(
        Sample = sample_mat[[sample_col]],
        Group  = sample_mat[[group_col]],
        PC1    = pca$x[, 1],
        PC2    = pca$x[, 2]
      )
      f <- file.path(asset_dir, "05_sample_pca.png")
      p <- ggplot2::ggplot(sample_pca_df,
                           ggplot2::aes(x = PC1, y = PC2, color = Group, label = Sample)) +
        ggplot2::geom_point(size = 4, alpha = 0.9) +
        ggplot2::geom_text(nudge_y = 0.05, size = 3, show.legend = FALSE) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::labs(title = "Sample-level PCA", x = "PC1", y = "PC2")
      .save_png_plot(p, f, height = 1600)
      plot_files$sample_pca <- normalizePath(f, winslash = "/", mustWork = FALSE)
    }
  }
  
  # 6. Trajectory
  trajectory_summary <- data.frame()
  if ("trajectory" %in% include_plots && !is.null(res_obj$Trajectory)) {
    trajectory_summary <- data.frame(
      Item  = c("Method", "Embedding", "Cluster column", "Start cluster", "Lineages"),
      Value = c(.summarize_scalar_info(res_obj$Trajectory$method),
                .summarize_scalar_info(res_obj$Trajectory$use_dim),
                .summarize_scalar_info(res_obj$Trajectory$cluster_col),
                .summarize_scalar_info(res_obj$Trajectory$start_cluster),
                .summarize_scalar_info(res_obj$Trajectory$n_lineages)),
      stringsAsFactors = FALSE
    )
    if (is.list(res_obj$Trajectory$files)) {
      pst_f <- res_obj$Trajectory$files$pseudotime_png
      dst_f <- res_obj$Trajectory$files$distribution_png
      if (!is.null(pst_f) && file.exists(pst_f))
        plot_files$trajectory <- normalizePath(pst_f, winslash = "/", mustWork = FALSE)
      if (!is.null(dst_f) && file.exists(dst_f))
        plot_files$trajectory_distribution <- normalizePath(dst_f, winslash = "/", mustWork = FALSE)
    }
  }
  
  # 7. Differential abundance
  stats_table      <- data.frame()
  top_stats        <- data.frame()
  abundance_bar_df <- data.frame()
  if ("abundance" %in% include_plots && !is.null(stats_res) &&
      is.data.frame(stats_res$Stats) && nrow(stats_res$Stats) > 0) {
    stats_table <- stats_res$Stats
    ord_col     <- if ("FDR" %in% colnames(stats_table)) "FDR" else "p_value"
    stats_table <- stats_table[order(stats_table[[ord_col]], na.last = TRUE), , drop = FALSE]
    top_stats   <- utils::head(stats_table, top_stats_n)
    
    mean_cols <- grep("^Mean_Prop_", colnames(top_stats), value = TRUE)
    if (length(mean_cols) == 2) {
      abundance_bar_df        <- top_stats
      abundance_bar_df$Effect <- abundance_bar_df[[mean_cols[2]]] - abundance_bar_df[[mean_cols[1]]]
      abundance_bar_df$Cluster <- factor(abundance_bar_df$Cluster,
                                         levels = rev(abundance_bar_df$Cluster))
      f <- file.path(asset_dir, "06_abundance_effect.png")
      p <- ggplot2::ggplot(abundance_bar_df, ggplot2::aes(x = Effect, y = Cluster)) +
        ggplot2::geom_col() +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::labs(title = sprintf("Top differential abundance (%s - %s)",
                                      mean_cols[2], mean_cols[1]),
                      x = "Mean proportion difference", y = "Cluster")
      .save_png_plot(p, f, width = 1800, height = 1500)
      plot_files$abundance_effect <- normalizePath(f, winslash = "/", mustWork = FALSE)
    }
  }
  
  # 8. Confusion matrix
  if ("confusion" %in% include_plots && template != "manuscript" &&
      true_label %in% colnames(df) && pred_label %in% colnames(df) &&
      true_label != pred_label) {
    tab <- table(Truth = df[[true_label]], Prediction = df[[pred_label]])
    if (all(dim(tab) > 1)) {
      confusion_df         <- as.data.frame(as.table(prop.table(tab, margin = 1) * 100))
      colnames(confusion_df) <- c("Truth", "Prediction", "Percent")
      f <- file.path(asset_dir, "07_confusion_matrix.png")
      p <- ggplot2::ggplot(confusion_df,
                           ggplot2::aes(x = Prediction, y = Truth, fill = Percent)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", Percent)), size = 3) +
        ggplot2::scale_fill_gradient(low = "white", high = "#2C3E50") +
        ggplot2::theme_classic(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "Confusion matrix", x = pred_label, y = true_label)
      .save_png_plot(p, f, height = 1800)
      plot_files$confusion <- normalizePath(f, winslash = "/", mustWork = FALSE)
    }
  }
  
  # 9. Metadata
  metadata_df <- data.frame(Parameter = character(), Value = character(), stringsAsFactors = FALSE)
  if ("metadata" %in% include_plots && template == "full") {
    if (!is.null(res_obj$Info)) {
      metadata_df <- data.frame(
        Parameter = names(res_obj$Info),
        Value     = vapply(res_obj$Info, .summarize_scalar_info, character(1)),
        stringsAsFactors = FALSE
      )
    }
    if (!is.null(res_obj$Trajectory)) {
      tmp <- data.frame(
        Parameter = paste0("Trajectory.", names(res_obj$Trajectory)),
        Value     = vapply(res_obj$Trajectory, .summarize_scalar_info, character(1)),
        stringsAsFactors = FALSE
      )
      metadata_df <- rbind(metadata_df, tmp)
    }
    metadata_df <- metadata_df[!is.na(metadata_df$Value), , drop = FALSE]
  }
  
  summary_lines <- c(
    sprintf("Cells analysed: %d", nrow(df)),
    sprintf("Samples: %d", length(unique(df[[sample_col]]))),
    if (has_group) sprintf("Groups: %d", length(unique(df[[group_col]]))) else NULL,
    if (!is.null(label_col)) sprintf("Primary label: %s (%d levels)", label_col,
                                     length(unique(df[[label_col]]))) else NULL,
    if (!is.null(embed_label)) sprintf("Embedding: %s", embed_label) else NULL,
    if ("Pseudotime" %in% colnames(df)) sprintf("Cells with pseudotime: %d",
                                                sum(is.finite(df$Pseudotime))) else NULL
  )
  summary_lines <- summary_lines[!is.na(summary_lines)]
  
  plot_files <- Filter(function(x) is.character(x) && length(x) == 1 &&
                         !is.na(x) && file.exists(x), plot_files)
  
  rmd <- c(
    "---",
    paste0('title: "', project_name, '"'),
    'date: "`r Sys.Date()`"',
    "output:",
    "  html_document:",
    "    theme: flatly",
    "    toc: true",
    "    toc_float: true",
    "    df_print: paged",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(DT)",
    "```",
    "",
    "## Executive Summary",
    "```{r}",
    "cat(paste0('- ', summary_lines, collapse='\\n'))",
    "```",
    ""
  )
  
  if (!is.null(plot_files$sample_counts) && file.exists(plot_files$sample_counts))
    rmd <- c(rmd, "## 1. Sample Overview",
             "```{r, out.width='100%'}", "knitr::include_graphics(plot_files$sample_counts)", "```", "")
  
  if (!is.null(plot_files$embedding) && file.exists(plot_files$embedding))
    rmd <- c(rmd, "## 2. Cell Atlas Overview",
             "```{r, out.width='100%'}", "knitr::include_graphics(plot_files$embedding)", "```", "")
  
  if (!is.null(plot_files$composition) && file.exists(plot_files$composition))
    rmd <- c(rmd, "## 3. Group Composition",
             "```{r, out.width='100%'}", "knitr::include_graphics(plot_files$composition)", "```", "")
  
  if (!is.null(plot_files$marker_dotplot) && file.exists(plot_files$marker_dotplot))
    rmd <- c(rmd, "## 4. Marker Identity Overview",
             "```{r, out.width='100%'}", "knitr::include_graphics(plot_files$marker_dotplot)", "```", "")
  
  if (!is.null(plot_files$sample_pca) && file.exists(plot_files$sample_pca))
    rmd <- c(rmd, "## 5. Sample-Level Structure",
             "```{r, out.width='100%'}", "knitr::include_graphics(plot_files$sample_pca)", "```", "")
  
  if (!is.null(plot_files$trajectory) && file.exists(plot_files$trajectory)) {
    rmd <- c(rmd, "## 6. Trajectory / Pseudotime",
             "```{r, out.width='100%'}", "knitr::include_graphics(plot_files$trajectory)", "```", "")
    if (!is.null(plot_files$trajectory_distribution) &&
        file.exists(plot_files$trajectory_distribution))
      rmd <- c(rmd, "```{r, out.width='100%'}",
               "knitr::include_graphics(plot_files$trajectory_distribution)", "```", "")
    if (nrow(trajectory_summary) > 0)
      rmd <- c(rmd, "```{r}",
               "datatable(trajectory_summary, options = list(dom='t', paging = FALSE))",
               "```", "")
  }
  
  if (nrow(top_stats) > 0) {
    rmd <- c(rmd, "## 7. Differential Abundance",
             "```{r}",
             "dt <- datatable(top_stats, options = list(pageLength = min(15, nrow(top_stats)), dom = 'tip'))",
             "cols_to_round <- intersect(c('p_value', 'FDR'), colnames(top_stats))",
             "if (length(cols_to_round) > 0) dt <- formatRound(dt, columns = cols_to_round, digits = 4)",
             "dt", "```", "")
    if (!is.null(plot_files$abundance_effect) && file.exists(plot_files$abundance_effect))
      rmd <- c(rmd, "```{r, out.width='90%'}",
               "knitr::include_graphics(plot_files$abundance_effect)", "```", "")
  }
  
  if (!is.null(plot_files$confusion) && file.exists(plot_files$confusion))
    rmd <- c(rmd, "## 8. Annotation / Projection QC",
             "```{r, out.width='100%'}", "knitr::include_graphics(plot_files$confusion)", "```", "")
  
  if (nrow(metadata_df) > 0)
    rmd <- c(rmd, "## 9. Pipeline Metadata",
             "```{r}",
             "datatable(metadata_df, options = list(pageLength = 12, dom = 'tip'))",
             "```", "")
  
  rmd <- c(rmd, "***", "*Final report generated automatically by EasyFlow2.*")
  
  rmd_file   <- tempfile(fileext = ".Rmd")
  writeLines(rmd, rmd_file)
  
  render_env <- new.env(parent = globalenv())
  assign("plot_files",         plot_files,         envir = render_env)
  assign("summary_lines",      summary_lines,      envir = render_env)
  assign("trajectory_summary", trajectory_summary, envir = render_env)
  assign("top_stats",          top_stats,          envir = render_env)
  assign("metadata_df",        metadata_df,        envir = render_env)
  
  rmarkdown::render(rmd_file, output_file = out_file, output_dir = getwd(),
                    envir = render_env, quiet = TRUE, knit_root_dir = getwd())
  message(">>> Final streamlined HTML report generated: ", file.path(getwd(), out_file))
  
  invisible(list(
    out_file         = file.path(getwd(), out_file),
    plot_files       = plot_files,
    selected_markers = selected_markers,
    template         = template,
    included         = include_plots
  ))
}

# ==============================================================================
# Module 6: Industrial-/publication-grade plotting and visualization toolbox (Publication-Ready)
# Filename: R/06_Visualization.R
# ==============================================================================

#' Unified publication-ready plotting theme
theme_publication <- function(base_size = 14, base_family = "sans") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      panel.background   = ggplot2::element_blank(),
      #   Change 1: size -> linewidth to eliminate deprecation warnings in ggplot2 >= 3.4
      panel.border       = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.2),
      panel.grid.major   = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      strip.background   = ggplot2::element_rect(fill = "#EFEFEF", color = "black", linewidth = 1.2),
      strip.text         = ggplot2::element_text(face = "bold", size = base_size),
      axis.text          = ggplot2::element_text(color = "black", size = base_size * 0.9),
      axis.title         = ggplot2::element_text(face = "bold", size = base_size),
      legend.background  = ggplot2::element_blank(),
      legend.key         = ggplot2::element_blank(),
      legend.text        = ggplot2::element_text(size = base_size * 0.9),
      plot.title         = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size * 1.2)
    )
}



# ==============================================================================
# 0. Global plotting order control
# ==============================================================================

#' Define global plotting order for all downstream figures
#'
#' @param res_obj Standard result object
#' @param order_list Named list. Example:
#'   list(Group = c("Control", "Model", "Treatment"),
#'        Cluster_Name = c("Naive T", "Effector T", "Treg"))
#'
#' @return Updated result object with factor levels stored in res_obj$Data and
#'   the original order definition stored in res_obj$Plot_Order
#' @export
Set_Plot_Order <- function(res_obj, order_list = list()) {
  if (length(order_list) == 0) return(res_obj)
  if (is.null(res_obj$Data) || !is.data.frame(res_obj$Data)) {
    stop("  res_obj$Data is missing or is not a data.frame")
  }

  df <- res_obj$Data
  clean_order_list <- list()

  for (nm in names(order_list)) {
    if (!(nm %in% colnames(df))) {
      warning(paste0("   Column not found, skipped ordering: ", nm))
      next
    }

    user_levels <- unique(as.character(order_list[[nm]]))
    existing_levels <- unique(as.character(df[[nm]]))
    final_levels <- unique(c(user_levels, setdiff(existing_levels, user_levels)))

    df[[nm]] <- factor(as.character(df[[nm]]), levels = final_levels)
    clean_order_list[[nm]] <- final_levels
  }

  res_obj$Data <- df
  if (is.null(res_obj$Plot_Order)) res_obj$Plot_Order <- list()
  res_obj$Plot_Order <- utils::modifyList(res_obj$Plot_Order, clean_order_list)
  return(res_obj)
}

#' Internal helper: apply stored plotting order to a data.frame
#' @keywords internal
Apply_Plot_Order <- function(df, res_obj = NULL, cols = NULL, order_list = NULL) {
  merged_order <- list()
  if (!is.null(res_obj$Plot_Order)) merged_order <- res_obj$Plot_Order
  if (!is.null(order_list)) merged_order <- utils::modifyList(merged_order, order_list)
  if (length(merged_order) == 0) return(df)

  target_cols <- names(merged_order)
  if (!is.null(cols)) target_cols <- intersect(target_cols, cols)
  target_cols <- intersect(target_cols, colnames(df))

  for (nm in target_cols) {
    user_levels <- unique(as.character(merged_order[[nm]]))
    existing_levels <- unique(as.character(df[[nm]]))
    final_levels <- unique(c(user_levels, setdiff(existing_levels, user_levels)))
    df[[nm]] <- factor(as.character(df[[nm]]), levels = final_levels)
  }

  return(df)
}

# ==============================================================================
# 1. Dimensionality-reduction visualization (UMAP / tSNE)
# ==============================================================================

#' tSNE density plot
#' @param res_obj Standard result object
#' @param group_by Column used for grouping/coloring; if NULL, show overall density
#' @param bins Number of hexagonal bins
#' @param alpha Transparency
Plot_tSNE_Density <- function(res_obj, group_by = NULL, bins = 60, alpha = 0.9) {
  #   Change 2: implement directly without relying on the missing Plot_2D_Density
  df <- Apply_Plot_Order(res_obj$Data, res_obj = res_obj, cols = group_by)
  if (!all(c("tSNE1", "tSNE2") %in% colnames(df))) stop("  tSNE1/tSNE2 coordinate columns are missing from the data")
  
  if (!is.null(group_by) && group_by %in% colnames(df)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = tSNE1, y = tSNE2)) +
      ggplot2::geom_hex(bins = bins, alpha = alpha) +
      ggplot2::scale_fill_viridis_c(option = "plasma") +
      ggplot2::facet_wrap(stats::as.formula(paste("~", group_by))) +
      theme_publication() +
      ggplot2::ggtitle(paste0("tSNE Density by ", group_by))
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = tSNE1, y = tSNE2)) +
      ggplot2::geom_hex(bins = bins, alpha = alpha) +
      ggplot2::scale_fill_viridis_c(option = "plasma") +
      theme_publication() +
      ggplot2::ggtitle("tSNE Density")
  }
  return(p)
}

#' Dimensionality-reduction scatter plot (UMAP / tSNE) with intelligent handling of anti-overplotting and layer order for million-cell data
Plot_Dim_Reduction <- function(res_obj, dim_method = c("UMAP", "tSNE"), color_by = "Cluster", split_by = NULL) {
  df         <- Apply_Plot_Order(res_obj$Data, res_obj = res_obj, cols = c(color_by, split_by))
  dim_method <- match.arg(dim_method)
  x_col      <- paste0(dim_method, "1")
  y_col      <- paste0(dim_method, "2")
  
  if (!(color_by %in% colnames(df))) stop(paste("  Coloring variable not found in the data:", color_by))
  
  # Intelligent layer sorting: send gray background cells to the bottom to avoid obscuring colored cells
  if ("Ungated" %in% unique(df[[color_by]])) {
    df <- df[order(df[[color_by]] != "Ungated"), ]
  } else {
    df <- df[sample(nrow(df)), ]
  }
  
  # Dynamically shrink point size according to cell count
  pt_size  <- ifelse(nrow(df) > 500000, 0.05, ifelse(nrow(df) > 100000, 0.2, 0.5))
  pt_alpha <- ifelse(nrow(df) > 500000, 0.3, 0.8)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_by]])) +
    ggplot2::geom_point(shape = 16, size = pt_size, alpha = pt_alpha) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  if (is.numeric(df[[color_by]])) {
    p <- p + ggplot2::scale_color_viridis_c(option = "plasma")
  } else {
    n_colors <- length(unique(df[[color_by]]))
    if ("Ungated" %in% unique(df[[color_by]])) {
      my_pal <- setNames(
        c("#E5E5E5", colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(n_colors - 1)),
        c("Ungated", setdiff(unique(df[[color_by]]), "Ungated"))
      )
      p <- p + ggplot2::scale_color_manual(values = my_pal)
    } else {
      p <- p + ggplot2::scale_color_manual(
        values = colorRampPalette(RColorBrewer::brewer.pal(min(12, n_colors), "Set3"))(n_colors)
      )
    }
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1)))
  }
  
  if (!is.null(split_by)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", split_by)))
  
  return(p + ggplot2::ggtitle(paste(dim_method, "Map colored by", color_by)))
}

# ==============================================================================
# 2. Population composition and abundance (Composition & Alluvial)
# ==============================================================================

#' Stacked bar plot of cell composition proportions
#' @param position "fill" = 100% stacked; "stack" = stack by absolute counts
#' @param correct_sample_bias Whether to calculate proportions by sample before taking group means (recommended TRUE to avoid sample-size bias)
Plot_Composition <- function(res_obj, x_var = "Group", fill_var = "Cluster_Name",
                             position = c("fill", "stack"),
                             correct_sample_bias = TRUE,
                             sample_col = "SampleID") {
  df       <- Apply_Plot_Order(res_obj$Data, res_obj = res_obj, cols = c(x_var, fill_var, sample_col))
  position <- match.arg(position)
  
  # stack stack mode always uses raw counts; fill mode can optionally apply sample-bias correction
  if (position == "stack" || !correct_sample_bias || !sample_col %in% colnames(df)) {
    count_df <- df |>
      dplyr::group_by(.data[[x_var]], .data[[fill_var]]) |>
      dplyr::summarise(Count = dplyr::n(), .groups = "drop")
    y_label <- if (position == "fill") "Proportion (100%)" else "Cell Count"
  } else {
    # fill + correct_sample_bias: first calculate proportions by sample and then take group means
    count_df <- df |>
      dplyr::count(.data[[sample_col]], .data[[x_var]], .data[[fill_var]], name = "n_cells") |>
      dplyr::group_by(.data[[sample_col]]) |>
      dplyr::mutate(prop = n_cells / sum(n_cells)) |>
      dplyr::ungroup() |>
      dplyr::group_by(.data[[x_var]], .data[[fill_var]]) |>
      dplyr::summarise(Count = mean(prop, na.rm = TRUE), .groups = "drop")
    y_label <- "Mean Sample Proportion"
  }
  
  count_df <- Apply_Plot_Order(count_df, res_obj = res_obj, cols = c(x_var, fill_var))

  p <- ggplot2::ggplot(count_df, ggplot2::aes(x = .data[[x_var]], y = Count, fill = .data[[fill_var]])) +
    ggplot2::geom_bar(stat = "identity", position = position, color = "black", linewidth = 0.3) +
    theme_publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = x_var, y = y_label, fill = fill_var) +
    ggplot2::ggtitle(paste("Cell Composition by", x_var))
  
  n_colors <- length(unique(count_df[[fill_var]]))
  p <- p + ggplot2::scale_fill_manual(
    values = colorRampPalette(RColorBrewer::brewer.pal(min(12, n_colors), "Set3"))(n_colors)
  )
  message("   [Figure 2] Cell composition plot generated")
  return(p)
}

#' Sankey/Alluvial plot (Alluvial Plot)
Plot_Alluvial <- function(res_obj, axis1 = "Group", axis2 = "Cluster_Name", fill_var = "Cluster_Name") {
  df <- Apply_Plot_Order(res_obj$Data, res_obj = res_obj, cols = c(axis1, axis2, fill_var))
  
  group_cols <- unique(c(axis1, axis2, fill_var))
  count_df <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(Freq = dplyr::n(), .groups = "drop")
  
  count_df <- Apply_Plot_Order(count_df, res_obj = res_obj, cols = c(axis1, axis2, fill_var))

  p <- ggplot2::ggplot(count_df, ggplot2::aes(y = Freq, axis1 = .data[[axis1]], axis2 = .data[[axis2]])) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = .data[[fill_var]]), width = 1/12, alpha = 0.7) +
    ggalluvial::geom_stratum(width = 1/8, color = "black", fill = "grey90") +
    ggplot2::geom_text(
      mapping  = ggplot2::aes(label = ggplot2::after_stat(stratum)),
      stat     = ggalluvial::StatStratum,
      size     = 3.5,
      fontface = "bold"
    ) +
    ggplot2::scale_x_discrete(limits = c(axis1, axis2), expand = c(.05, .05)) +
    theme_publication() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank()) +
    ggplot2::ggtitle("Alluvial / Sankey Flow Plot")
  
  n_colors <- length(unique(count_df[[fill_var]]))
  p <- p + ggplot2::scale_fill_manual(
    values = colorRampPalette(RColorBrewer::brewer.pal(min(12, n_colors), "Set3"))(n_colors)
  )
  message("   [Figure 3] Sankey/Alluvial plot generated")
  return(p)
}

# ==============================================================================
# 3. Expression-profile visualization (Dotplot & Heatmap & Violin)
# ==============================================================================

#' Standard single-cell dot plot (Dotplot: color = expression level, size = positive proportion)
Plot_Dotplot <- function(res_obj, markers = NULL, group_by = "Cluster_Name") {
  df <- Apply_Plot_Order(res_obj$Data, res_obj = res_obj, cols = group_by)
  
  meta_cols <- c("tSNE1", "tSNE2", "UMAP1", "UMAP2", "Cluster", "Group", "SampleID",
                 "Time", "TIME", "CellType", "id", "Cluster_Name", "Gate_Status", "Pseudotime")
  if (is.null(markers)) markers <- setdiff(colnames(df)[sapply(df, is.numeric)], meta_cols)
  
  scaled_df <- df
  for (m in markers) scaled_df[[m]] <- as.numeric(scale(scaled_df[[m]]))
  
  plot_data <- scaled_df |>
    dplyr::select(dplyr::all_of(c(group_by, markers))) |>
    tidyr::pivot_longer(cols = dplyr::all_of(markers), names_to = "Marker", values_to = "Expression") |>
    dplyr::group_by(.data[[group_by]], Marker) |>
    dplyr::summarise(
      Mean_Exp = mean(Expression, na.rm = TRUE),
      Pct_Exp  = sum(Expression > 0, na.rm = TRUE) / dplyr::n() * 100,
      .groups  = "drop"
    )

  plot_data <- Apply_Plot_Order(plot_data, res_obj = res_obj, cols = group_by)
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Marker, y = .data[[group_by]])) +
    ggplot2::geom_point(ggplot2::aes(size = Pct_Exp, color = Mean_Exp)) +
    ggplot2::scale_size_continuous(range = c(0, 6), limits = c(0, 100), name = "% Expressing") +
    ggplot2::scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3",
                                   midpoint = 0, name = "Mean Z-Score") +
    theme_publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
                   panel.grid.major = ggplot2::element_line(color = "grey90", linetype = "dashed")) +
    ggplot2::labs(x = "Markers", y = group_by, title = "Expression Dotplot (Bubble Plot)")
  
  message("   [Figure 4] Core-marker dot plot generated")
  return(p)
}


#' All-in-one violin plot grouped by cluster, then group
#'
#' Places all violins in a single panel, ordered as Cluster1_Group1, Cluster1_Group2, ...
Plot_Violin_AllInOne <- function(res_obj, marker, cluster_by = "Cluster_Name", group_by = "Group",
                                 fill_by = NULL, add_boxplot = TRUE, trim = FALSE,
                                 add_stats = TRUE, stat_method = "wilcox.test",
                                 stat_label = "p.signif", comparisons = NULL,
                                 hide_ns = FALSE, violin_alpha = 0.7,
                                 logfc_digits = 2, logfc_prefix = "log2FC=",
                                 p_y_offset = 0.08, logfc_y_offset = 0.16,
                                 point_size = NULL) {
  df <- Apply_Plot_Order(res_obj$Data, res_obj = res_obj, cols = c(cluster_by, group_by, fill_by))

  if (!(marker %in% colnames(df))) stop(paste("  Marker not found in the data:", marker))
  if (!(cluster_by %in% colnames(df))) stop(paste("  Cluster column not found in the data:", cluster_by))
  if (!(group_by %in% colnames(df))) stop(paste("  Group column not found in the data:", group_by))

  if (is.null(fill_by)) fill_by <- group_by
  if (!(fill_by %in% colnames(df))) stop(paste("  Fill column not found in the data:", fill_by))

  df <- df[!is.na(df[[marker]]) & !is.na(df[[cluster_by]]) & !is.na(df[[group_by]]), , drop = FALSE]

  cluster_levels <- if (is.factor(df[[cluster_by]])) levels(df[[cluster_by]]) else unique(as.character(df[[cluster_by]]))
  group_levels <- if (is.factor(df[[group_by]])) levels(df[[group_by]]) else unique(as.character(df[[group_by]]))

  df$PlotGroup <- paste(df[[cluster_by]], df[[group_by]], sep = "\n")
  plot_levels <- unlist(lapply(cluster_levels, function(cl) paste(cl, group_levels, sep = "\n")))
  plot_levels <- plot_levels[plot_levels %in% unique(as.character(df$PlotGroup))]
  df$PlotGroup <- factor(as.character(df$PlotGroup), levels = plot_levels)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[["PlotGroup"]], y = .data[[marker]], fill = .data[[fill_by]])
  ) +
    ggplot2::geom_violin(trim = trim, alpha = violin_alpha, scale = "width") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    ) +
    ggplot2::labs(
      title = paste(marker, "Expression by", cluster_by, "and", group_by),
      x = paste(cluster_by, "and", group_by),
      y = marker
    )

  if (!is.null(point_size)) {
    p <- p + ggplot2::geom_jitter(width = 0.12, size = point_size, alpha = 0.5)
  }

  if (add_boxplot) {
    p <- p + ggplot2::geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA)
  }

  if (add_stats) {
    if (length(group_levels) != 2) {
      stop("  Statistical annotation currently supports exactly 2 groups per cluster.")
    }

    g1 <- group_levels[1]
    g2 <- group_levels[2]
    global_range <- diff(range(df[[marker]], na.rm = TRUE))
    if (!is.finite(global_range) || global_range == 0) global_range <- max(abs(df[[marker]]), na.rm = TRUE)
    if (!is.finite(global_range) || global_range == 0) global_range <- 1

    stat_list <- lapply(cluster_levels, function(cl) {
      sub <- df[df[[cluster_by]] == cl, , drop = FALSE]
      sub <- sub[!is.na(sub[[marker]]) & !is.na(sub[[group_by]]), , drop = FALSE]
      vals1 <- sub[sub[[group_by]] == g1, marker, drop = TRUE]
      vals2 <- sub[sub[[group_by]] == g2, marker, drop = TRUE]
      if (length(vals1) < 1 || length(vals2) < 1) return(NULL)

      test_res <- tryCatch(
        switch(stat_method,
               "t.test" = stats::t.test(vals2, vals1),
               "wilcox.test" = stats::wilcox.test(vals2, vals1),
               stop("  stat_method only supports 'wilcox.test' or 't.test'")),
        error = function(e) NULL
      )
      if (is.null(test_res)) return(NULL)

      pval <- test_res$p.value
      p_label <- if (stat_label == "p.format") {
        paste0("p=", formatC(pval, format = "e", digits = 2))
      } else {
        if (pval <= 0.0001) "****" else if (pval <= 0.001) "***" else if (pval <= 0.01) "**" else if (pval <= 0.05) "*" else "ns"
      }
      if (hide_ns && identical(p_label, "ns")) return(NULL)

      mean1 <- mean(vals1, na.rm = TRUE)
      mean2 <- mean(vals2, na.rm = TRUE)
      logfc <- log2((mean2 + 1e-9) / (mean1 + 1e-9))

      y_base <- max(sub[[marker]], na.rm = TRUE)
      xmin <- paste(cl, g1, sep = "\n")
      xmax <- paste(cl, g2, sep = "\n")
      xmid <- mean(match(c(xmin, xmax), levels(df$PlotGroup)))

      data.frame(
        cluster = cl,
        group1 = xmin,
        group2 = xmax,
        y.position = y_base + global_range * p_y_offset,
        p.label = p_label,
        logFC = logfc,
        logFC_label = paste0(logfc_prefix, sprintf(paste0("%.", logfc_digits, "f"), logfc)),
        xmid = xmid,
        stringsAsFactors = FALSE
      )
    })

    stat_df <- do.call(rbind, stat_list)

    if (!is.null(stat_df) && nrow(stat_df) > 0) {
      p <- p +
        ggpubr::stat_pvalue_manual(
          stat_df,
          label = "p.label",
          xmin = "group1",
          xmax = "group2",
          y.position = "y.position",
          tip.length = 0.01,
          bracket.size = 0.4,
          inherit.aes = FALSE
        ) +
        ggplot2::geom_text(
          data = stat_df,
          ggplot2::aes(x = xmid, y = y.position + global_range * (logfc_y_offset - p_y_offset), label = logFC_label),
          inherit.aes = FALSE,
          size = 3.2,
          fontface = "bold"
        ) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0.02, max(0.08, logfc_y_offset + 0.04)))
        )
    }
  }

  return(p)
}

#' Combined violin plot + box plot
Compare_Populations <- function(res_obj, marker, group_by = "CellType") {
  df <- Apply_Plot_Order(res_obj$Data, res_obj = res_obj, cols = group_by)

  if (!(group_by %in% colnames(df))) stop(paste("  Grouping column not found in the data:", group_by))
  if (!(marker   %in% colnames(df))) stop(paste("  Marker not found in the data:", marker))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[group_by]], y = .data[[marker]], fill = .data[[group_by]])) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.7) +
    ggplot2::geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggpubr::stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
    ggplot2::ggtitle(paste0("Comparison of ", marker))

  return(p)
}

#' Barplot comparison with LogFC annotation
#'
#' Creates a barplot showing mean expression levels across groups with LogFC annotations
#'
#' @param res_obj Result object containing Data
#' @param marker Marker name to compare
#' @param group_by Column name for grouping (default: "CellType")
#' @param reference_group Reference group for LogFC calculation (default: first group)
#' @param show_error_bars Show error bars (default: TRUE, displays SEM)
#' @param logfc_threshold Minimum |LogFC| to display annotation (default: 0.5)
#' @param p_threshold P-value threshold for significance (default: 0.05)
#' @param facet_by Optional column name for faceting (default: NULL, no faceting)
#'
#' @return ggplot object
#' @export
Compare_Populations_Barplot <- function(res_obj, marker, group_by = "CellType",
                                        reference_group = NULL, show_error_bars = TRUE,
                                        logfc_threshold = 0.5, p_threshold = 0.05,
                                        facet_by = NULL) {
  df <- Apply_Plot_Order(res_obj$Data, res_obj = res_obj, cols = group_by)

  if (!(group_by %in% colnames(df))) stop(paste("  Grouping column not found in the data:", group_by))
  if (!(marker   %in% colnames(df))) stop(paste("  Marker not found in the data:", marker))
  if (!is.null(facet_by) && !(facet_by %in% colnames(df))) stop(paste("  Facet column not found in the data:", facet_by))

  # If faceting, calculate stats for each facet group separately
  if (!is.null(facet_by)) {
    facet_levels <- unique(df[[facet_by]])
    stats_list <- list()

    for (facet_val in facet_levels) {
      df_sub <- df[df[[facet_by]] == facet_val, ]

      # Calculate mean and SEM for each group within this facet
      stats_sub <- aggregate(df_sub[[marker]], by = list(Group = df_sub[[group_by]]),
                             FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                 sem = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))),
                                                 n = sum(!is.na(x))))
      stats_sub <- do.call(data.frame, stats_sub)
      colnames(stats_sub) <- c("Group", "Mean", "SEM", "N")
      stats_sub[[facet_by]] <- facet_val

      # Set reference group
      ref_group <- if (is.null(reference_group)) as.character(stats_sub$Group[1]) else reference_group
      if (!(ref_group %in% stats_sub$Group)) {
        warning(paste("Reference group", ref_group, "not found in facet", facet_val, "- using first group"))
        ref_group <- as.character(stats_sub$Group[1])
      }
      ref_mean <- stats_sub$Mean[stats_sub$Group == ref_group]

      # Calculate LogFC and p-values
      stats_sub$LogFC <- log2((stats_sub$Mean + 0.01) / (ref_mean + 0.01))
      stats_sub$P_value <- NA

      for (i in seq_len(nrow(stats_sub))) {
        if (stats_sub$Group[i] != ref_group) {
          group_data <- df_sub[[marker]][df_sub[[group_by]] == stats_sub$Group[i]]
          ref_data <- df_sub[[marker]][df_sub[[group_by]] == ref_group]
          test_result <- tryCatch(t.test(group_data, ref_data), error = function(e) NULL)
          if (!is.null(test_result)) stats_sub$P_value[i] <- test_result$p.value
        }
      }

      stats_list[[as.character(facet_val)]] <- stats_sub
    }

    stats <- do.call(rbind, stats_list)
    rownames(stats) <- NULL

  } else {
    # No faceting - original logic
    stats <- aggregate(df[[marker]], by = list(Group = df[[group_by]]),
                       FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                           sem = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))),
                                           n = sum(!is.na(x))))
    stats <- do.call(data.frame, stats)
    colnames(stats) <- c("Group", "Mean", "SEM", "N")

    # Set reference group
    if (is.null(reference_group)) {
      reference_group <- as.character(stats$Group[1])
    }
    if (!(reference_group %in% stats$Group)) {
      stop(paste("  Reference group not found:", reference_group))
    }
    ref_mean <- stats$Mean[stats$Group == reference_group]

    # Calculate LogFC and p-values
    stats$LogFC <- log2((stats$Mean + 0.01) / (ref_mean + 0.01))
    stats$P_value <- NA

    for (i in seq_len(nrow(stats))) {
      if (stats$Group[i] != reference_group) {
        group_data <- df[[marker]][df[[group_by]] == stats$Group[i]]
        ref_data <- df[[marker]][df[[group_by]] == reference_group]
        test_result <- tryCatch(t.test(group_data, ref_data), error = function(e) NULL)
        if (!is.null(test_result)) stats$P_value[i] <- test_result$p.value
      }
    }
  }

  stats <- Apply_Plot_Order(stats, res_obj = res_obj, cols = c("Group", facet_by),
                            order_list = setNames(list(levels(df[[group_by]])), "Group"))

  # Create significance labels
  stats$Significance <- ifelse(is.na(stats$P_value), "",
                               ifelse(stats$P_value < 0.001, "***",
                                      ifelse(stats$P_value < 0.01, "**",
                                             ifelse(stats$P_value < 0.05, "*", "ns"))))

  # Create LogFC labels (only for significant and above threshold)
  stats$LogFC_label <- ifelse(
    !is.na(stats$P_value) & stats$P_value < p_threshold & abs(stats$LogFC) > logfc_threshold,
    sprintf("LogFC=%.2f\n%s", stats$LogFC, stats$Significance),
    ""
  )

  # Mark reference group (only if not faceting, to avoid confusion)
  if (is.null(facet_by)) {
    stats$Group_label <- ifelse(stats$Group == reference_group,
                                paste0(as.character(stats$Group), "
(Ref)"),
                                as.character(stats$Group))
  } else {
    stats$Group_label <- as.character(stats$Group)
  }

  # Preserve plotting order on x-axis after relabeling
  stats$Group_label <- factor(stats$Group_label, levels = unique(stats$Group_label))
  stats$Group <- factor(as.character(stats$Group), levels = unique(as.character(stats$Group)))

  # Create barplot
  p <- ggplot2::ggplot(stats, ggplot2::aes(x = Group_label, y = Mean, fill = Group)) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3)

  # Add error bars if requested
  if (show_error_bars) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - SEM, ymax = Mean + SEM),
                                     width = 0.3, linewidth = 0.5)
  }

  # Add LogFC annotations
  p <- p +
    ggplot2::geom_text(ggplot2::aes(label = LogFC_label, y = Mean + SEM + max(Mean) * 0.05),
                       size = 3.5, fontface = "bold", vjust = 0) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    ggplot2::labs(
      x = group_by,
      y = paste0("Mean ", marker, " Expression"),
      title = paste0("Comparison of ", marker, " across ", group_by)
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15)))

  # Add faceting if requested
  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)), scales = "free_y")
  }

  return(p)
}


#' Expression-profile heatmap
Plot_Heatmap <- function(res_obj, group_by = "Cluster", scale = "row",
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         color_palette = NULL, show_border = TRUE, fontsize = 10) {
  df <- Apply_Plot_Order(res_obj$Data, res_obj = res_obj, cols = group_by)
  
  meta_cols <- c("tSNE1", "tSNE2", "UMAP1", "UMAP2", "Cluster", "Group", "SampleID",
                 "Time", "TIME", "CellType", "id", "Cluster_Name", "Gate_Status", "Pseudotime")
  mk <- setdiff(colnames(df)[sapply(df, is.numeric)], meta_cols)
  
  hm_data <- aggregate(df[, mk], by = list(Group = df[[group_by]]), FUN = median, na.rm = TRUE)
  hm_data <- Apply_Plot_Order(hm_data, res_obj = res_obj, cols = "Group",
                              order_list = setNames(list(levels(df[[group_by]])), "Group"))
  hm_data <- hm_data[order(hm_data$Group), , drop = FALSE]
  rownames(hm_data) <- hm_data$Group
  mat <- t(as.matrix(hm_data[, -1]))
  
  pal <- if (is.null(color_palette)) {
    colorRampPalette(c("navy", "white", "firebrick3"))(100)
  } else {
    color_palette
  }
  
  p <- pheatmap::pheatmap(
    mat,
    scale        = scale,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    color        = pal,
    border_color = if (show_border) "white" else NA,
    fontsize     = fontsize,
    angle_col    = 45,
    main         = paste0("Heatmap (Scale: ", scale, ") - ", group_by)
  )
  message("   [Figure 6] Aggregated subpopulation expression heatmap generated")
  return(list(Plot = p, Matrix = mat))
}

# ==============================================================================
# 4. Model evaluation and quality control (Confusion Matrix & QC)
# ==============================================================================

#' Confusion-matrix heatmap
Plot_Confusion_Matrix <- function(res_obj, true_label = "CellType", pred_label = "Cluster_Name") {
  df <- Apply_Plot_Order(res_obj$Data, res_obj = res_obj, cols = c(true_label, pred_label))
  if (!(true_label %in% colnames(df)) || !(pred_label %in% colnames(df))) {
    stop("  Unable to find the corresponding label columns")
  }
  
  conf_mat      <- table(Truth = df[[true_label]], Prediction = df[[pred_label]])
  conf_mat_prop <- prop.table(conf_mat, margin = 1) * 100
  melted_conf   <- reshape2::melt(conf_mat_prop)
  melted_conf$Truth <- factor(as.character(melted_conf$Truth), levels = rownames(conf_mat_prop))
  melted_conf$Prediction <- factor(as.character(melted_conf$Prediction), levels = colnames(conf_mat_prop))
  
  p <- ggplot2::ggplot(melted_conf, ggplot2::aes(x = Prediction, y = Truth, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", value)),
                       color = ifelse(melted_conf$value > 50, "white", "black"), size = 3.5) +
    ggplot2::scale_fill_gradient(low = "white", high = "#2C3E50", name = "% Match") +
    theme_publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggplot2::labs(x = paste("Predicted:", pred_label),
                  y = paste("Ground Truth:", true_label),
                  title = "Confusion Matrix")
  
  message("   [Figure 7] AI-predicted classification confusion-matrix heatmap generated")
  return(p)
}

#' Ridge plot / Joy plot (Ridge Plot)
Plot_Density_Line <- function(res_obj, marker, group_by = "Cluster_Name", custom_colors = NULL) {
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    stop("Package 'ggridges' is required. Install it with install.packages('ggridges').")
  }
  
  df <- Apply_Plot_Order(res_obj$Data, res_obj = res_obj, cols = group_by)
  if (!(marker   %in% colnames(df))) stop(paste("  Marker not found in the data:", marker))
  if (!(group_by %in% colnames(df))) stop(paste("  Grouping column not found in the data:", group_by))
  
  df <- df[!is.na(df[[group_by]]), ]
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[marker]], y = .data[[group_by]], fill = .data[[group_by]])) +
    ggridges::geom_density_ridges(alpha = 0.8, scale = 1.5, rel_min_height = 0.01,
                                  color = "white", linewidth = 0.5) +
    theme_publication() +
    ggplot2::theme(legend.position = "none", panel.grid.major.y = ggplot2::element_blank()) +
    ggplot2::labs(x = paste(marker, "Expression (arcsinh)"), y = group_by,
                  title = paste("Ridge Plot of", marker, "by", group_by))
  
  if (!is.null(custom_colors)) {
    p <- p + ggplot2::scale_fill_manual(values = custom_colors)
  } else {
    n_colors <- length(unique(df[[group_by]]))
    safe_n   <- max(3, min(9, n_colors))
    pal      <- colorRampPalette(RColorBrewer::brewer.pal(safe_n, "Set1"))(n_colors)
    p        <- p + ggplot2::scale_fill_manual(values = pal)
  }
  message("   [Figure 8] Ridge plot generated")
  return(p)
}

#' Scatter plot with marginal density histograms
Plot_Scatterhist_Line <- function(res_obj,
                                  x_marker,
                                  y_marker,
                                  group_by = "Cluster_Name",
                                  custom_colors = NULL,
                                  point_size = 0.5,
                                  point_alpha = 0.6,
                                  density_linewidth = 0.6,
                                  margin_ratio = 0.15,   #          0.15         15%
                                  density_height = 0.8)  #         0.8            80%
{
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package 'patchwork' is required.")
  if (!requireNamespace("cowplot", quietly = TRUE)) stop("Package 'cowplot' is required.")
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Package 'RColorBrewer' is required.")
  
  df <- res_obj$Data
  df <- Apply_Plot_Order(df, res_obj = res_obj, cols = c(group_by))
  
  if (!(x_marker %in% colnames(df))) stop(paste("  x_marker not found:", x_marker))
  if (!(y_marker %in% colnames(df))) stop(paste("  y_marker not found:", y_marker))
  if (!(group_by %in% colnames(df))) stop(paste("  group_by not found:", group_by))
  
  df <- df[!is.na(df[[x_marker]]) & !is.na(df[[y_marker]]) & !is.na(df[[group_by]]), , drop = FALSE]
  if (nrow(df) == 0) stop("  No valid rows remaining after removing NA values.")
  
  grp_levels <- if (is.factor(df[[group_by]])) {
    levels(df[[group_by]])
  } else {
    unique(as.character(df[[group_by]]))
  }
  grp_levels <- grp_levels[grp_levels %in% unique(as.character(df[[group_by]]))]
  df[[group_by]] <- factor(as.character(df[[group_by]]), levels = grp_levels)
  
  n_colors <- length(grp_levels)
  pal <- if (!is.null(custom_colors)) {
    custom_colors
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(min(12, n_colors), "Set3"))(n_colors)
  }
  
  if (is.null(names(pal))) {
    names(pal) <- grp_levels[seq_len(min(length(pal), length(grp_levels)))]
  }
  pal <- pal[grp_levels]
  
  # ==================== 1.      ====================
  p_main_legend <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[x_marker]],
      y = .data[[y_marker]],
      color = .data[[group_by]]
    )
  ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_color_manual(values = pal, name = group_by) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(x = x_marker, y = y_marker)
  
  legend <- cowplot::get_legend(p_main_legend)
  p_main <- p_main_legend + ggplot2::theme(legend.position = "none")
  
  # ==================== 2.        ====================
  p_top <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[x_marker]],
      y = ggplot2::after_stat(scaled),
      color = .data[[group_by]]
    )
  ) +
    ggplot2::stat_density(
      geom = "line",
      position = "identity",
      linewidth = density_linewidth
    ) +
    ggplot2::scale_color_manual(values = pal, guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, 1 / density_height), expand = c(0, 0)) +
    ggplot2::theme_void()
  
  # ==================== 3.        (       ) ====================
  p_right <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[y_marker]],             #         X     
      y = ggplot2::after_stat(scaled),   # Y          
      color = .data[[group_by]]
    )
  ) +
    ggplot2::stat_density(
      geom = "line",
      position = "identity",
      linewidth = density_linewidth
    ) +
    ggplot2::scale_color_manual(values = pal, guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, 1 / density_height), expand = c(0, 0)) +
    ggplot2::coord_flip() +              #                
    ggplot2::theme_void()
  
  # ==================== 4.        ====================
  # A  p_main, B  p_top, C  p_right, #   
  layout_design <- "
    B#
    AC
  "
  p_core <- p_main + p_top + p_right +
    patchwork::plot_layout(
      design = layout_design,
      widths = c(1, margin_ratio),      #        C   A    
      heights = c(margin_ratio, 1)      #        B   A    
    )
  
  # ==================== 5.      ====================
  p_final <- cowplot::plot_grid(
    p_core,
    legend,
    ncol = 1,
    rel_heights = c(1, 0.12)
  )
  
  message("   [Figure] Dual-marker marginal-density scatter plot generated successfully")
  return(p_final)
}
# ==============================================================================
# 5. Macro sample-level analysis (Sample-Level PCA)
# ==============================================================================

#' Sample-level PCA visualization
Plot_Sample_PCA <- function(res_obj, color_by = "Group", label_by = "SampleID") {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required. Install it with install.packages('ggrepel').")
  }
  
  df <- Apply_Plot_Order(res_obj$Data, res_obj = res_obj, cols = c(color_by, label_by))
  if (!(color_by %in% colnames(df))) stop(paste("  Grouping column not found in the data:", color_by))
  if (!(label_by %in% colnames(df))) stop(paste("  Sample identifier column not found in the data:", label_by))
  
  message("\n>>> Aggregating single-cell data and calculating sample-level features...")
  
  meta_cols    <- c("tSNE1", "tSNE2", "UMAP1", "UMAP2", "Cluster", "Group", "SampleID",
                    "Time", "TIME", "CellType", "id", "Cluster_Name", "Gate_Status",
                    "Pseudotime", "PC1", "PC2")
  numeric_cols <- colnames(df)[sapply(df, is.numeric)]
  markers      <- setdiff(numeric_cols, meta_cols)
  
  sample_df <- df |>
    dplyr::group_by(.data[[label_by]], .data[[color_by]]) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(markers), stats::median, na.rm = TRUE),
                     .groups = "drop")
  
  if (nrow(sample_df) < 3) {
    warning("   Fewer than 3 samples; unable to produce a meaningful PCA plot.")
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("Not enough samples for PCA"))
  }
  
  message(">>> Running principal component analysis (PCA)...")
  mat         <- as.matrix(sample_df[, markers])
  pca_res     <- stats::prcomp(mat, scale. = TRUE)
  var_explained <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 1)
  
  plot_data <- data.frame(
    Sample = sample_df[[label_by]],
    Group  = sample_df[[color_by]],
    PC1    = pca_res$x[, 1],
    PC2    = pca_res$x[, 2]
  )
  plot_data <- Apply_Plot_Order(plot_data, res_obj = res_obj, cols = c("Group", "Sample"),
                                order_list = setNames(list(levels(df[[color_by]]), levels(df[[label_by]])), c("Group", "Sample")))
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    ggplot2::geom_point(size = 5, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 4, show.legend = FALSE, fontface = "bold") +
    theme_publication() +
    ggplot2::labs(
      x     = paste0("PC1 (", var_explained[1], "%)"),
      y     = paste0("PC2 (", var_explained[2], "%)"),
      title = "Sample-Level PCA"
    )
  
  n_colors <- length(unique(plot_data$Group))
  safe_n   <- max(3, min(9, n_colors))
  pal      <- colorRampPalette(RColorBrewer::brewer.pal(safe_n, "Set1"))(n_colors)
  p        <- p + ggplot2::scale_color_manual(values = pal)
  
  message("   [Figure 10] Sample-level PCA plot generated")
  return(p)
}

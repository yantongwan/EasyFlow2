#' 快速读取并合并 FCS 文件
#'
#' @param work_dir FCS 文件夹路径
#' @param group_keywords 分组关键词 (如 c("Ctrl", "Treat"))
#' @param n_cells_per_sample 每个样本抽样数 (建议 5000-50000)
#' @import flowCore
#' @import dplyr
#' @export
Prepare_Raw_Data <- function(work_dir, group_keywords, n_cells_per_sample=5000) {
  files <- list.files(work_dir, pattern = ".fcs$", full.names = TRUE)
  if(length(files) == 0) stop("❌ 未找到 FCS 文件")

  combined_data <- data.frame()
  message(">>> 正在合并 ", length(files), " 个文件...")

  for (file in files) {
    fcs <- flowCore::read.FCS(file, transformation = FALSE, truncate_max_range = FALSE)
    raw_exprs <- flowCore::exprs(fcs)
    exprs_data_trans <- asinh(raw_exprs / 150)

    use_cols <- grep("-A$", colnames(exprs_data_trans), value = TRUE)
    use_cols <- use_cols[!grepl("Time|TIME", use_cols)]
    data_clean <- as.data.frame(exprs_data_trans[, use_cols])

    if (nrow(data_clean) > n_cells_per_sample) {
      data_sub <- data_clean[sample(nrow(data_clean), n_cells_per_sample), ]
    } else {
      data_sub <- data_clean
    }

    data_sub$SampleID <- basename(file)
    data_sub$Group <- "Other"
    for (key in group_keywords) {
      if (grepl(key, basename(file), ignore.case = TRUE)) {
        data_sub$Group <- key; break
      }
    }
    combined_data <- dplyr::bind_rows(combined_data, data_sub)
  }
  row.names(combined_data) <- paste0("Cell_", 1:nrow(combined_data))
  return(combined_data)
}

#' 应用蛋白名称映射 (批量重命名通道)
#'
#' 将原始的荧光通道名 (如 PE-Cy7-A) 替换为生物学蛋白名 (如 CD4)。
#'
#' @param df 原始数据框 (由 Prepare_Raw_Data 生成)
#' @param panel_list 一个命名向量，格式为 c("通道名" = "蛋白名", ...)
#' @return 返回列名已更新的数据框
#' @export
Apply_Panel_Mapping <- function(df, panel_list) {
  message(">>> 🔄 正在应用 Panel 映射...")

  # 获取当前列名
  current_cols <- colnames(df)
  new_cols <- current_cols

  found_count <- 0

  # 遍历映射表进行替换
  for (channel in names(panel_list)) {
    # 精确匹配
    if (channel %in% current_cols) {
      idx <- which(current_cols == channel)
      protein_name <- panel_list[[channel]]
      new_cols[idx] <- protein_name
      found_count <- found_count + 1
    } else {
      warning(paste("⚠️ 没找到通道:", channel, "，跳过。"))
    }
  }

  colnames(df) <- new_cols
  message(paste0(">>> ✅ 成功替换了 ", found_count, " 个通道名称！"))
  return(df)
}

# ==============================================================================
# 模块 2: 模式 A - 全自动分析
# ==============================================================================

#' 启动全自动分析流程
#' @param work_dir FCS 文件夹路径
#' @param group_keywords 分组关键词
#' @param n_cells_per_sample 每个样本抽样数
#' @param k_neighbors 聚类精度
#' @export
Run_Automated_Analysis <- function(work_dir, group_keywords, n_cells_per_sample=3000, k_neighbors=30) {
  message(">>> 🤖 [模式 A] 启动全自动分析流程...")
  raw_data <- Prepare_Raw_Data(work_dir, group_keywords, n_cells_per_sample)
  results <- Run_Analysis_On_Gated_Data(raw_data, k_neighbors)
  results$Info$Method <- "Automated_Random_Sampling"
  message(">>> 🎉 全自动分析完成！")
  return(results)
}

# ==============================================================================
# 模块 3: 模式 B - 单群过滤
# ==============================================================================

#' 启动单群过滤 APP
#' @param merged_df 原始数据
#' @import shiny
#' @import plotly
#' @import miniUI
#' @export
Run_Interactive_Gating <- function(merged_df) {

  numeric_cols <- colnames(merged_df)[sapply(merged_df, is.numeric)]
  numeric_cols <- numeric_cols[!grepl("tSNE|Cluster", numeric_cols)]

  ui <- miniPage(
    gadgetTitleBar("模式 B: 单群过滤 (Single Filter)"),
    miniTabstripPanel(
      miniTabPanel("Gating", icon = icon("filter"),
                   miniContentPanel(
                     fillCol(flex = c(1, 4),
                             fillRow(
                               div(selectInput("x_axis", "X 轴:", choices = numeric_cols, selected = numeric_cols[1]),
                                   selectInput("y_axis", "Y 轴:", choices = numeric_cols, selected = numeric_cols[2]),
                                   style = "padding-right: 10px;"),
                               div(h4(textOutput("cell_count_text")),
                                   actionButton("filter_btn", "✂️ 保留圈中细胞", class = "btn-danger"),
                                   actionButton("reset_btn", "🔄 重置", class = "btn-warning"),
                                   helpText("提示：此模式用于层层筛选，最后只保留剩下的一种细胞。"),
                                   style = "padding-left: 10px;")
                             ),
                             plotlyOutput("scatter_plot", height = "100%")
                     )
                   )
      )
    )
  )

  server <- function(input, output, session) {
    values <- reactiveValues(original_data = merged_df, current_data = merged_df)

    output$scatter_plot <- renderPlotly({
      req(input$x_axis, input$y_axis)
      plot_data <- values$current_data
      plot_ly(data = plot_data, x = ~get(input$x_axis), y = ~get(input$y_axis), key = ~row.names(plot_data),
              type = 'scattergl', mode = 'markers', marker = list(size = 2, opacity = 0.5, color = '#3366CC'), source = "subset_plot") %>%
        layout(xaxis = list(title = input$x_axis), yaxis = list(title = input$y_axis), dragmode = "lasso") %>% toWebGL()
    })

    output$cell_count_text <- renderText({
      paste0("当前剩余: ", nrow(values$current_data), " (", round(nrow(values$current_data)/nrow(values$original_data)*100, 1), "%)")
    })

    observeEvent(input$filter_btn, {
      event_data <- event_data("plotly_selected", source = "subset_plot")
      if (is.null(event_data)) { showNotification("⚠️ 请先圈选！", type = "warning"); return() }
      values$current_data <- values$current_data[row.names(values$current_data) %in% as.character(event_data$key), ]
      showNotification(paste0("✅ 已筛选，剩余 ", nrow(values$current_data)), type = "message")
    })

    observeEvent(input$reset_btn, { values$current_data <- values$original_data })
    observeEvent(input$done, { stopApp(values$current_data) })
  }
  runGadget(ui, server, viewer = dialogViewer("Single Gating", width = 1000, height = 800))
}

# ==============================================================================
# 模块 3.5: 模式 C - 多群选择
# ==============================================================================

#' 启动多群选择 APP
#' @param merged_df 原始数据
#' @export
Run_Multiclass_Gating <- function(merged_df) {

  numeric_cols <- colnames(merged_df)[sapply(merged_df, is.numeric)]
  numeric_cols <- numeric_cols[!grepl("tSNE|Cluster", numeric_cols)]

  ui <- miniPage(
    gadgetTitleBar("模式 C: 多群选择 (Multi-Class Selection)"),
    miniTabstripPanel(
      miniTabPanel("Gating", icon = icon("crosshairs"),
                   miniContentPanel(
                     fillCol(flex = c(1, 4),
                             fillRow(
                               div(selectInput("x_axis", "X 轴:", choices = numeric_cols, selected = numeric_cols[1]),
                                   selectInput("y_axis", "Y 轴:", choices = numeric_cols, selected = numeric_cols[2]),
                                   style = "padding-right: 5px;"),
                               div(textInput("group_name", "群体名称 (如 CD4):", value = "Population_1"),
                                   actionButton("add_btn", "➕ 添加", class = "btn-success"),
                                   actionButton("filter_btn", "✂️ 仅保留", class = "btn-danger"),
                                   actionButton("reset_view_btn", "🔄 重置视图", class = "btn-warning"),
                                   actionButton("clear_list_btn", "🗑️ 清空列表", class = "btn-default"),
                                   style = "padding-left: 5px;")
                             ),
                             fillRow(flex = c(3, 1),
                                     plotlyOutput("scatter_plot", height = "100%"),
                                     div(h4("已保存的群体:"), helpText("点击'重置视图'不会删除这里的列表"), tableOutput("gated_list_table"),
                                         style = "overflow-y: scroll; background-color: #f5f5f5; padding: 10px;")
                             )
                     )
                   )
      )
    )
  )

  server <- function(input, output, session) {
    values <- reactiveValues(original_data = merged_df, current_data = merged_df, saved_populations = list(), saved_info = data.frame(Name=character(), Count=integer()))

    output$scatter_plot <- renderPlotly({
      req(input$x_axis, input$y_axis)
      plot_data <- values$current_data
      plot_ly(data = plot_data, x = ~get(input$x_axis), y = ~get(input$y_axis), key = ~row.names(plot_data),
              type = 'scattergl', mode = 'markers', marker = list(size = 2, opacity = 0.5, color = '#999999'), source = "subset_plot") %>%
        layout(xaxis = list(title = input$x_axis), yaxis = list(title = input$y_axis), dragmode = "lasso") %>% toWebGL()
    })

    observeEvent(input$add_btn, {
      event_data <- event_data("plotly_selected", source = "subset_plot")
      if (is.null(event_data) || input$group_name == "") { showNotification("⚠️ 请圈选并输入名称", type = "warning"); return() }
      subset_df <- values$current_data[row.names(values$current_data) %in% as.character(event_data$key), ]
      subset_df$CellType <- input$group_name
      values$saved_populations[[length(values$saved_populations) + 1]] <- subset_df
      values$saved_info <- rbind(values$saved_info, data.frame(Name = input$group_name, Count = nrow(subset_df)))
      showNotification(paste0("✅ 已添加: ", input$group_name), type = "message")
    })

    observeEvent(input$filter_btn, {
      event_data <- event_data("plotly_selected", source = "subset_plot")
      if (!is.null(event_data)) {
        values$current_data <- values$current_data[row.names(values$current_data) %in% as.character(event_data$key), ]
        showNotification("✂️ 视图已过滤", type = "warning")
      }
    })

    observeEvent(input$reset_view_btn, { values$current_data <- values$original_data; showNotification("🔄 视图已复原", type = "message") })
    observeEvent(input$clear_list_btn, { values$saved_populations <- list(); values$saved_info <- data.frame(Name=character(), Count=integer()); showNotification("🗑️ 已清空", type = "warning") })
    output$gated_list_table <- renderTable({ values$saved_info })
    observeEvent(input$done, {
      if (length(values$saved_populations) == 0) stopApp(NULL)
      else stopApp(dplyr::bind_rows(values$saved_populations))
    })
  }
  runGadget(ui, server, viewer = dialogViewer("Multi-Class Gating", width = 1100, height = 800))
}

# ==============================================================================
# 模块 4: 核心分析引擎
# ==============================================================================

#' 自动分析圈好的数据
#' @param gated_data 圈门后的数据
#' @param k_neighbors 聚类精度
#' @import Rtsne
#' @import igraph
#' @import RANN
#' @export
Run_Analysis_On_Gated_Data <- function(gated_data, k_neighbors = 30) {
  message(">>> 🚀 [分析引擎] 正在分析 ", nrow(gated_data), " 个细胞...")

  ignore_cols <- c("tSNE1", "tSNE2", "Cluster", "SampleID", "Group", "Time", "TIME", "CellType")
  numeric_cols <- colnames(gated_data)[sapply(gated_data, is.numeric)]
  marker_cols <- setdiff(numeric_cols, ignore_cols)

  message(paste0(">>> 🕸️ 构建 KNN 图 (k=", k_neighbors, ")..."))
  cluster_mat <- gated_data[, marker_cols]
  knn_res <- RANN::nn2(data = cluster_mat, k = k_neighbors)
  edges <- matrix(c(rep(1:nrow(cluster_mat), each = k_neighbors), as.vector(t(knn_res$nn.idx))), ncol = 2)
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  set.seed(123); louvain_res <- igraph::cluster_louvain(g)

  gated_data$Cluster <- as.factor(louvain_res$membership)
  message(">>> 🔍 识别出 Cluster 数量: ", length(levels(gated_data$Cluster)))

  message(">>> ☕️ 正在运行 tSNE...")
  set.seed(42)
  perplex <- ifelse(nrow(gated_data) < 91, floor((nrow(gated_data) - 1) / 3), 30)
  tsne_out <- Rtsne(cluster_mat, dims = 2, perplexity = perplex, verbose = FALSE, check_duplicates = FALSE)
  gated_data$tSNE1 <- tsne_out$Y[, 1]
  gated_data$tSNE2 <- tsne_out$Y[, 2]

  message(">>> 🎉 分析完成！")
  return(list(Data = gated_data, Info = list(K = k_neighbors)))
}

# ==============================================================================
# 模块 5: 绘图函数集 (双阳图 / 多群对比)
# ==============================================================================

#' 双阳图绘制
#' @import ggplot2
#' @import scales
#' @export
plot_dual_color <- function(res_obj, marker1, marker2) {
  df <- res_obj$Data
  check <- function(m) { if(m %in% colnames(df)) return(m); m2<-gsub("-",".",m); if(m2 %in% colnames(df)) return(m2); stop("❌ 找不到 Marker") }
  m1 <- check(marker1); m2 <- check(marker2)

  v1 <- scales::rescale(df[[m1]], to=c(0,1))^1.2
  v2 <- scales::rescale(df[[m2]], to=c(0,1))^1.2
  cols <- rgb(red = v1, green = 0, blue = v2, alpha = 0.6)

  ggplot(df, aes(tSNE1, tSNE2)) + geom_point(size=0.5, color=cols) + theme_bw() +
    theme(panel.grid=element_blank()) + ggtitle(paste0(m1, " (Red) vs ", m2, " (Blue)"))
}

#' 多群比较绘图
#' @import ggplot2
#' @import ggpubr
#' @import RColorBrewer
#' @export
Compare_Populations <- function(gated_data, marker) {
  if(!"CellType" %in% colnames(gated_data)) stop("❌ 无 CellType 列")
  check <- function(m) { if(m %in% colnames(gated_data)) return(m); m2<-gsub("-",".",m); if(m2 %in% colnames(gated_data)) return(m2); stop("❌ 找不到 Marker") }
  m <- check(marker)

  ggplot(gated_data, aes(x = CellType, y = .data[[m]], fill = CellType)) +
    geom_violin(trim=FALSE, alpha=0.7) + geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
    theme_bw() + scale_fill_brewer(palette="Set2") + stat_compare_means(label="p.signif", method="t.test", ref.group=".all.") +
    ggtitle(paste0("Comparison of ", m))
}

# ==============================================================================
# 模块 6: 峰峦图与散点图 (重叠密度与双向散点)
# ==============================================================================

#' 绘制重叠密度图 (Y轴为细胞量)
#' @import ggplot2
#' @export
Plot_Density_Line <- function(res_obj, marker, group_by = "Cluster", custom_colors = NULL) {
  df <- res_obj$Data
  check <- function(m) { if(m %in% colnames(df)) return(m); m2<-gsub("-",".",m); if(m2 %in% colnames(df)) return(m2); stop("❌ 找不到 Marker") }
  real_m <- check(marker)

  safe_m <- gsub("-", "_", real_m)
  df[[safe_m]] <- df[[real_m]]

  p <- ggplot(df, aes(x = .data[[safe_m]], color = .data[[group_by]])) +
    geom_density(size = 1, adjust = 1.5) + theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(color = "black"), legend.position = "right") +
    labs(y = "Density (Cell Amount)", x = real_m, title = paste0("Density: ", real_m, " by ", group_by))

  if (!is.null(custom_colors)) p <- p + scale_color_manual(values = custom_colors)
  return(p)
}

#' 绘制双向边缘密度散点图 (无填充)
#' @import ggpubr
#' @export
Plot_Scatterhist_Line <- function(res_obj, x_marker, y_marker, group_by = "Cluster", custom_colors = NULL) {
  df <- res_obj$Data
  check <- function(m) { if(m %in% colnames(df)) return(m); m2<-gsub("-",".",m); if(m2 %in% colnames(df)) return(m2); stop("❌ 找不到 Marker") }
  real_mx <- check(x_marker); real_my <- check(y_marker)

  safe_mx <- gsub("-", "_", real_mx)
  safe_my <- gsub("-", "_", real_my)
  df[[safe_mx]] <- df[[real_mx]]
  df[[safe_my]] <- df[[real_my]]

  p <- ggpubr::ggscatterhist(
    df, x = safe_mx, y = safe_my, color = group_by,
    size = 0.5, alpha = 0.6,
    palette = if(is.null(custom_colors)) "jco" else custom_colors,
    margin.plot = "density",
    margin.params = list(fill = NA, color = group_by, size = 0.8),
    xlab = real_mx, ylab = real_my,
    ggtheme = theme_bw() + theme(panel.grid = element_blank())
  )
  return(p)
}

# ==============================================================================
# 模块 7: 热图绘制
# ==============================================================================

#' 绘制全能热图 (支持自定义参数)
#' @import dplyr
#' @import pheatmap
#' @import RColorBrewer
#' @export
Plot_Heatmap <- function(res_obj, group_by = "Cluster", scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, color_palette = NULL, show_border = TRUE, fontsize = 10) {
  df <- res_obj$Data
  meta_cols <- c("tSNE1", "tSNE2", "Cluster", "Group", "SampleID", "Time", "TIME", "CellType", "id")
  numeric_cols <- colnames(df)[sapply(df, is.numeric)]
  marker_cols <- setdiff(numeric_cols, meta_cols)

  message(">>> 🌡️ 生成热图 Markers: ", paste(marker_cols, collapse=", "))
  safe_markers <- paste0("`", marker_cols, "`")
  formula_str <- paste("cbind(", paste(safe_markers, collapse=","), ") ~", group_by)
  heatmap_data <- aggregate(as.formula(formula_str), data = df, FUN = median)

  rownames(heatmap_data) <- heatmap_data[[1]]
  heatmap_mat <- as.matrix(heatmap_data[, -1])
  heatmap_mat_t <- t(heatmap_mat)

  if (is.null(color_palette)) {
    if (scale == "none") color_palette <- colorRampPalette(c("#F7F7F7", "#D6604D", "#B2182B"))(100)
    else color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  }
  border_col <- if(show_border) "white" else NA

  p <- pheatmap::pheatmap(
    heatmap_mat_t, scale = scale,
    cluster_rows = cluster_rows, cluster_cols = cluster_cols,
    color = color_palette, border_color = border_col, fontsize = fontsize, angle_col = 45,
    main = paste0("Heatmap (Scale: ", scale, ") - ", group_by)
  )
  return(list(Plot = p, Matrix = heatmap_mat_t))
}



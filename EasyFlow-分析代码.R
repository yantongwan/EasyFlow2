# 设置工作目录 (请修改路径)
getwd()
setwd("/Users/wanyantong/Desktop/EasyFlow2")
# 1. 安装 Bioconductor 包 (flowCore)
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("flowCore")) BiocManager::install("flowCore", update = FALSE, ask = FALSE)
# 2. 安装本地包 (EasyFlow2)
if (!require("EasyFlow2")) install.packages("EasyFlow2_0.0.0.9000.tar.gz", repos = NULL)
# 3. 批量安装 CRAN 包
pkgs <- c("shiny", "plotly", "miniUI", "dplyr", "Rtsne", "ggplot2", 
          "stringr", "ggpubr", "scales", "RColorBrewer", "igraph", "RANN","ggridges","pheatmap")
# 自动筛选未安装的包并下载
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)
# 4. 一键加载所有包
invisible(lapply(c(pkgs, "flowCore", "EasyFlow2"), library, character.only = TRUE))

# 5. 设置路径
my_folder <- "/Users/wanyantong/Desktop/20251129/LLY 20251129 cytokine PBM/colon"

# 6. 注释数据
my_panel <- c(
  "BV605-A" = "CD4",
  "eFluor450-A" = "CD11b",
  "BV711-A" = "CD45",
  "BV510-A" = "CD11c",
  "PE-Dazzle594-A" = "Ly6c",
  "PerCP-Cy5.5-A" = "Ly6g",
  "PE-A" = "GM-CSF",
  "PE-Cy7-A" = "IL17a",
  "BV421-A" = "TNFa",
  "PE-Cy5-A" = "T-bet",
  "APC-A" = "IFNg",
  "[live/dead fv-ef780]-A" = "L/D" 
)
#模式A，如果做模式A可以不用下面的代码，直接绘图（如果做这一步就不用进行下面的分析了）
final_results <- Run_Automated_Analysis(
  work_dir = my_folder, 
  group_keywords = c("pc", "pl"), 
  n_cells_per_sample = 3000
)


# 读取数据 (设置大一点以免漏掉稀有细胞) -> 下面模式B和C选择其中一个分析
raw_data <- Prepare_Raw_Data(my_folder, c("Con", "PBM"), n_cells_per_sample = 30000)
raw_data <- Apply_Panel_Mapping(raw_data, my_panel)
# 模式B分析 -> 运行 APP -> 圈单个细胞群
single_gated_cells <- Run_Interactive_Gating(raw_data)

# 模式B分析
final_results <- Run_Analysis_On_Gated_Data(single_gated_cells, 
                                            k_neighbors = 30) #k值越大聚类越少

# 模式C分析 -> 运行 APP -> 圈多个细胞群
multi_gated_cells <- Run_Multiclass_Gating(raw_data)

# 模式C分析
final_results <- Run_Analysis_On_Gated_Data(multi_gated_cells, 
                                            k_neighbors = 30) #k值越大聚类越少

# ===========================================================
# 🎨 1. t-SNE 降维图系列
# ===========================================================

# [1.1] 按 Cluster 上色 (查看亚群)
colnames(final_results$Data)
n_clusters <- length(unique(final_results$Data$Cluster))
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_clusters)

p <- ggplot(final_results$Data, aes(tSNE1, tSNE2, color = Cluster)) +
  geom_point(size = 0.1, alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = my_colors) +
  ggtitle("t-SNE by Cluster") +
  guides(color = guide_legend(override.aes = list(size = 4)))

p

ggsave("Cluster_tSNE.pdf", plot = p, width = 6, height = 6)

# [1.2] 按 Group 上色 (查看组间分布差异)
p <- ggplot(final_results$Data, aes(tSNE1, tSNE2, color = Group)) +
  geom_point(size = 0.1, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Group) + # 分面展示更清晰
  ggtitle("t-SNE by Group (Split View)") +
  guides(color = guide_legend(override.aes = list(size = 4)))

p

ggsave("Group_tSNE.pdf", plot = p, width = 10, height = 6)

# [1.3] 按样本 (Sample) 上色 (检查批次效应)
p <- ggplot(final_results$Data, aes(tSNE1, tSNE2, color = SampleID)) +
  geom_point(size = 0.1, alpha = 0.5) +
  theme_bw() +
  ggtitle("t-SNE by Sample")

p

ggsave("Sample_tSNE.pdf", plot = p, width = 7, height = 6)

# [1.4] 按基因/Marker 上色 (FeaturePlot)
# 请将 'PE-Cy7-A' 替换为你真实的 Marker 名字
target_marker <- "BV421-A" 
ggplot(final_results$Data, aes(tSNE1, tSNE2, color = .data[[target_marker]])) +
  geom_point(size = 0.1) +
  scale_color_viridis_c(option = "plasma") + # 热图色
  theme_bw() +
  ggtitle(paste0("FeaturePlot: ", target_marker))

#如果想换颜色的话
ggplot(final_results$Data, aes(tSNE1, tSNE2, color = .data[[target_marker]])) +
  geom_point(size = 0.1) +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
theme_bw() +
  theme(panel.grid = element_blank()) + 
  ggtitle(paste0("FeaturePlot: ", target_marker))

# [1.5] 按 CellType 上色 (仅适用于多群模式 B)
if("CellType" %in% colnames(final_results$Data)) {
  ggplot(final_results$Data, aes(tSNE1, tSNE2, color = CellType)) +
    geom_point(size = 0.1) +
    theme_bw() +
    ggtitle("t-SNE by User-Defined CellType")
}

# [1.6] 双阳图 (Dual Color)
# 红 + 蓝 = 紫 (双阳)
p <- plot_dual_color(final_results, "BV421-A", "APC-A")
p
ggsave("TNFa_IFNg_tSNE.pdf", plot = p, width = 6, height = 6)
# ===========================================================
# 📊 2. 统计图系列
# ===========================================================

# [2.1] 小提琴图：Cluster vs Marker (看亚群定义)
ggplot(final_results$Data, aes(x = Cluster, y = `APC-A`, fill = Cluster)) +
  geom_violin(scale = "width", trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none") +
  ggtitle("Expression by Cluster")

ggplot(final_results$Data, aes(x = Group, y = `PE-Cy5-A`, fill = Group)) + 
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  ggtitle("PE Expression by Group") + 
  stat_compare_means()
# [2.2] 箱线图：Group vs Marker (按 Cluster 分面，看组间差异)
ggplot(final_results$Data, aes(x = Group, y = `PE-Cy5-A`, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.1, alpha = 0.2) + # 显示散点
  theme_bw() +
  facet_wrap(~Cluster, scales = "free_y") + # 每个Cluster单独一张图
  stat_compare_means(label = "p.signif") +  # 自动加星号 (*, **)
  ggtitle("Marker Diff between Groups (per Cluster)")

# [2.3] 柱状图：细胞比例构成 (Composition)
# 看不同组里，各个 Cluster 的比例变化
p <- ggplot(final_results$Data, aes(x = Group, fill = Cluster)) +
  geom_bar(position = "fill") + # fill = 百分比堆叠
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  scale_fill_manual(values = my_colors) +
  labs(y = "Percentage", title = "Cluster Composition per Group")
p
ggsave("cell_percentage_Cluster.pdf", plot = p, width = 6, height = 6)
# [2.4] 多群比较图 (仅适用于多群模式 B)
# 比较你自己命名的群体 (如 CD4 vs CD8)
if("CellType" %in% colnames(final_results$Data)) {
  Compare_Populations(final_results$Data, "BV421-A")
}

# [2.5] 比较不同 Group 中 PE-Cy7 的表达峰值
# 场景 A: 如果你按 Cluster 分组
n_clusters <- length(unique(final_results$Data$Cluster))
colors_for_cluster <- colorRampPalette(brewer.pal(12, "Set3"))(n_clusters)
# 场景 B: 如果你按 Group 分组 (比如 pc, pl)
n_groups <- length(unique(final_results$Data$Group))
colors_for_group <- colorRampPalette(brewer.pal(8, "Set1"))(n_groups)
# fill = NA 已经写在函数里了
Plot_Density_Line(
  final_results, 
  marker = "PE-Cy7-A", 
  group_by = "Cluster", 
  custom_colors = colors_for_cluster
)
# [2.6] 看 Cluster 分布 (X=PE-Cy7, Y=APC)
Plot_Scatterhist_Line(
  final_results, 
  x_marker = "PE-Cy7-A", 
  y_marker = "APC-A", 
  group_by = "Cluster", 
  custom_colors = colors_for_cluster
)

Plot_Scatterhist_Line(
  final_results, 
  x_marker = "PE-Cy7-A", 
  y_marker = "APC-A", 
  group_by = "Group", 
  custom_colors = colors_for_group
)

# ===========================================================
# 📊 3. 热图系列
# ===========================================================
# [3.1] 热图数据导出
# 运行函数并把结果存到 hm 对象里
hm <- Plot_Heatmap(final_results, group_by = "Cluster", scale = "row")
# 这里的 hm$Matrix 就是你要的【未标准化】的原始矩阵！
print(head(hm$Matrix)) 
write.csv(hm$Matrix, "My_Heatmap_Raw_Data.csv")

# [3.1] 基础热图绘制
Plot_Heatmap(
  final_results, 
  group_by = "Cluster", 
  scale = "none",         # <--- 1. 不标准化 (展示原始ArcSinh值)
  cluster_cols = FALSE,   # <--- 2. 列不聚类 (保持 1,2,3...顺序)
  cluster_rows = TRUE,    #      行(Marker)依然聚类，把相似Marker放一起
  show_border = TRUE      # <--- 4. 显示网格线
)

# [3.2] 换颜色？标准化？
my_purple_yellow <- colorRampPalette(c("blue", "white", "red"))(100)
p <- Plot_Heatmap(
  final_results, 
  group_by = "Cluster",
  scale = "row",
  cluster_cols = FALSE, 
  color_palette = my_purple_yellow, # <--- 1. 自定义颜色
  show_border = FALSE               # <--- 2. 去掉格子白线
)
p
ggsave("heatmap_Cluster.pdf", plot = p$Plot, width = 8, height = 6)

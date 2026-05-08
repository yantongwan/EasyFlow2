# ==============================================================================
# 1. Environment setup, dependency loading, and module import
# ==============================================================================
getwd()
setwd("/Volumes/Extreme/lly go colitis facs")

pkgs <- c("shiny", "plotly", "miniUI", "dplyr", "Rtsne", "ggplot2",
          "stringr", "ggpubr", "scales", "RColorBrewer", "igraph", "RANN",
          "ggridges", "pheatmap", "httr2", "jsonlite", "hexbin", "viridis",
          "flowCore", "uwot", "ranger", "BiocManager", "flowStats",
          "reshape2", "tidyr")

new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

if (!requireNamespace("FlowSOM",            quietly = TRUE)) BiocManager::install("FlowSOM")
if (!requireNamespace("slingshot",          quietly = TRUE)) BiocManager::install(c("slingshot", "SingleCellExperiment"))
if (!requireNamespace("ggalluvial",         quietly = TRUE)) install.packages("ggalluvial")
if (!requireNamespace("ggrepel",            quietly = TRUE)) install.packages("ggrepel")

invisible(lapply(c(pkgs, "FlowSOM"), library, character.only = TRUE))

# Automatically load all core modules from the R function directory.
function_dir <- "/Volumes/1 TB/lly go colitis facs/EasyFlow/R_20260506"
r_files <- list.files(function_dir, pattern = "\\.[Rr]$", full.names = TRUE)

message(">>> Loading EasyFlow2 core modules...")
for (file in sort(r_files)) {
  source(file)
  message(sprintf("  Loaded: %s", basename(file)))
}
message(">>> EasyFlow2 engine is ready.")

# ==============================================================================
# Configuration
# ==============================================================================
my_folder <- "/Volumes/Extreme/lly go colitis facs/20240619 colitis cytokine/mln"

my_panel <- c(
  "BV421-A"              = "Cd206",
  "eFluor450-A"          = "Cd11b",
  "PacBlue-A"            = "Cd80",
  "BV510-A"              = "Cd11c",
  "BV605-A"              = "Cd4",
  "BV650-A"              = "Cd8",
  "BV711-A"              = "Cd45",
  "FITC-A"               = "Cd45.1",
  "PerCP-Cy5.5-A"        = "Ly6g",
  "PE-A"                 = "Cd64",
  "PE-Dazzle594-A"       = "Ly6c",
  "PE-Cy5-A"             = "MHC-II",
  "PE-Cy7-A"             = "F4/80",
  "PE-Fire810-A"         = "Cd44",
  "AF647-A"              = "Siglecf",
  "APC-A"                = "Cd62l",
  "AF700-A"              = "MER-TK",
  "LiveDeadFixableBlue-A" = "L/D"
)

my_panel <- c(
  "[live/dead fv-ef780]-A"              = "L/D",
  "BV605-A"            = "CD4",
  "BV510-A"            = "CD103",
  "FITC-A"               = "GM-CSF",
  "PE-Cy7-A"           = "IL17a",
  "BV421-A"            = "TNFa",
  "PE-A"                 = "IL6",
  "AF700-A"            = "IL10",
  "APC-A"              = "IFNg",
  "eFluor450-A"        = "Foxp3"
)
# ==============================================================================
# Network proxy and API configuration
# ==============================================================================
if (!nzchar(Sys.getenv("OPENROUTER_API_KEY"))) message("Set OPENROUTER_API_KEY before running LLM annotation.")
Sys.setenv(http_proxy  = "http://127.0.0.1:7890")
Sys.setenv(https_proxy = "http://127.0.0.1:7890")

tryCatch({
  res <- jsonlite::fromJSON("https://ipapi.co/json/")
  cat("Current R-session IP:", res$ip, "\n")
  cat("Current R-session location:", res$city, ",", res$country_name, "\n")
}, error = function(e) cat("Network request failed. Check the proxy or internet connection.\n"))

# ==============================================================================
# Step 1: Data import, QC, and preprocessing
# ==============================================================================
raw_data <- Prepare_Raw_Data(
  work_dir = my_folder, 
  group_keywords = c("wt", "ko"), 
  n_cells_per_sample = 50000, 
  do_QC = TRUE,                 # Enable debris and doublet removal
  do_compensation = TRUE,       # Apply the compensation matrix embedded in FCS files
  transform_method = "arcsinh", # Choose "arcsinh" or "logicle"
  cofactor = 150                # Typical values: 150 for flow cytometry, 5 for CyTOF
)

raw_data <- Apply_Panel_Mapping(raw_data, my_panel)
head(raw_data)
# ==============================================================================
# Step 2: Interactive gating
# ==============================================================================
gated_cells <- Run_Universal_Gating(raw_data)
Plot_Gating_Tree()
head(raw_data)
# ==============================================================================
# Step 3: Dimensionality reduction and clustering
# ==============================================================================
ref_results <- Run_Analysis_On_Gated_Data(
  gated_data       = gated_cells,
  k_neighbors      = 30,
  cluster_mode     = "B",
  do_scale         = TRUE,
  cluster_method   = "louvain",
  flowsom_k        = 20,
  tsne_perplexity  = NULL,   # NULL = infer automatically
  umap_n_neighbors = 15,     # UMAP parameter
  umap_min_dist    = 0.2     # UMAP parameter
)

# ==============================================================================
# Step 4: Dual-track AI annotation
# ==============================================================================
custom_rules <- list(
  "CD4 T cell"        = c("Cd4"),
  "Treg"              = c("Cd4", "Foxp3"),
  "Activated"         = c("Cd44"),
  "Naive"             = c("Cd62l"),
  "Memory"            = c("Cd44"),
  "IFNg+"             = c("Ifng"),
  "TNFa+"             = c("Tnfa")
)


my_ontology <- c(
  "CD4 T cell",
  "Unknown"
)

Normalize_Annotation_Label <- Normalize_Annotation_Label_Enhanced

final_results <- Annotate_Clusters_Dual_Track(
  res_obj                 = ref_results,
  rules                   = custom_rules,
  model                   = "openai/gpt-5.4",
  use_llm_only            = TRUE,
  ontology                = my_ontology,
  max_clusters_per_batch  = 20,
  base_population_context = paste(
    "All cells in this analysis were pre-gated as CD4+ T cells",
    "from mLN.",
    "All annotations must stay within the CD4 T-cell universe."
  ),
  allow_functional_prefix = TRUE,
  unknown_threshold       = 35,
  annotation_granularity  = "fine"
)

# If the AI annotations need manual curation, edit the report and rerun this block:
corrected_report <- read.csv("Cluster_Annotation_Report.csv")
final_results$Data$Cluster_Name <- corrected_report$fine_name[
match(final_results$Data$Cluster, corrected_report$cluster)]

# ==============================================================================
# Step 4b (optional): Large-scale machine-learning projection
# ==============================================================================
batch1_folder <- "/Volumes/1 TB/lly go colitis facs/20240619 colitis cytokine/mln"
batch2_folder <- "/Volumes/1 TB/lly go colitis facs/situmulate mlns"

batch1_label <- "Batch1_20240619_mLN"
batch2_label <- "Batch2_stimulated_mLN"
gating_template_path <- "/Volumes/Extreme_Pro/Codex/Gating_Thresholds_Template.csv"

batch1_raw_data <- Prepare_Raw_Data(
  work_dir         = batch1_folder,
  group_keywords   = c("wt", "ko"),
  n_cells_per_sample = 1500000,
  do_QC            = TRUE,
  do_compensation  = TRUE,
  transform_method = "arcsinh",
  cofactor         = 150
)
batch1_raw_data <- Apply_Panel_Mapping(batch1_raw_data, my_panel)
batch1_raw_data$Batch <- batch1_label

batch2_raw_data <- Prepare_Raw_Data(
  work_dir         = batch2_folder,
  group_keywords   = c("wt", "ko"),
  n_cells_per_sample = 1500000,
  do_QC            = TRUE,
  do_compensation  = TRUE,
  transform_method = "arcsinh",
  cofactor         = 150
)
batch2_raw_data <- Apply_Panel_Mapping(batch2_raw_data, my_panel)
batch2_raw_data$Batch <- batch2_label

# Apply the gating template. The margin parameter tightens boundaries and reduces outlier carryover.
batch1_raw_data <- Apply_Gating_Template(
  df            = batch1_raw_data,
  target_node   = "cell_live_cell_CD4",
  template_path = gating_template_path,
  margin        = 0.02
)

batch2_raw_data <- Apply_Gating_Template(
  df            = batch2_raw_data,
  target_node   = "cell_live_cell_CD4",
  template_path = gating_template_path,
  margin        = 0.02
)

massive_raw_data <- dplyr::bind_rows(batch1_raw_data, batch2_raw_data)
massive_raw_data$Batch <- factor(massive_raw_data$Batch, levels = c(batch1_label, batch2_label))

final_massive_results <- Project_Atlas_To_Massive_Data(                                                                                                               
  ref_obj              = final_results,                                                                                                                               
  massive_data         = massive_raw_data,                                                                                                                            
  target_label         = "Cluster_Name",                                                                                                                              
  cv_folds             = 5,                                                                                                                                           
  rejection_threshold  = 0.55,                                                                                                                                        
  use_class_weights    = TRUE,   # Enable class weighting                                                                                                               
  rare_cell_boost      = 2.0,    # Rare-cell protection multiplier
  batch_correction     = "adtnorm_lite",
  batch_col            = "Batch",
  reference_batch      = batch1_label,
  batch_diagnostic_pdf = "ADTnorm_Batch_Correction_Diagnostics.pdf"
)   

final_massive_results <- Set_Plot_Order(
  final_massive_results,
  order_list = list(
    Group = c("wt", "ko"),
    Batch = c(batch1_label, batch2_label)
  )
)

projection_dim <- "UMAP"
projection_color <- "Cluster_Name"
projection_split <- "Batch"

p_dim <- Plot_Dim_Reduction(final_massive_results, dim_method = projection_dim,
                            color_by = projection_color, split_by = projection_split)
out_dim_name <- paste0("01_", projection_dim, "_by_", projection_color,
                       "_split_", projection_split, "_ADTnormLite.pdf")
ggsave(out_dim_name, plot = p_dim, width = 12, height = 6)

# Optional rapid projection for a new sample. This reuses the trained model and skips CV.
if (exists("new_sample_data")) {
  new_sample_result <- Project_New_Sample(
    ref_obj             = final_massive_results,
    new_data            = new_sample_data,
    rejection_threshold = 0.55
  )
} else {
  message(">>> new_sample_data not found; skipping the optional rapid new-sample projection.")
}

# ==============================================================================
# Step 5: Central control panel
# ==============================================================================
head(final_results)
# [A] Dimensionality-reduction visualization
my_dim   <- "UMAP"
my_color <- "Cluster"
my_split <- NULL

# [B] Marker-expression comparison
my_target_marker <- "IL17a"
my_group_var     <- "Group"
my_fill_var      <- "Group"
my_facet_var     <- "Cluster_Name"

# [C] Cell-population composition
my_bar_x        <- "Group"
my_bar_fill     <- "Cluster_Name"
my_sankey_axis1 <- "Group"
my_sankey_axis2 <- "Cluster_Name"

# [D] Dual-marker co-expression
my_x_marker    <- "IL17a"
my_y_marker    <- "Foxp3"
my_density_grp <- "Cluster_Name"

# [E] Heatmap and dot plot
my_exp_grp     <- "Cluster_Name"
my_dot_markers <- NULL

# [F] Model evaluation
my_conf_true <- "CellType"
my_conf_pred <- "Cluster_Name"

# [G] Sample-level PCA
my_sample_group <- "Group"
my_sample_id    <- "SampleID"

final_results <- Set_Plot_Order(
  final_results,
  order_list = list(
    Group = c("wt", "ko")
  )
)

# ==============================================================================
# Step 6: Generate publication-ready figures
# ==============================================================================
message("\n>>> Generating publication-ready figures...")

# [6.1] Main dimensionality-reduction plot
p_dim <- Plot_Dim_Reduction(final_results, dim_method = my_dim,
                            color_by = my_color, split_by = my_split)
out_dim_name <- paste0("01_", my_dim, "_by_", my_color,
                       ifelse(is.null(my_split), "", paste0("_split_", my_split)), ".pdf")
ggsave(out_dim_name, plot = p_dim, width = 6, height = 6)

# [6.2] Cell-population composition plots (correct_sample_bias = TRUE reduces unequal sample-size bias)
p_comp_fill  <- Plot_Composition(final_results, x_var = my_bar_x,
                                 fill_var = my_bar_fill, position = "fill",
                                 correct_sample_bias = TRUE)
p_comp_stack <- Plot_Composition(final_results, x_var = my_bar_x,
                                 fill_var = my_bar_fill, position = "stack",
                                 correct_sample_bias = TRUE)
ggsave(paste0("02_Composition_Fill_by_",  my_bar_x, ".pdf"), plot = p_comp_fill,  width = 6, height = 6)
ggsave(paste0("02_Composition_Count_by_", my_bar_x, ".pdf"), plot = p_comp_stack, width = 6, height = 6)

# [6.3] Alluvial / Sankey flow plot
p_sankey <- Plot_Alluvial(final_results, axis1 = my_sankey_axis1,
                          axis2 = my_sankey_axis2, fill_var = my_bar_fill)
ggsave("03_Alluvial_Sankey_Flow.pdf", plot = p_sankey, width = 12, height = 6)

# [6.4] Dot plot
p_dot <- Plot_Dotplot(final_results, markers = my_dot_markers, group_by = my_exp_grp)
ggsave(paste0("04_Dotplot_Expression_by_", my_exp_grp, ".pdf"), plot = p_dot, width = 10, height = 6)

# [6.5] Violin, box, and bar plots
p_v <- Plot_Violin_AllInOne(
  res_obj        = final_results,
  marker         = my_target_marker,
  cluster_by     = my_facet_var,
  group_by       = my_group_var,
  fill_by        = my_fill_var,
  add_boxplot    = TRUE,
  trim           = FALSE,
  add_stats      = TRUE,
  stat_method    = "wilcox.test",   # or "t.test"
  stat_label     = "p.signif",      # or "p.format"
  hide_ns        = FALSE,
  p_y_offset     = 0.08,            # p-value label height
  logfc_y_offset = 0.16,            # logFC label height
  logfc_digits   = 2
)

ggsave(
  paste0("05_Violin_", my_target_marker, "_allinone_by_", my_facet_var, "_and_", my_group_var, ".pdf"),
  plot = p_v,
  width = 12,
  height = 5
)

p_bar <- Compare_Populations_Barplot(
  final_results,
  marker = my_target_marker,
  group_by = my_group_var,
  reference_group = "wt",
  show_error_bars = TRUE,
  logfc_threshold = 0.1,
  p_threshold = 0.05,
  facet_by = my_facet_var
)                                                                                                                                                                

if (!is.null(my_facet_var)) {                                                                                                                                         
  ggsave(paste0("06_Barplot_LogFC_", my_target_marker, "_by_", my_group_var, "_split_", my_facet_var, ".pdf"),
         plot = p_bar, width = 7, height = 10)                                                                                                                        
} else {                                                                                                                                                              
  ggsave(paste0("06_Barplot_LogFC_", my_target_marker, "_by_", my_group_var, ".pdf"),                                                                                 
         plot = p_bar, width = 8, height = 6)                                                                                                                         
}         


# [6.6] Expression heatmap
hm_res <- Plot_Heatmap(final_results, group_by = my_exp_grp, scale = "row",
                       cluster_rows = TRUE, cluster_cols = TRUE)
pdf(paste0("06_Heatmap_Expression_by_", my_exp_grp, ".pdf"), width = 8, height = 6)
print(hm_res$Plot)
invisible(dev.off())

# [6.7] Confusion matrix
if (my_conf_true %in% colnames(final_results$Data) &&
    my_conf_pred %in% colnames(final_results$Data)) {
  p_conf <- Plot_Confusion_Matrix(final_results, true_label = my_conf_true, pred_label = my_conf_pred)
  ggsave("07_Confusion_Matrix_Evaluation.pdf", plot = p_conf, width = 8, height = 7)
}

# [6.8] Ridge plot
p_density <- Plot_Density_Line(final_results, marker = my_target_marker, group_by = my_group_var)
ggsave(paste0("08_Ridge_", my_target_marker, "_by_", my_group_var, ".pdf"),
       plot = p_density, width = 6, height = 5)

# [6.9] Dual marginal-density scatter plot
p_scatterhist <- Plot_Scatterhist_Line(
  final_results,
  x_marker = my_x_marker,
  y_marker = my_y_marker,
  group_by = my_density_grp
)

ggplot2::ggsave(
  filename = paste0("09_Scatterhist_", my_x_marker, "_vs_", my_y_marker, ".pdf"),
  plot = p_scatterhist,
  width = 7,
  height = 7
)

# [6.10] Sample-level PCA
p_sample_pca <- Plot_Sample_PCA(final_results, color_by = my_sample_group, label_by = my_sample_id)
ggsave(paste0("10_Sample_Level_PCA_by_", my_sample_group, ".pdf"),
       plot = p_sample_pca, width = 7, height = 6)

message("\n>>> All publication-ready figures have been generated.")



# ==============================================================================
# Step 7: Statistical analysis and report output
# ==============================================================================

my_stat_target <- "Cluster_Name"
my_stat_group  <- "Group"
my_stat_method <- "wilcoxon"
my_stat_trans  <- "clr"

# [7.1] Differential abundance statistics
stats_results <- Perform_Abundance_Stats(
  res_obj      = final_results,
  target_label = my_stat_target,
  group_col    = my_stat_group,
  method       = my_stat_method,
  transform    = my_stat_trans
)

# [7.2] Optional marker differential-expression analysis
# de_results <- Perform_Marker_DE(
#   res_obj     = final_results,
#   group_col   = my_stat_group,
#   cluster_col = "Cluster_Name",
#   method      = "wilcoxon"
# )

# [7.3] Trajectory inference
traj_res     <- Infer_Trajectory(res_obj = final_results, start_cluster = NULL, use_dim = "UMAP")
final_results <- traj_res$res_obj

# [7.4] Export FCS files
Export_To_FCS(final_results, out_dir = "09_Final_FCS_Output")

# [7.5] Generate HTML reports
Generate_Final_HTML_Report(
  res_obj      = final_results,
  stats_res    = stats_results,
  project_name = "EasyFlow2 Final Report",
  out_file     = "EasyFlow2_Final_Report.html",
  template     = "manuscript"
)

Generate_Final_HTML_Report(
  res_obj   = final_results,
  stats_res = stats_results,
  out_file  = "EasyFlow2_Review_Report.html",
  template  = "review"
)

Generate_Final_HTML_Report(
  res_obj       = final_results,
  stats_res     = stats_results,
  out_file      = "EasyFlow2_Custom_Report.html",
  include_plots = c("sample_counts", "embedding", "composition",
                    "marker_dotplot", "trajectory", "abundance")
)

# ==============================================================================
# Module 3: Core Dimensionality Reduction and Clustering Engine (Scaling + Leiden/FlowSOM + Benchmark)
# Filename: R/03_DimReduction.R
# ==============================================================================

#' Core dimensionality reduction and clustering engine
#' @param gated_data Input cell data frame
#' @param k_neighbors Number of KNN neighbors (used for Louvain/Leiden graph construction)
#' @param cluster_mode "A" = "A" = cluster labeled cells only; "B" = cluster all gated-in cells (including Ungated)
#' @param do_scale Whether to perform Z-score scaling before clustering/dimensionality reduction
#' @param cluster_method Select clustering algorithm: "louvain", "leiden", "flowsom"
#' @param flowsom_k If FlowSOM is selected, the expected number of meta-clusters (default 20)
#' @param tsne_perplexity t-SNE Perplexity parameter; if NULL, automatically inferred from cell count
#' @param umap_n_neighbors UMAP Number of neighbors (default 15)
#' @param umap_min_dist UMAP Minimum-distance parameter controlling embedding compactness (default 0.2)
Run_Analysis_On_Gated_Data <- function(gated_data,
                                       k_neighbors = 30,
                                       cluster_mode = "B",
                                       do_scale = FALSE,
                                       cluster_method = c("louvain", "leiden", "flowsom"),
                                       flowsom_k = 20,
                                       tsne_perplexity = NULL,
                                       umap_n_neighbors = 15,
                                       umap_min_dist = 0.2) {
  
  cluster_method <- match.arg(cluster_method)
  benchmarks <- list()
  total_start <- Sys.time()
  
  message("\n=======================================================")
  message(sprintf(">>> [Analysis engine start] mode: %s | algorithm: %s | Scaling: %s", cluster_mode, toupper(cluster_method), do_scale))
  
  # --- 0. Data-mode filtering ---
  if (cluster_mode == "A") {
    message(">>> [Mode A] Perform in-depth clustering analysis only on specifically named subpopulations...")
    analysis_data <- gated_data[!gated_data$CellType %in% c("Ungated", "Unknown"), ]
  } else {
    message(">>> [Mode B] Perform global clustering on all filtered cells (including Ungated)...")
    analysis_data <- gated_data
  }
  if (nrow(analysis_data) < 50) stop("  Too few eligible cells for dimensionality reduction and clustering.")
  
  # Extract pure feature matrix
  ignore_cols <- c("tSNE1", "tSNE2", "UMAP1", "UMAP2", "Cluster", "SampleID", "Group",
                   "Time", "TIME", "CellType", "Pseudotime", "Gate_Status")
  numeric_cols <- colnames(analysis_data)[sapply(analysis_data, is.numeric)]
  marker_cols <- setdiff(numeric_cols, ignore_cols)
  cluster_mat <- analysis_data[, marker_cols]
  
  # --- 1. Scaling (Z-score Scaling) ---
  step_start <- Sys.time()
  if (do_scale) {
    message(">>> Performing feature Z-score scaling...")
    cluster_mat <- as.data.frame(scale(cluster_mat))
  }
  benchmarks[["Scaling"]] <- round(as.numeric(difftime(Sys.time(), step_start, units = "secs")), 2)
  
  # --- 2. Clustering ---
  step_start <- Sys.time()
  message(sprintf(">>> Running clustering algorithm [%s]...", toupper(cluster_method)))
  
  if (cluster_method %in% c("louvain", "leiden")) {
    message(paste0("    Building KNN topology graph (k = ", k_neighbors, ")..."))
    knn_res <- RANN::nn2(data = cluster_mat, k = k_neighbors)
    edges <- matrix(c(rep(1:nrow(cluster_mat), each = k_neighbors), as.vector(t(knn_res$nn.idx))), ncol = 2)
    g <- igraph::graph_from_edgelist(edges, directed = FALSE)
    
    set.seed(123)
    if (cluster_method == "louvain") {
      analysis_data$Cluster <- as.factor(igraph::cluster_louvain(g)$membership)
    } else {
      analysis_data$Cluster <- as.factor(igraph::cluster_leiden(g, objective_function = "modularity")$membership)
    }
    
  } else if (cluster_method == "flowsom") {
    if (!requireNamespace("FlowSOM", quietly = TRUE)) stop("  Missing the 'FlowSOM' package; please run BiocManager::install('FlowSOM')")
    
    #   Change 1: dynamically adjust SOM grid size according to cell count to avoid overfitting on small datasets
    n_cells <- nrow(cluster_mat)
    # Target number of micro-cluster nodes is approximately min(cell_count/50, 400); take the square root to obtain the grid side length
    som_nodes <- max(25, min(400, floor(n_cells / 50)))
    som_dim   <- ceiling(sqrt(som_nodes))
    message(sprintf("    Training FlowSOM self-organizing map (grid: %dx%d = %d nodes, target meta-clusters = %d)...",
                    som_dim, som_dim, som_dim^2, flowsom_k))
    
    set.seed(42)
    fsom_input    <- FlowSOM::ReadInput(as.matrix(cluster_mat))
    fsom_model    <- FlowSOM::BuildSOM(fsom_input, xdim = som_dim, ydim = som_dim, silent = TRUE)
    fsom_model    <- FlowSOM::BuildMST(fsom_model, silent = TRUE)
    meta_clustering <- FlowSOM::metaClustering_consensus(fsom_model$map$codes, k = flowsom_k)
    analysis_data$Cluster <- as.factor(meta_clustering[fsom_model$map$mapping[, 1]])
  }
  benchmarks[["Clustering"]] <- round(as.numeric(difftime(Sys.time(), step_start, units = "secs")), 2)
  
  # --- 3. tSNE dimensionality reduction ---
  step_start <- Sys.time()
  message(">>> Running t-SNE dimensionality reduction...")
  set.seed(42)
  
  #   Change 2: support user-defined perplexity; if NULL, infer automatically
  if (is.null(tsne_perplexity)) {
    tsne_perplexity <- ifelse(nrow(analysis_data) < 91, floor((nrow(analysis_data) - 1) / 3), 30)
  } else {
    # Safety check: perplexity cannot exceed (n-1)/3
    max_perp <- floor((nrow(analysis_data) - 1) / 3)
    if (tsne_perplexity > max_perp) {
      warning(sprintf("tsne_perplexity=%d is outside the safe range and has been automatically adjusted to %d", tsne_perplexity, max_perp))
      tsne_perplexity <- max_perp
    }
  }
  message(sprintf("    t-SNE perplexity = %d", tsne_perplexity))
  
  tsne_out <- Rtsne::Rtsne(cluster_mat, dims = 2, perplexity = tsne_perplexity,
                           verbose = FALSE, check_duplicates = FALSE)
  analysis_data$tSNE1 <- tsne_out$Y[, 1]
  analysis_data$tSNE2 <- tsne_out$Y[, 2]
  benchmarks[["tSNE"]] <- round(as.numeric(difftime(Sys.time(), step_start, units = "secs")), 2)
  
  # --- 4. UMAP dimensionality reduction and manifold saving ---
  step_start <- Sys.time()
  message(sprintf(">>> Running UMAP (n_neighbors=%d, min_dist=%.2f)...", umap_n_neighbors, umap_min_dist))
  set.seed(42)
  umap_model <- NULL
  if (requireNamespace("uwot", quietly = TRUE)) {
    #   Change 3: use user-supplied n_neighbors and min_dist parameters
    umap_model <- uwot::umap(cluster_mat,
                             n_neighbors = umap_n_neighbors,
                             min_dist    = umap_min_dist,
                             ret_model   = TRUE,
                             verbose     = FALSE)
    analysis_data$UMAP1 <- umap_model$embedding[, 1]
    analysis_data$UMAP2 <- umap_model$embedding[, 2]
  }
  benchmarks[["UMAP"]] <- round(as.numeric(difftime(Sys.time(), step_start, units = "secs")), 2)
  
  # --- 5. performance benchmark report ---
  benchmarks[["Total"]] <- round(as.numeric(difftime(Sys.time(), total_start, units = "secs")), 2)
  
  message("\n=======================================================")
  message("[Benchmark] Module performance benchmark report (elapsed time):")
  if (do_scale) message(sprintf("   Scaling (Scaling):  %.2f  sec", benchmarks[["Scaling"]]))
  message(sprintf("   %s clustering:            %.2f sec", toupper(cluster_method), benchmarks[["Clustering"]]))
  message(sprintf("   t-SNE dimensionality reduction:         %.2f  sec", benchmarks[["tSNE"]]))
  message(sprintf("   UMAP dimensionality reduction:          %.2f  sec", benchmarks[["UMAP"]]))
  message("   --------------------------------------")
  message(sprintf("   Total engine time:         %.2f  sec", benchmarks[["Total"]]))
  message("=======================================================\n")
  
  return(list(
    Data = analysis_data,
    Info = list(
      Mode           = cluster_mode,
      Scaled         = do_scale,
      Cluster_Method = cluster_method,
      UMAP_Model     = umap_model,
      UMAP_Params    = list(n_neighbors = umap_n_neighbors, min_dist = umap_min_dist),
      tSNE_Params    = list(perplexity = tsne_perplexity),
      Benchmarks     = benchmarks
    )
  ))
}

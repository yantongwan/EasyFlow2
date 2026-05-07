---
title: EasyFlow2
---

# EasyFlow2

EasyFlow2 is an R package for streamlined flow cytometry analysis, from FCS preprocessing and interactive gating to clustering, annotation, visualization, statistics, trajectory inference, FCS export, and HTML reporting.

## Install

```r
install.packages(c("remotes", "BiocManager"))

BiocManager::install(c(
  "flowCore",
  "flowStats",
  "FlowSOM",
  "Biobase",
  "S4Vectors",
  "SingleCellExperiment",
  "slingshot"
))

install.packages(c("broom.mixed", "lmerTest"))  # optional, for method = "lmm"

remotes::install_github(
  "yantongwan/EasyFlow2",
  dependencies = TRUE,
  build_vignettes = FALSE
)
```

```r
library(EasyFlow2)
```

## Quick Start

```r
raw_data <- Prepare_Raw_Data(
  work_dir = "path/to/fcs_files",
  group_keywords = c("Control", "Treatment"),
  n_cells_per_sample = 5000,
  do_QC = TRUE,
  do_compensation = TRUE,
  transform_method = "arcsinh"
)

gated_data <- Run_Universal_Gating(raw_data)

res <- Run_Analysis_On_Gated_Data(
  gated_data,
  cluster_mode = "B",
  cluster_method = "louvain"
)

Plot_Dim_Reduction(res, dim_method = "UMAP", color_by = "Cluster")
```

## Included Workflow Manual

```r
system.file(
  "workflows",
  "EasyFlow2_Executable_Analysis_Workflow_Integrated.Rmd",
  package = "EasyFlow2"
)
```

The repository also includes the rendered integrated PDF manual under `inst/workflows/`.

## Core Capabilities

- FCS import, QC, compensation, transformation, downsampling, and panel mapping.
- Interactive Shiny/Plotly gating with gating-tree tracking.
- Louvain, Leiden, and FlowSOM clustering with t-SNE and UMAP.
- Rule-based and LLM-assisted cluster annotation.
- Atlas projection and new-sample projection.
- Publication-ready flow cytometry visualizations.
- Differential abundance, marker DE, trajectory inference, FCS export, and automated reports.

# EasyFlow2

EasyFlow2 is an R package for end-to-end flow cytometry analysis. It integrates FCS preprocessing, compensation, transformation, interactive gating, gating-tree export, dimensionality reduction, clustering, annotation, atlas projection, visualization, abundance statistics, trajectory inference, FCS export, and automated HTML reporting.

## Installation

Install the package directly from GitHub:

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

Then load it:

```r
library(EasyFlow2)
```

If you only need the core package first, install GitHub dependencies with:

```r
remotes::install_github("yantongwan/EasyFlow2", build_vignettes = FALSE)
```

Install the optional Bioconductor packages above before using FlowSOM clustering, FCS export, or trajectory inference.

## What EasyFlow2 Does

- Reads FCS files and performs quality control, compensation, transformation, downsampling, and panel-name mapping.
- Provides an interactive Shiny and Plotly gating workstation with gating-tree tracking and threshold-template export.
- Runs graph-based clustering, FlowSOM clustering, t-SNE, and UMAP.
- Supports rule-based and LLM-assisted cell-type annotation.
- Projects manually curated reference atlases onto large-scale or new flow cytometry datasets.
- Produces publication-ready UMAP/t-SNE, composition, alluvial, dotplot, violin, ridge, heatmap, PCA, and scatter-density plots.
- Performs differential abundance analysis, marker differential expression, trajectory inference, FCS export, and final HTML report generation.

## Basic Workflow

```r
library(EasyFlow2)

raw_data <- Prepare_Raw_Data(
  work_dir = "path/to/fcs_files",
  group_keywords = c("Control", "Treatment"),
  n_cells_per_sample = 5000,
  do_QC = TRUE,
  do_compensation = TRUE,
  transform_method = "arcsinh",
  cofactor = 150
)

panel <- c(
  "FITC-A" = "CD4",
  "PE-A" = "CD8",
  "APC-A" = "CD45"
)
raw_data <- Apply_Panel_Mapping(raw_data, panel)

gated_data <- Run_Universal_Gating(raw_data)

res <- Run_Analysis_On_Gated_Data(
  gated_data,
  cluster_mode = "B",
  cluster_method = "louvain",
  do_scale = FALSE
)

Plot_Dim_Reduction(res, dim_method = "UMAP", color_by = "Cluster")
Plot_Composition(res, x_var = "Group", fill_var = "Cluster")
```

## Integrated Workflow Manual

The integrated executable workflow manual is installed with the package:

```r
system.file(
  "workflows",
  "EasyFlow2_Executable_Analysis_Workflow_Integrated.Rmd",
  package = "EasyFlow2"
)

system.file(
  "workflows",
  "EasyFlow2_Executable_Analysis_Workflow_Integrated.pdf",
  package = "EasyFlow2"
)
```

Template files are also available:

```r
system.file("templates", "Gating_Thresholds_Template.csv", package = "EasyFlow2")
system.file("templates", "Gating_Tree_Structure.csv", package = "EasyFlow2")
```

## Main Function Modules

| Module | Role |
|---|---|
| `01_Preprocessing.R` | FCS import, QC, compensation, transformation, downsampling, panel mapping, template gating |
| `02_Interactive_UI.R` | interactive gating, terminal-gate assignment, gating-tree visualization |
| `03_DimReduction.R` | clustering, t-SNE, UMAP, benchmarking |
| `04_Annotation.R` | rule-based and LLM-assisted annotation |
| `05_MachineLearning.R` | ADTnorm-lite correction and atlas/new-sample projection |
| `06_Visualization.R` | publication-ready plots |
| `07_Stats_and_Export.R` | abundance statistics, marker DE, trajectory inference, FCS export, HTML report |

## Notes

EasyFlow2 is designed for research workflows with potentially large single-cell flow cytometry outputs. Large intermediate `.rds`, FCS export, figure, and report-output folders should be regenerated locally and are intentionally excluded from the package source.

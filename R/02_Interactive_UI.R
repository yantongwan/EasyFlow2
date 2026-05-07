# ==============================================================================
# Module 2: Interactive Gating System (intelligent trimming + explicit gating-tree topology tracking/export + density view)
# Filename: R/02_Interactive_UI.R
# ==============================================================================

.Infer_Terminal_Gate_Nodes <- function(gate_labels, tree_df = NULL) {
  gate_labels <- unique(as.character(gate_labels))
  gate_labels <- gate_labels[!is.na(gate_labels) & nzchar(gate_labels)]
  gate_labels <- setdiff(gate_labels, c("Ungated", "Unknown"))
  
  if (!is.null(tree_df) &&
      all(c("Node", "Parent") %in% colnames(tree_df)) &&
      nrow(tree_df) > 0) {
    nodes <- unique(as.character(tree_df$Node))
    parents <- unique(as.character(tree_df$Parent))
    terminal_nodes <- setdiff(nodes, parents)
    terminal_nodes <- terminal_nodes[!is.na(terminal_nodes) & nzchar(terminal_nodes)]
    return(sort(terminal_nodes))
  }
  
  # Fallback for cached objects without a tree CSV: infer internal nodes from
  # EasyFlow's parent_child naming convention.
  internal_nodes <- gate_labels[vapply(gate_labels, function(label) {
    any(startsWith(gate_labels, paste0(label, "_")))
  }, logical(1))]
  
  sort(setdiff(gate_labels, internal_nodes))
}

#' Assign final CellType labels from terminal gates only
#'
#' Cells assigned to an internal gate, Ungated, or no gate are labeled Unknown.
#' The original gate path is preserved in Gate_Path for traceability.
Assign_Terminal_Gate_CellType <- function(df,
                                          tree_df = NULL,
                                          tree_path = "Gating_Tree_Structure.csv",
                                          gate_col = NULL,
                                          output_col = "CellType",
                                          unknown_label = "Unknown",
                                          keep_gate_path = TRUE) {
  if (is.null(gate_col)) {
    if ("Gate_Status" %in% colnames(df)) {
      gate_col <- "Gate_Status"
    } else if (output_col %in% colnames(df)) {
      gate_col <- output_col
    } else {
      df[[output_col]] <- unknown_label
      if (keep_gate_path && !"Gate_Path" %in% colnames(df)) df$Gate_Path <- NA_character_
      attr(df, "Terminal_Gates") <- character(0)
      message(sprintf(
        ">>> No gate label column was found; all %d cells were labeled %s.",
        nrow(df),
        unknown_label
      ))
      return(df)
    }
  }
  
  if (!gate_col %in% colnames(df)) {
    stop("Gate label column not found: ", gate_col)
  }
  
  if (is.null(tree_df) && !is.null(tree_path) && file.exists(tree_path)) {
    tree_df <- utils::read.csv(tree_path, stringsAsFactors = FALSE)
  }
  
  gate_labels <- as.character(df[[gate_col]])
  terminal_nodes <- .Infer_Terminal_Gate_Nodes(gate_labels, tree_df = tree_df)
  
  if (keep_gate_path && !"Gate_Path" %in% colnames(df)) {
    df$Gate_Path <- gate_labels
  }
  
  df[[output_col]] <- ifelse(
    !is.na(gate_labels) & gate_labels %in% terminal_nodes,
    gate_labels,
    unknown_label
  )
  
  message(sprintf(
    ">>> Terminal-gate CellType assignment completed: %d terminal gates; %d cells labeled %s.",
    length(terminal_nodes),
    sum(df[[output_col]] == unknown_label, na.rm = TRUE),
    unknown_label
  ))
  
  attr(df, "Terminal_Gates") <- terminal_nodes
  df
}

Run_Universal_Gating <- function(merged_df) {
  
  numeric_cols <- colnames(merged_df)[sapply(merged_df, is.numeric)]
  numeric_cols <- numeric_cols[!grepl("tSNE|UMAP|Cluster", numeric_cols)]
  
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("EasyFlow2 Gating Tree Workstation)"),
    miniUI::miniContentPanel(padding = 0,
                             # Increase height to accommodate the new checkbox
                             shiny::div(style = "padding: 15px; height: 245px; background-color: #fcfcfc; border-bottom: 1px solid #ddd;",
                                        shiny::fluidRow(
                                          shiny::column(3, 
                                                        shiny::selectInput("x_axis", "X Axis:", choices = numeric_cols, selected = numeric_cols[1]),
                                                        #   Change 1: add a checkbox to control density display
                                                        shiny::checkboxInput("show_density", "Show density", value = TRUE)
                                          ),
                                          shiny::column(3, shiny::selectInput("y_axis", "Y Axis:", choices = numeric_cols, selected = numeric_cols[2])),
                                          shiny::column(3, shiny::selectInput("parent_pop", "  Current level:", choices = c("All_Cells"), selected = "All_Cells")),
                                          shiny::column(3, shiny::div(
                                            shiny::h5(shiny::textOutput("cell_count_text"), style = "color:#2c3e50; font-weight:bold; margin-top: 20px; text-align:right;"),
                                            shiny::h6(shiny::textOutput("history_step_text"), style = "color:#7f8c8d; text-align:right; margin-top: 2px;")
                                          ))
                                        ),
                                        shiny::fluidRow(
                                          shiny::column(4,
                                                        shiny::div(style = "background-color: #fdf2e9; padding: 10px; border-radius: 5px; border: 1px solid #e67e22; height: 125px;",
                                                                   shiny::strong("1. Data cleaning (global)", style = "color:#d35400;"), shiny::br(),
                                                                   shiny::div(
                                                                     style = "margin-top: 5px; display:flex; flex-wrap:wrap; gap:8px;",
                                                                     shiny::actionButton("keep_btn", "   Set as global gate", class = "btn-primary", style = "padding: 4px 8px;"),
                                                                     shiny::actionButton("remove_btn", "   Remove selection", class = "btn-danger", style = "padding: 4px 8px;"),
                                                                     shiny::actionButton("undo_btn", "   Undo action", class = "btn-warning", style = "padding: 4px 8px;")
                                                                   )
                                                        )
                                          ),
                                          shiny::column(4,
                                                        shiny::div(style = "background-color: #e8f8f5; padding: 10px; border-radius: 5px; border: 1px solid #1abc9c; height: 125px;",
                                                                   shiny::strong("2. Subpopulation extraction (inherit parent path)", style = "color:#16a085;"), shiny::br(),
                                                                   shiny::div(style = "display: flex; gap: 10px; margin-top: 5px;",
                                                                              shiny::textInput("group_name", NULL, value = "CD45", width = "120px"),
                                                                              shiny::actionButton("mark_btn", "  Annotate selection", class = "btn-success", style = "height: 34px;")
                                                                   )
                                                        )
                                          ),
                                          shiny::column(4,
                                                        shiny::div(style = "background-color: #f4f6f7; padding: 8px; border-radius: 5px; border: 1px solid #bdc3c7; height: 125px;",
                                                                   shiny::div(style = "display:flex; justify-content:space-between; align-items:center;",
                                                                              shiny::strong("  Tree branches", style = "color:#2c3e50;"),
                                                                              shiny::actionButton("export_btn", "  Export template", class = "btn-default btn-xs", style = "padding: 2px 5px;")
                                                                   ),
                                                                   shiny::div(style = "height: 50px; overflow-y: auto; margin-top: 5px; border-top: 1px dashed #ccc; padding-top: 5px;",
                                                                              shiny::uiOutput("legend_ui")
                                                                   )
                                                        )
                                          )
                                        )
                             ),
                             # Adjust the plot area height accordingly
                             shiny::div(style = "height: calc(100% - 260px); padding: 10px;",
                                        plotly::plotlyOutput("scatter_plot", height = "100%")
                             )
    )
  )
  
  server <- function(input, output, session) {
    init_df <- merged_df
    if (!"Gate_Status" %in% colnames(init_df)) init_df$Gate_Status <- "Ungated"
    
    values <- shiny::reactiveValues(
      data = init_df,
      #   Optimization 1: save only the Gate_Status column in history (with row names), reducing size by 95%
      history = list(init_df[, "Gate_Status", drop = FALSE]),
      tree = data.frame(Node = character(), Parent = character(), stringsAsFactors = FALSE)
    )
    
    my_palette <- c("#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#17BECF", "#BCBD22", "#1F77B4", "#8C564B", "black")
    
    shiny::observe({
      gates <- unique(values$data$Gate_Status)
      valid_gates <- sort(setdiff(gates, "Ungated"))
      new_choices <- c("All_Cells", valid_gates)
      curr_sel <- if (input$parent_pop %in% new_choices) input$parent_pop else "All_Cells"
      shiny::updateSelectInput(session, "parent_pop", choices = new_choices, selected = curr_sel)
    })
    
    output$cell_count_text <- shiny::renderText({
      if (input$parent_pop == "All_Cells") {
        paste0("Global total cells: ", nrow(values$data), " / initial: ", nrow(merged_df))
      } else {
        pop_count <- sum(values$data$Gate_Status == input$parent_pop | startsWith(values$data$Gate_Status, paste0(input$parent_pop, "_")))
        paste0("[", input$parent_pop, "] node cell count: ", pop_count)
      }
    })
    
    output$history_step_text <- shiny::renderText({
      current_step <- length(values$history)
      max_steps <- 15
      paste0("History: ", current_step, "/", max_steps)
    })
    
    output$legend_ui <- shiny::renderUI({
      gated_idx <- values$data$Gate_Status != "Ungated"
      if (!any(gated_idx)) return(shiny::p("No annotated populations yet", style = "color: #7f8c8d; font-size: 12px;"))
      
      unique_gates <- sort(unique(values$data$Gate_Status[gated_idx]))
      color_map <- setNames(rep(my_palette, length.out = length(unique_gates)), unique_gates)
      
      tags_list <- lapply(unique_gates, function(g) {
        #   Change: calculate cumulative cell counts including this node and all child nodes
        # Child nodes start with "parent_"; use this property for matching
        count <- sum(values$data$Gate_Status == g | startsWith(values$data$Gate_Status, paste0(g, "_")))
        
        shiny::tags$span(
          style = paste0("background-color:", color_map[g], "; color: white; padding: 2px 6px; border-radius: 8px; margin: 2px; display: inline-block; font-size: 11px; box-shadow: 1px 1px 2px rgba(0,0,0,0.2);"),
          paste0(g, " (", count, ")")
        )
      })
      shiny::tagList(tags_list)
    })
    
    output$scatter_plot <- plotly::renderPlotly({
      shiny::req(input$x_axis, input$y_axis)
      
      # 1. Slice the data for the current gating level
      plot_data <- values$data
      if (input$parent_pop != "All_Cells") {
        target_idx <- plot_data$Gate_Status == input$parent_pop | startsWith(plot_data$Gate_Status, paste0(input$parent_pop, "_"))
        plot_data <- plot_data[target_idx, ]
      }
      if (nrow(plot_data) == 0) return(plotly::plot_ly() %>% plotly::layout(title = "No cells at this level"))
      
      x_name <- input$x_axis
      y_name <- input$y_axis
      x_val <- plot_data[[x_name]]
      y_val <- plot_data[[y_name]]
      
      # 2. Minimal DataFrame
      mini_plot_data <- data.frame(
        x = x_val,
        y = y_val,
        Gate_Status = as.character(plot_data$Gate_Status),
        Cell_ID = row.names(plot_data),
        stringsAsFactors = FALSE
      )
      
      # 3. Calculate density heatmap gradient colors (blue-cyan-yellow-red)
      point_colors <- if (length(x_val) > 2 && var(x_val) > 0) {
        grDevices::densCols(x_val, y_val, colramp = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))
      } else {
        rep("blue", length(x_val))
      }
      
      # Override with dedicated colors for already gated subpopulations
      gated_idx <- mini_plot_data$Gate_Status != "Ungated"
      if (any(gated_idx)) {
        global_unique <- sort(unique(values$data$Gate_Status[values$data$Gate_Status != "Ungated"]))
        color_map <- setNames(rep(my_palette, length.out = length(global_unique)), global_unique)
        point_colors[gated_idx] <- color_map[mini_plot_data$Gate_Status[gated_idx]]
      }
      
      # 4. Adaptive point size and opacity (preserving your larger point-size settings)
      pt_size <- ifelse(nrow(mini_plot_data) > 200000, 3, 4)
      pt_opac <- ifelse(nrow(mini_plot_data) > 200000, 0.08, 0.25) 
      
      # Initialize the canvas
      p <- plotly::plot_ly(source = "subset_plot")
      
      #   Draw the base-layer WebGL scatter first
      p <- p %>% plotly::add_trace(
        data = mini_plot_data,
        x = ~x, 
        y = ~y,
        key = ~Cell_ID,
        type = 'scattergl',
        mode = 'markers',
        marker = list(size = pt_size, color = point_colors, opacity = pt_opac),
        # Force opacity to remain fixed so surrounding cells do not disappear during selection
        selected = list(marker = list(opacity = pt_opac)),
        unselected = list(marker = list(opacity = pt_opac))
      )
      
      #   Key breakthrough: compute pure paths and overlay WebGL lines!
      if (isTRUE(input$show_density)) {
        if (nrow(mini_plot_data) <= 50000) {
          
          density_data <- mini_plot_data
          if (nrow(density_data) > 25000) {
            density_data <- density_data[sample(1:nrow(density_data), 25000), ]
          }
          
          # Step A: compute 2D kernel density in the R backend
          kd <- MASS::kde2d(density_data$x, density_data$y, n = 50)
          
          # Step B: extract contour coordinate matrices (nlevels=12 approximates the previously full-looking distribution)
          clines <- grDevices::contourLines(x = kd$x, y = kd$y, z = kd$z, nlevels = 12)
          
          if (length(clines) > 0) {
            # Step C: connect fragmented contour pieces into continuous polylines (using NA breaks; the key to WebGL rendering)
            line_x <- unlist(lapply(clines, function(cl) c(cl$x, NA)))
            line_y <- unlist(lapply(clines, function(cl) c(cl$y, NA)))
            
            # Step D: add the lines as scattergl too! Elements added later in the same engine stay on top!
            p <- p %>% plotly::add_trace(
              x = line_x, 
              y = line_y,
              type = 'scattergl',  # Convert everything to the accelerated WebGL engine
              mode = 'lines',
              line = list(color = 'black', width = 2),
              hoverinfo = "none",
              showlegend = FALSE
            )
          }
          
        } else {
          shiny::showNotification("   Cell count exceeds 50,000; contour overlays were automatically disabled to prevent lag.", type = "warning", duration = 4)
        }
      }
      
      # Finally set the overall layout
      p <- p %>% plotly::layout(
        xaxis = list(title = x_name),
        yaxis = list(title = y_name),
        dragmode = "lasso",
        showlegend = FALSE,
        hovermode = FALSE 
      )
      
      p
    })
    
    save_history <- function() {
      #   Optimization 1: save only the lightweight status column
      mini_state <- values$data[, "Gate_Status", drop = FALSE]
      values$history <- c(values$history, list(mini_state))
      if (length(values$history) > 15) values$history <- values$history[-1]
    }
    
    shiny::observeEvent(input$keep_btn, {
      ed <- plotly::event_data("plotly_selected", source = "subset_plot")
      if (is.null(ed)) return()
      save_history()
      values$data <- values$data[row.names(values$data) %in% as.character(ed$key), ]
      shiny::showNotification("  Set as the global gate!", type = "message")
    })
    
    shiny::observeEvent(input$remove_btn, {
      ed <- plotly::event_data("plotly_selected", source = "subset_plot")
      if (is.null(ed)) return()
      save_history()
      values$data <- values$data[!(row.names(values$data) %in% as.character(ed$key)), ]
      shiny::showNotification("   Selected debris cells were removed!", type = "message")
    })
    
    shiny::observeEvent(input$undo_btn, {
      if (length(values$history) > 1) {
        values$history <- values$history[-length(values$history)]
        last_state <- values$history[[length(values$history)]]
        
        #   Optimization 1: on undo, restore full data directly from init_df using row names
        values$data <- init_df[rownames(last_state), ]
        values$data$Gate_Status <- last_state$Gate_Status
        
        shiny::showNotification("   Undo successful", type = "message")
      } else {
        shiny::showNotification("   Already at the earliest history state", type = "warning")
      }
    })
    
    shiny::observeEvent(input$mark_btn, {
      ed <- plotly::event_data("plotly_selected", source = "subset_plot")
      if (is.null(ed) || input$group_name == "") {
        shiny::showNotification("   Please select cells and enter a name first!", type = "warning")
        return()
      }
      save_history()
      keys <- as.character(ed$key)
      
      new_label <- input$group_name
      if (input$parent_pop != "All_Cells") {
        clean_name <- gsub("^_|_$", "", input$group_name)
        new_label <- paste0(input$parent_pop, "_", clean_name)
      }
      
      values$data[keys, "Gate_Status"] <- new_label
      
      if (!(new_label %in% values$tree$Node)) {
        values$tree <- rbind(values$tree, data.frame(Node = new_label, Parent = input$parent_pop, stringsAsFactors = FALSE))
      }
      
      shiny::showNotification(paste0("  Tree node created: ", new_label), type = "message")
    })
    
    shiny::observeEvent(input$export_btn, {
      gated_cells <- values$data[values$data$Gate_Status != "Ungated", ]
      if (nrow(gated_cells) == 0) {
        shiny::showNotification("   No annotated subpopulations available for export!", type = "warning")
        return()
      }
      
      res_list <- list()
      unique_gates <- unique(gated_cells$Gate_Status)
      for (g in unique_gates) {
        sub_g <- gated_cells[gated_cells$Gate_Status == g, numeric_cols, drop = FALSE]
        if (nrow(sub_g) > 0) {
          q_low  <- apply(sub_g, 2, quantile, probs = 0.01, na.rm = TRUE)
          q_high <- apply(sub_g, 2, quantile, probs = 0.99, na.rm = TRUE)
          tmp_df <- data.frame(
            Tree_Node     = g,
            Channel       = numeric_cols,
            Min_Threshold = round(q_low,  3),
            Max_Threshold = round(q_high, 3),
            stringsAsFactors = FALSE
          )
          res_list[[g]] <- tmp_df
        }
      }
      tpl_df <- do.call(rbind, res_list)
      write.csv(tpl_df, "Gating_Thresholds_Template.csv", row.names = FALSE)
      shiny::showNotification("  Threshold template exported to the working directory! (using the 1%/99% quantiles)", type = "message")
    })
    
    shiny::observeEvent(input$done, {
      final_res <- values$data
      final_gate_status <- as.character(final_res$Gate_Status)
      
      message("\n=======================================================")
      message(">>>   Manual gating completed!")
      
      tree_out <- NULL
      if (nrow(values$tree) > 0) {
        tree_out <- values$tree
        
        #   Change: when exporting, also calculate cumulative cell counts including all child nodes
        tree_out$Cell_Count <- sapply(tree_out$Node, function(x) {
          sum(final_gate_status == x | startsWith(final_gate_status, paste0(x, "_")))
        })
        
        write.csv(tree_out, "Gating_Tree_Structure.csv", row.names = FALSE)
        message(">>>   [Hierarchy topology] Gating tree structure exported to: Gating_Tree_Structure.csv")
        
        message("\n  Cell Gating Hierarchy Tree (Gating Tree Topology):")
        message("All_Cells (", nrow(final_res), ")")
        
        print_tree <- function(parent, prefix = "") {
          children <- tree_out[tree_out$Parent == parent, ]
          if (nrow(children) > 0) {
            for (k in 1:nrow(children)) {
              is_last  <- (k == nrow(children))
              connector <- ifelse(is_last, "    ", "    ")
              message(prefix, connector, children$Node[k], " (", children$Cell_Count[k], ")")
              new_prefix <- paste0(prefix, ifelse(is_last, "    ", "    "))
              print_tree(children$Node[k], new_prefix)
            }
          }
        }
        print_tree("All_Cells")
        
      } else {
        message(">>>    No structured subpopulation nodes were generated.")
      }
      
      final_res <- Assign_Terminal_Gate_CellType(
        final_res,
        tree_df = tree_out,
        gate_col = "Gate_Status",
        output_col = "CellType",
        unknown_label = "Unknown"
      )
      final_res$Gate_Status <- NULL
      
      message("\n      Unknown/non-terminal cells: ", sum(final_res$CellType == "Unknown"), " cells")
      message("=======================================================\n")
      
      shiny::stopApp(final_res)
    })
  }
  
  shiny::runGadget(ui, server, viewer = shiny::dialogViewer("Hierarchy Gating Station", width = 1450, height = 920))
}


#' Draw and export the visualized gating hierarchy tree (PDF)
Plot_Gating_Tree <- function(csv_path = "Gating_Tree_Structure.csv", save_path = "00_Gating_Tree_Plot.pdf") {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("  Missing the igraph package")
  
  if (!file.exists(csv_path)) {
    message("   Gating topology file not found; this gating session may not have generated child levels.")
    return(NULL)
  }
  
  tree_df <- read.csv(csv_path)
  if (nrow(tree_df) == 0) return(NULL)
  
  edges <- tree_df[, c("Parent", "Node")]
  g <- igraph::graph_from_data_frame(edges, directed = TRUE)
  
  node_names <- igraph::V(g)$name
  node_labels <- sapply(node_names, function(x) {
    if (x %in% tree_df$Node) {
      count <- tree_df$Cell_Count[tree_df$Node == x]
      return(paste0(x, "\n(", count, ")"))
    } else {
      return(paste0(x, "\n(Root)"))
    }
  })
  
  pdf(save_path, width = 8, height = 6)
  par(mar = c(1, 1, 3, 1))
  
  plot(g,
       layout            = igraph::layout_as_tree(g, flip.y = FALSE),
       vertex.label      = node_labels,
       vertex.color      = "#E8F8F5",
       vertex.frame.color = "#1ABC9C",
       vertex.shape      = "rectangle",
       vertex.size       = 45,
       vertex.size2      = 25,
       vertex.label.color = "#2C3E50",
       vertex.label.cex  = 0.8,
       vertex.label.font = 2,
       edge.arrow.size   = 0.4,
       edge.color        = "#95A5A6",
       edge.width        = 2,
       main              = "Cell Gating Hierarchy Tree"
  )
  
  dev.off()
  message(">>>   Visualized gating hierarchy tree saved to: ", save_path)
}

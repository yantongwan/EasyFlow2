# ==============================================================================
# Module 4: AI Dual-Track Annotation Center - INTEGRATED VERSION
# Filename: R/04_Annotation.R
#
# Features:
# - Negative marker support in custom_rules
# - Extended state_name logic for multiple cell types
# - Fixed normalization (trusts LLM output, no forced overrides)
# - Duplicate word removal in fine_name
# ==============================================================================

Call_OpenRouter_Chat <- function(messages, model, api_key) {
  if (!nzchar(api_key)) stop("  OPENROUTER_API_KEY is not set")
  req <- httr2::request("https://openrouter.ai/api/v1/chat/completions") |>
    httr2::req_headers(
      Authorization = paste("Bearer", api_key),
      `Content-Type` = "application/json"
    ) |>
    httr2::req_body_json(
      list(model = model, messages = messages, temperature = 0.1),
      auto_unbox = TRUE
    ) |>
    httr2::req_error(is_error = function(resp) FALSE)

  resp <- httr2::req_perform(req)
  txt  <- httr2::resp_body_string(resp)
  obj  <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE), error = function(e) NULL)

  if (is.null(obj) || is.null(obj$choices)) {
    warning("   LLM API request failed: ", txt)
    return(NULL)
  }
  obj$choices[[1]]$message$content
}

.extract_json_array <- function(text) {
  text <- trimws(gsub("^```json\\n?|^```\\n?|\\n?```$", "", text))

  start <- regexpr("\\[", text)
  if (start == -1) return(NULL)

  depth <- 0L
  chars <- strsplit(text, "")[[1]]
  end   <- NA_integer_

  for (i in seq(start, length(chars))) {
    if (chars[i] == "[") depth <- depth + 1L
    if (chars[i] == "]") {
      depth <- depth - 1L
      if (depth == 0L) {
        end <- i
        break
      }
    }
  }

  if (is.na(end)) return(NULL)
  substr(text, start, end)
}

`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !all(is.na(a)) && any(nzchar(trimws(as.character(a))))) a else b
}

.clean_text_scalar <- function(x, default = "") {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(default)
  x <- as.character(x)[1]
  x <- trimws(gsub("\\s+", " ", x))
  if (!nzchar(x)) default else x
}

.to_title_case_loose <- function(x) {
  x <- .clean_text_scalar(x, "")
  if (!nzchar(x)) return("")
  parts <- strsplit(tolower(x), " ")[[1]]
  parts <- vapply(parts, function(p) {
    if (!nzchar(p)) return("")
    paste0(toupper(substr(p, 1, 1)), substr(p, 2, nchar(p)))
  }, character(1))
  gsub("\\s+", " ", trimws(paste(parts, collapse = " ")))
}

.contains_any <- function(text, patterns) {
  text <- tolower(.clean_text_scalar(text, ""))
  if (!nzchar(text)) return(FALSE)
  any(vapply(patterns, function(p) grepl(p, text, perl = TRUE), logical(1)))
}

Build_Cluster_Summary_Values <- function(df, marker_cols, group_col = "Cluster") {
  global_mean <- vapply(marker_cols, function(m) mean(df[[m]], na.rm = TRUE), numeric(1))
  global_sd   <- vapply(marker_cols, function(m) sd(df[[m]], na.rm = TRUE), numeric(1))
  global_q75  <- vapply(marker_cols, function(m) quantile(df[[m]], 0.75, na.rm = TRUE), numeric(1))

  clusters    <- sort(unique(df[[group_col]]))
  out         <- vector("list", length(clusters))
  names(out)  <- as.character(clusters)

  for (cl in clusters) {
    sub <- df[df[[group_col]] == cl, , drop = FALSE]
    cl_mean <- vapply(marker_cols, function(m) mean(sub[[m]], na.rm = TRUE), numeric(1))
    marker_stats <- list()
    for (m in marker_cols) {
      z_score <- if (global_sd[m] > 0) {
        (cl_mean[m] - global_mean[m]) / global_sd[m]
      } else {
        0
      }

      fold_change <- (cl_mean[m] + 0.01) / (global_mean[m] + 0.01)
      pct_positive <- mean(sub[[m]] > global_q75[m], na.rm = TRUE) * 100

      expr_level <- if (z_score > 2) {
        "very_high"
      } else if (z_score > 1) {
        "high"
      } else if (z_score > -1) {
        "medium"
      } else {
        "low"
      }

      marker_stats[[m]] <- list(
        cluster_avg = round(cl_mean[m], 3),
        global_avg  = round(global_mean[m], 3),
        z_score = round(z_score, 2),
        fold_change = round(fold_change, 2),
        pct_positive = round(pct_positive, 1),
        expr_level = expr_level
      )
    }
    out[[as.character(cl)]] <- list(
      cluster = as.character(cl),
      cell_count = nrow(sub),
      expression_profile = marker_stats
    )
  }
  out
}

# ==============================================================================
#   NEW: Enhanced rule matching with NEGATIVE markers
# ==============================================================================

Match_Rules_With_Negatives <- function(cluster_data, rules) {
  # cluster_data: expression_profile from Build_Cluster_Summary_Values
  # rules: list with positive and negative markers

  for (rule_name in names(rules)) {
    rule <- rules[[rule_name]]

    # Handle old format (just a vector of markers)
    if (is.character(rule)) {
      positive_markers <- rule
      negative_markers <- character(0)
    } else {
      positive_markers <- rule$positive %||% character(0)
      negative_markers <- rule$negative %||% character(0)
    }

    # Check positive markers
    positive_match <- if (length(positive_markers) > 0) {
      all(vapply(positive_markers, function(m) {
        if (!is.null(cluster_data[[m]])) {
          cluster_data[[m]]$cluster_avg > cluster_data[[m]]$global_avg
        } else {
          FALSE
        }
      }, logical(1)))
    } else {
      TRUE
    }

    # Check negative markers (should be LOW)
    negative_match <- if (length(negative_markers) > 0) {
      all(vapply(negative_markers, function(m) {
        if (!is.null(cluster_data[[m]])) {
          cluster_data[[m]]$cluster_avg <= cluster_data[[m]]$global_avg ||
            cluster_data[[m]]$z_score < 0
        } else {
          TRUE
        }
      }, logical(1)))
    } else {
      TRUE
    }

    if (positive_match && negative_match) {
      return(rule_name)
    }
  }

  return(NULL)
}

# ==============================================================================
#   IMPROVED: Normalize_Annotation_Label
# - Trusts LLM's major_lineage output (no forced overrides)
# - Extended state_name logic for multiple cell types
# - Removes duplicate words in fine_name
# ==============================================================================

Normalize_Annotation_Label <- function(
    major_lineage,
    minor_lineage = "",
    state_name = "",
    biomarker = "",
    reasoning = "",
    base_population_context = NULL,
    allow_functional_prefix = TRUE,
    ontology = NULL,
    unknown_threshold = NULL,
    confidence = NA_real_
) {
  major_lineage <- .clean_text_scalar(major_lineage, "Unknown")
  minor_lineage <- .clean_text_scalar(minor_lineage, "")
  state_name    <- .clean_text_scalar(state_name, "")
  biomarker     <- .clean_text_scalar(biomarker, "")
  reasoning     <- .clean_text_scalar(reasoning, "")
  base_population_context <- .clean_text_scalar(base_population_context, "")

  raw_concat <- paste(major_lineage, minor_lineage, state_name, biomarker, reasoning, base_population_context)
  raw_low    <- tolower(raw_concat)
  major_low  <- tolower(major_lineage)
  minor_low  <- tolower(minor_lineage)
  state_low  <- tolower(state_name)

  # Low confidence handling
  conf_val <- suppressWarnings(as.numeric(confidence))
  if (!is.na(conf_val) && !is.null(unknown_threshold) && conf_val < unknown_threshold) {
    if (!.contains_any(raw_low, c("foxp3", "cd3", "cd4", "cd8", "nk", "cd68", "f4_80", "ly6c", "cd11b", "cd11c", "mrc1", "cd206", "ifng", "tnfa"))) {
      major_lineage <- "Unknown"
      minor_lineage <- ""
      state_name    <- ""
    }
  }

  # Basic cleaning
  if (major_low %in% c("na", "null", "", "unknown", "unassigned")) major_lineage <- "Unknown"
  if (minor_low %in% c("na", "null", "unknown", "unassigned")) minor_lineage <- ""
  if (state_low %in% c("na", "null", "unknown", "unassigned")) state_name <- ""

  #   REMOVED: Auto-standardization logic that forced major_lineage changes
  # Now we trust LLM's major_lineage output

  major_low <- tolower(major_lineage)

  # Helper function to extract z-score
  extract_z_score <- function(marker_name, reasoning_text) {
    pattern <- paste0(marker_name, ".*?\\(.*?z[_=]?\\s*[:-]?\\s*(-?\\d+\\.?\\d*)")
    if (grepl(pattern, reasoning_text, ignore.case = TRUE, perl = TRUE)) {
      z_str <- sub(paste0(".*", pattern, ".*"), "\\1", reasoning_text, ignore.case = TRUE, perl = TRUE)
      z_val <- suppressWarnings(as.numeric(z_str))
      if (!is.na(z_val)) return(z_val)
    }
    return(NA)
  }

  # ===========================================================================
  #   CELL-TYPE SPECIFIC LOGIC (only processes minor_lineage and state_name)
  # ===========================================================================

  state_tokens <- character(0)

  # ---------------------------------------------------------------------------
  # CD4 T CELLS
  # ---------------------------------------------------------------------------
  if (major_low == "cd4 t cell") {
    # Subtype determination
    foxp3_z <- extract_z_score("foxp3", reasoning)
    ifng_z <- extract_z_score("ifng", reasoning)
    il17a_z <- extract_z_score("il17a", reasoning)
    il4_z <- extract_z_score("il4", reasoning)

    if (!is.na(il17a_z) && il17a_z > 1.5) {
      minor_lineage <- "Th17"
    } else if (!is.na(ifng_z) && ifng_z > 1.5) {
      minor_lineage <- "Th1"
    } else if (!is.na(il4_z) && il4_z > 1.5) {
      minor_lineage <- "Th2"
    } else if (!is.na(foxp3_z) && foxp3_z > 1.5) {
      minor_lineage <- "Treg"
    } else if (tolower(minor_lineage) == "conventional cd4 t cell") {
      minor_lineage <- ""
    }

    # Functional markers
    if (isTRUE(allow_functional_prefix)) {
      if (!is.na(ifng_z) && ifng_z > 1.0) state_tokens <- c(state_tokens, "IFNg+")
      if (!is.na(il17a_z) && il17a_z > 1.0) state_tokens <- c(state_tokens, "IL17+")
      if (!is.na(il4_z) && il4_z > 1.0) state_tokens <- c(state_tokens, "IL4+")

      tnfa_z <- extract_z_score("tnfa", reasoning)
      il2_z <- extract_z_score("il2", reasoning)
      if (!is.na(tnfa_z) && tnfa_z > 1.0) state_tokens <- c(state_tokens, "TNFa+")
      if (!is.na(il2_z) && il2_z > 1.0) state_tokens <- c(state_tokens, "IL2+")
    }

    # Differentiation states
    if (.contains_any(raw_low, c("effector memory", "tem", "cd44.*high.*cd62l.*low"))) {
      state_tokens <- c(state_tokens, "Effector Memory")
    } else if (.contains_any(raw_low, c("central memory", "tcm", "cd44.*high.*cd62l.*high"))) {
      state_tokens <- c(state_tokens, "Central Memory")
    } else if (.contains_any(raw_low, c("memory", "cd44.*high")) && !.contains_any(raw_low, c("naive"))) {
      state_tokens <- c(state_tokens, "Memory")
    }

    if (.contains_any(raw_low, c("naive", "cd62l.*high", "cd44.*low")) && !.contains_any(raw_low, c("memory"))) {
      state_tokens <- c(state_tokens, "Naive")
    }

    # Activation/exhaustion
    if (.contains_any(raw_low, c("activat", "cd69", "cd25.*high"))) {
      state_tokens <- c(state_tokens, "Activated")
    }
    if (.contains_any(raw_low, c("exhaust", "pd1", "pd-1", "lag3", "tim3", "tigit"))) {
      state_tokens <- c(state_tokens, "Exhausted")
    }

    # Remove Naive from Treg
    if (identical(minor_lineage, "Treg")) {
      state_tokens <- state_tokens[!state_tokens %in% c("Naive")]
    }
  }

  # ---------------------------------------------------------------------------
  # CD8 T CELLS
  # ---------------------------------------------------------------------------
  else if (major_low == "cd8 t cell") {
    # Subtype determination
    if (.contains_any(raw_low, c("exhaust", "pd1", "lag3", "tim3", "tigit"))) {
      if (!nzchar(minor_lineage)) minor_lineage <- "Exhausted CD8 T cell"
    } else if (.contains_any(raw_low, c("gzmb", "granzyme.*b", "prf1", "perforin", "cytotoxic"))) {
      if (!nzchar(minor_lineage)) minor_lineage <- "Cytotoxic CD8 T cell"
    } else if (tolower(minor_lineage) == "conventional cd8 t cell") {
      minor_lineage <- ""
    }

    # Functional markers
    if (isTRUE(allow_functional_prefix)) {
      gzmb_z <- extract_z_score("gzmb|granzyme.?b", reasoning)
      prf1_z <- extract_z_score("prf1|perforin", reasoning)
      ifng_z <- extract_z_score("ifng", reasoning)
      tnfa_z <- extract_z_score("tnfa", reasoning)

      if (!is.na(gzmb_z) && gzmb_z > 1.0) state_tokens <- c(state_tokens, "GzmB+")
      if (!is.na(prf1_z) && prf1_z > 1.0) state_tokens <- c(state_tokens, "Perforin+")
      if (!is.na(ifng_z) && ifng_z > 1.0) state_tokens <- c(state_tokens, "IFNg+")
      if (!is.na(tnfa_z) && tnfa_z > 1.0) state_tokens <- c(state_tokens, "TNFa+")
    }

    # Differentiation states
    if (.contains_any(raw_low, c("effector memory", "tem"))) {
      state_tokens <- c(state_tokens, "Effector Memory")
    } else if (.contains_any(raw_low, c("central memory", "tcm"))) {
      state_tokens <- c(state_tokens, "Central Memory")
    } else if (.contains_any(raw_low, c("memory", "cd44"))) {
      state_tokens <- c(state_tokens, "Memory")
    }

    if (.contains_any(raw_low, c("naive", "cd62l"))) {
      state_tokens <- c(state_tokens, "Naive")
    }

    if (.contains_any(raw_low, c("activat", "cd69"))) {
      state_tokens <- c(state_tokens, "Activated")
    }
    if (.contains_any(raw_low, c("exhaust", "pd1", "lag3", "tim3"))) {
      state_tokens <- c(state_tokens, "Exhausted")
    }
  }

  # ---------------------------------------------------------------------------
  # MACROPHAGES
  # ---------------------------------------------------------------------------
  else if (major_low == "macrophage") {
    # M1/M2 polarization (only if minor_lineage is empty)
    if (!nzchar(minor_lineage)) {
      if (.contains_any(raw_low, c("m2", "cd206", "mrc1", "arg1", "il4", "chil3", "ym1"))) {
        minor_lineage <- "M2 macrophage"
      } else if (.contains_any(raw_low, c("m1", "il1b", "tnfa", "ifng", "nos2", "inos", "cd86"))) {
        minor_lineage <- "M1 macrophage"
      }
    }

    # Functional markers
    if (isTRUE(allow_functional_prefix)) {
      inos_z <- extract_z_score("inos|nos2", reasoning)
      arg1_z <- extract_z_score("arg1", reasoning)
      il12_z <- extract_z_score("il12", reasoning)
      il10_z <- extract_z_score("il10", reasoning)

      if (!is.na(inos_z) && inos_z > 1.0) state_tokens <- c(state_tokens, "iNOS+")
      if (!is.na(arg1_z) && arg1_z > 1.0) state_tokens <- c(state_tokens, "Arg1+")
      if (!is.na(il12_z) && il12_z > 1.0) state_tokens <- c(state_tokens, "IL12+")
      if (!is.na(il10_z) && il10_z > 1.0) state_tokens <- c(state_tokens, "IL10+")
    }

    # States
    if (.contains_any(raw_low, c("activat", "cd80", "cd86"))) {
      state_tokens <- c(state_tokens, "Activated")
    }
    if (.contains_any(raw_low, c("inflammat"))) {
      state_tokens <- c(state_tokens, "Inflammatory")
    }
  }

  # ---------------------------------------------------------------------------
  # MONOCYTES
  # ---------------------------------------------------------------------------
  else if (major_low == "monocyte") {
    # Classical vs non-classical (only if minor_lineage is empty)
    if (!nzchar(minor_lineage)) {
      if (.contains_any(raw_low, c("ly6c.*high", "ly6chi", "classical"))) {
        minor_lineage <- "Classical monocyte"
      } else if (.contains_any(raw_low, c("ly6c.*low", "ly6clo", "non.?classical"))) {
        minor_lineage <- "Non-classical monocyte"
      }
    }

    # States
    if (.contains_any(raw_low, c("inflammat"))) {
      state_tokens <- c(state_tokens, "Inflammatory")
    }
    if (.contains_any(raw_low, c("activat"))) {
      state_tokens <- c(state_tokens, "Activated")
    }
  }

  # ---------------------------------------------------------------------------
  # NEUTROPHILS
  # ---------------------------------------------------------------------------
  else if (major_low == "neutrophil") {
    if (.contains_any(raw_low, c("activat", "cd11b.*high"))) {
      state_tokens <- c(state_tokens, "Activated")
    }
    if (.contains_any(raw_low, c("inflammat"))) {
      state_tokens <- c(state_tokens, "Inflammatory")
    }
  }

  # ---------------------------------------------------------------------------
  # EOSINOPHILS
  # ---------------------------------------------------------------------------
  else if (major_low == "eosinophil") {
    if (.contains_any(raw_low, c("activat"))) {
      state_tokens <- c(state_tokens, "Activated")
    }
    if (.contains_any(raw_low, c("inflammat"))) {
      state_tokens <- c(state_tokens, "Inflammatory")
    }
  }

  # ---------------------------------------------------------------------------
  # DENDRITIC CELLS
  # ---------------------------------------------------------------------------
  else if (major_low == "dendritic cell") {
    # DC subtypes (only if minor_lineage is empty)
    if (!nzchar(minor_lineage)) {
      if (.contains_any(raw_low, c("cdc1", "cd8a", "xcr1"))) {
        minor_lineage <- "cDC1"
      } else if (.contains_any(raw_low, c("cdc2", "cd11b.*high"))) {
        minor_lineage <- "cDC2"
      } else if (.contains_any(raw_low, c("pdc", "plasmacytoid"))) {
        minor_lineage <- "pDC"
      }
    }

    # States
    if (.contains_any(raw_low, c("mature", "cd80.*high", "cd86.*high"))) {
      state_tokens <- c(state_tokens, "Mature")
    } else if (.contains_any(raw_low, c("immature"))) {
      state_tokens <- c(state_tokens, "Immature")
    }
    if (.contains_any(raw_low, c("activat"))) {
      state_tokens <- c(state_tokens, "Activated")
    }
  }

  # ---------------------------------------------------------------------------
  # B CELLS
  # ---------------------------------------------------------------------------
  else if (major_low == "b cell") {
    # Subtypes (only if minor_lineage is empty)
    if (!nzchar(minor_lineage)) {
      if (.contains_any(raw_low, c("plasma", "cd138"))) {
        minor_lineage <- "Plasma cell"
      } else if (.contains_any(raw_low, c("germinal center", "gc b"))) {
        minor_lineage <- "Germinal center B cell"
      }
    }

    # States
    if (.contains_any(raw_low, c("memory", "cd27"))) {
      state_tokens <- c(state_tokens, "Memory")
    }
    if (.contains_any(raw_low, c("naive"))) {
      state_tokens <- c(state_tokens, "Naive")
    }
    if (.contains_any(raw_low, c("activat"))) {
      state_tokens <- c(state_tokens, "Activated")
    }
  }

  # ---------------------------------------------------------------------------
  # NK CELLS
  # ---------------------------------------------------------------------------
  else if (major_low == "nk cell") {
    if (.contains_any(raw_low, c("activat"))) {
      state_tokens <- c(state_tokens, "Activated")
    }
    if (.contains_any(raw_low, c("cytotoxic", "gzmb", "perforin"))) {
      state_tokens <- c(state_tokens, "Cytotoxic")
    }
  }

  # ---------------------------------------------------------------------------
  # UNIVERSAL: Ki67+ (proliferation)
  # ---------------------------------------------------------------------------
  if (isTRUE(allow_functional_prefix)) {
    ki67_z <- extract_z_score("ki67", reasoning)
    if (!is.na(ki67_z) && ki67_z > 1.0) {
      state_tokens <- c(state_tokens, "Ki67+")
    }
  }

  # ===========================================================================
  # Finalize state_name
  # ===========================================================================
  state_tokens <- unique(state_tokens[nzchar(state_tokens)])

  # Prevent contradictory states
  if ("Naive" %in% state_tokens && any(c("Memory", "Effector Memory", "Central Memory") %in% state_tokens)) {
    state_tokens <- state_tokens[state_tokens != "Naive"]
  }
  if ("Effector Memory" %in% state_tokens || "Central Memory" %in% state_tokens) {
    state_tokens <- state_tokens[state_tokens != "Memory"]
  }
  if ("Naive" %in% state_tokens && "Activated" %in% state_tokens) {
    state_tokens <- state_tokens[state_tokens != "Activated"]
  }

  # Preferred order
  preferred_state_order <- c(
    "Ki67+",
    "IFNg+", "TNFa+", "IL2+", "IL4+", "IL17+", "GzmB+", "Perforin+",
    "iNOS+", "Arg1+", "IL12+", "IL10+",
    "Activated", "Inflammatory", "Mature", "Immature",
    "Exhausted", "Cytotoxic",
    "Effector Memory", "Central Memory", "Memory", "Naive"
  )

  state_tokens <- preferred_state_order[preferred_state_order %in% state_tokens]

  # Construct final names
  coarse_name <- major_lineage
  base_fine   <- if (nzchar(minor_lineage)) minor_lineage else major_lineage
  state_name  <- paste(state_tokens, collapse = " ")
  state_name  <- trimws(gsub("\\s+", " ", state_name))

  fine_name <- if (nzchar(state_name)) paste(state_name, base_fine) else base_fine
  fine_name <- trimws(gsub("\\s+", " ", fine_name))

  #   Remove duplicate words in fine_name
  fine_name_words <- strsplit(fine_name, " ")[[1]]
  fine_name_lower <- tolower(fine_name_words)
  if (any(duplicated(fine_name_lower))) {
    # Keep first occurrence of each word (case-insensitive)
    keep_idx <- !duplicated(fine_name_lower)
    fine_name <- paste(fine_name_words[keep_idx], collapse = " ")
  }

  # Ontology validation
  if (!is.null(ontology) && length(ontology) > 0 && !(coarse_name %in% c(ontology, "Unknown"))) {
    coarse_name <- "Unknown"
    if (!nzchar(minor_lineage)) minor_lineage <- ""
    if (!nzchar(state_name)) state_name <- ""
    fine_name <- if (nzchar(state_name)) paste(state_name, coarse_name) else coarse_name
  }

  list(
    major_lineage = coarse_name,
    minor_lineage = minor_lineage,
    state_name    = state_name,
    coarse_name   = coarse_name,
    fine_name     = fine_name
  )
}

Run_Agentic_LLM_Loop <- function(
    cluster_payload,
    model,
    api_key,
    ontology = NULL,
    base_population_context = NULL,
    max_clusters_per_batch = 20,
    allow_functional_prefix = TRUE,
    unknown_threshold = 35
) {
  onto_prompt <- ""
  if (!is.null(ontology) && length(ontology) > 0) {
    onto_prompt <- paste(
      "CRITICAL ONTOLOGY RULE:",
      "You MUST select the major_lineage STRICTLY from the following controlled vocabulary list:",
      paste0("['", paste(ontology, collapse = "', '"), "']"),
      "minor_lineage can be more refined, but it must remain biologically compatible with the selected major_lineage.",
      "If completely unidentifiable, output major_lineage = 'Unknown'.",
      sep = "\n"
    )
  }

  context_prompt <- ""
  if (!is.null(base_population_context) && nzchar(base_population_context)) {
    context_prompt <- paste(
      "=========================================================================",
      paste0("ABSOLUTE GROUND TRUTH (BASE CONTEXT): The user has already sorted/pre-gated this entire dataset as [ ", base_population_context, " ]."),
      "Therefore, you MUST assume EVERY single cluster in this dataset inherently belongs to this gated biological universe.",
      "Do not reject a lineage merely because some markers are below global average.",
      "Use absolute cluster_avg, marker combinations, and biological compatibility.",
      "=========================================================================",
      sep = "\n"
    )
  }

  payload_list <- tryCatch(
    jsonlite::fromJSON(cluster_payload, simplifyVector = FALSE),
    error = function(e) {
      warning("   cluster_payload failed to parse; sending as one whole batch.")
      NULL
    }
  )

  if (!is.null(payload_list) && length(payload_list) > max_clusters_per_batch) {
    n_batches <- ceiling(length(payload_list) / max_clusters_per_batch)
    message(sprintf(
      "     Detected %d clusters; sending in %d batches (each batch <= %d)...",
      length(payload_list), n_batches, max_clusters_per_batch
    ))

    all_results <- list()
    batch_names <- names(payload_list)

    for (b in seq_len(n_batches)) {
      idx_start <- (b - 1) * max_clusters_per_batch + 1
      idx_end   <- min(b * max_clusters_per_batch, length(payload_list))
      batch_list <- payload_list[idx_start:idx_end]
      batch_payload <- jsonlite::toJSON(batch_list, auto_unbox = TRUE, pretty = TRUE)

      message(sprintf(
        "   [Batch %d/%d] processing cluster %s ~ %s...",
        b, n_batches, batch_names[idx_start], batch_names[idx_end]
      ))

      batch_result <- .run_single_llm_batch(
        cluster_payload         = batch_payload,
        model                   = model,
        api_key                 = api_key,
        onto_prompt             = onto_prompt,
        context_prompt          = context_prompt,
        base_population_context = base_population_context,
        allow_functional_prefix = allow_functional_prefix,
        unknown_threshold       = unknown_threshold,
        ontology                = ontology
      )

      if (!is.null(batch_result)) {
        all_results[[length(all_results) + 1]] <- batch_result
      }
    }

    if (length(all_results) == 0) return(NULL)
    return(do.call(rbind, all_results))
  }

  .run_single_llm_batch(
    cluster_payload         = cluster_payload,
    model                   = model,
    api_key                 = api_key,
    onto_prompt             = onto_prompt,
    context_prompt          = context_prompt,
    base_population_context = base_population_context,
    allow_functional_prefix = allow_functional_prefix,
    unknown_threshold       = unknown_threshold,
    ontology                = ontology
  )
}

.run_single_llm_batch <- function(
    cluster_payload,
    model,
    api_key,
    onto_prompt = "",
    context_prompt = "",
    base_population_context = NULL,
    allow_functional_prefix = TRUE,
    unknown_threshold = 35,
    ontology = NULL
) {
  base_context_text <- if (!is.null(base_population_context) && nzchar(base_population_context)) {
    paste0(
      "USER-PROVIDED GATING CONTEXT:\n",
      base_population_context, "\n",
      "You MUST treat this as a strong biological prior. ",
      "Do not assign cell types outside this gated universe unless the marker profile is overwhelmingly contradictory.\n"
    )
  } else {
    paste0(
      "USER-PROVIDED GATING CONTEXT:\n",
      "No additional gated context was provided. Use the marker profile itself.\n"
    )
  }

  message("   [Agent 1: Primary annotator] Generating hierarchical draft labels...")

  sys_actor <- paste(
    "You are an AI assistant specialized in high-dimensional flow cytometry annotation.",
    "",
    "  CRITICAL: Each marker now includes these metrics:",
    "- cluster_avg: mean expression in this cluster",
    "- global_avg: mean expression across all cells",
    "- z_score: standardized score (PRIMARY decision metric)",
    "- fold_change: cluster_avg / global_avg",
    "- pct_positive: % of cells above 75th percentile",
    "",
    "Your task is to annotate each cluster using THREE layers:",
    "1) major_lineage: broad cell class such as CD4 T cell, CD8 T cell, NK cell, Macrophage, Monocyte",
    "2) minor_lineage: biologically refined subtype within the major lineage, such as Treg, Th1-like CD4 T cell, M1 macrophage, M2 macrophage",
    "3) state_name: functional/state modifier such as IFNg+, TNFa+, Ki67+, Activated, Memory, Exhausted",
    "",
    "  DECISION RULES (MOST IMPORTANT):",
    "1. Use z_score as PRIMARY metric:",
    "   - z_score > 2.0 = VERY HIGH (strong positive marker)",
    "   - z_score > 1.5 = HIGH (positive marker)",
    "   - z_score > 1.0 = MODERATE (weak positive)",
    "   - z_score < 0 = NEGATIVE",
    "",
    "2. For PRIMARY lineage markers (CD3, CD4, CD8, CD19, F4/80, CD68, CD11b, etc.):",
    "   - REQUIRE z_score > 1.5",
    "   - If z_score < 1.5, DO NOT assign that lineage",
    "",
    "3. For SUBTYPE markers (Foxp3, Tbet, CD206, Arg1, etc.):",
    "   - REQUIRE z_score > 1.0",
    "",
    "4.   CRITICAL - EXCLUSION/LOW-EXPRESSION markers (equally important as high markers):",
    "   Many cell types are defined by BOTH high AND low markers:",
    "   ",
    "   CD4 T cell subtypes:",
    "   - Naive CD4: CD62L z > 1.0, CD44 z < 0, CD25 z < 0",
    "   - Memory CD4: CD44 z > 1.0, CD62L z < 0 (for effector memory)",
    "   - Treg (canonical): Foxp3 z > 1.5, IFNg z < 0.5, IL17a z < 0.5",
    "   - Th1: IFNg z > 1.5, IL17a z < 0.5, Foxp3 z < 1.5",
    "   - Th17: IL17a z > 1.5, IFNg z < 0.5, Foxp3 z < 1.5",
    "   ",
    "   Macrophage subtypes:",
    "   - M1: iNOS z > 1.0, IL12 z > 1.0, CD206 z < 0, Arg1 z < 0",
    "   - M2: CD206 z > 1.0, Arg1 z > 1.0, iNOS z < 0, IL12 z < 0",
    "   ",
    "   Lineage exclusion:",
    "   - CD4 T cell: CD8 z < 0",
    "   - CD8 T cell: CD4 z < 0",
    "   - T cell: CD19 z < 0",
    "   - B cell: CD3 z < 0",
    "   - Macrophage: CD3 z < 0",
    "   ",
    "      If exclusion markers are NOT low, reduce confidence by 30-50 points or reconsider lineage",
    "",
    "5. Check EXCLUSION markers (legacy rule, see rule 4 for details):",
    "   - CD4 T cells MUST have CD8 z_score < 0",
    "   - CD8 T cells MUST have CD4 z_score < 0",
    "   - T cells MUST have CD19 z_score < 0",
    "   - B cells MUST have CD3 z_score < 0",
    "   - Macrophages MUST have CD3 z_score < 0",
    "",
    "6. If NO marker has z_score > 1.5:",
    "   - Output major_lineage = 'Unknown'",
    "   - Set confidence < 40",
    "",
    "  EXAMPLES (note both HIGH and LOW markers):",
    "",
    "Example 1 - Clear Treg (note low IFNg/IL17a):",
    "CD3: {z_score: 3.2, expr_level: 'very_high'}",
    "CD4: {z_score: 2.8, expr_level: 'very_high'}",
    "Foxp3: {z_score: 2.5, expr_level: 'very_high'}",
    "CD25: {z_score: 1.8, expr_level: 'high'}",
    "CD8: {z_score: -0.5, expr_level: 'low'}",
    "IFNg: {z_score: -0.3, expr_level: 'low'}",
    "IL17a: {z_score: -0.4, expr_level: 'low'}",
    "  major_lineage='CD4 T cell', minor_lineage='Treg', state_name='', confidence=95",
    "Reasoning: Strong CD3+CD4+Foxp3+ signature with CD25. CD8/IFNg/IL17a all low confirms canonical Treg",
    "",
    "Example 2 - M2 Macrophage (note low iNOS):",
    "F4_80: {z_score: 3.5, expr_level: 'very_high'}",
    "CD68: {z_score: 2.2, expr_level: 'very_high'}",
    "CD206: {z_score: 2.8, expr_level: 'very_high'}",
    "Arg1: {z_score: 1.5, expr_level: 'high'}",
    "iNOS: {z_score: -0.5, expr_level: 'low'}",
    "  major_lineage='Macrophage', minor_lineage='M2 macrophage', state_name='', confidence=90",
    "Reasoning: F4/80+CD68+ macrophage with high CD206 and Arg1, low iNOS indicates M2 polarization",
    "",
    "Example 3 - Activated Th1:",
    "CD3: {z_score: 3.0, expr_level: 'very_high'}",
    "CD4: {z_score: 2.5, expr_level: 'very_high'}",
    "Tbet: {z_score: 2.0, expr_level: 'very_high'}",
    "IFNg: {z_score: 2.3, expr_level: 'very_high'}",
    "Ki67: {z_score: 1.8, expr_level: 'high'}",
    "Foxp3: {z_score: -0.8, expr_level: 'low'}",
    "  major_lineage='CD4 T cell', minor_lineage='Th1', state_name='IFNg+ Ki67+', confidence=92",
    "Reasoning: CD4+ T cell with Tbet and IFNg expression indicates Th1, Ki67+ shows proliferation",
    "",
    "Example 4 - Th17 (NOT Treg despite being CD4+):",
    "CD3: {z_score: 2.8, expr_level: 'very_high'}",
    "CD4: {z_score: 2.2, expr_level: 'very_high'}",
    "IL17a: {z_score: 2.7, expr_level: 'very_high'}",
    "Foxp3: {z_score: 0.7, expr_level: 'medium'}",
    "IFNg: {z_score: 0.1, expr_level: 'medium'}",
    "  major_lineage='CD4 T cell', minor_lineage='Th17', state_name='IL17a+', confidence=90",
    "Reasoning: CD4+ T cell with very high IL17a. Foxp3 z_score < 1.5 so NOT Treg. This is Th17.",
    "",
    "Example 5 - CD4 T cell without clear subtype (LEAVE minor_lineage EMPTY):",
    "CD3: {z_score: 2.5, expr_level: 'very_high'}",
    "CD4: {z_score: 1.8, expr_level: 'high'}",
    "Foxp3: {z_score: 0.3, expr_level: 'medium'}",
    "IFNg: {z_score: -0.4, expr_level: 'low'}",
    "IL17a: {z_score: -0.3, expr_level: 'low'}",
    "  major_lineage='CD4 T cell', minor_lineage='', state_name='', confidence=70",
    "Reasoning: CD4+ T cell without strong subtype markers. Foxp3 z < 1.5 so NOT Treg. No Th1/Th17 signature. Leave minor_lineage empty.",
    "",
    "Example 6 - Naive CD4 T cell (note LOW CD44 and CD25):",
    "CD3: {z_score: 2.8, expr_level: 'very_high'}",
    "CD4: {z_score: 2.5, expr_level: 'very_high'}",
    "CD62L: {z_score: 2.2, expr_level: 'very_high'}",
    "CD44: {z_score: -1.5, expr_level: 'low'}",
    "CD25: {z_score: -0.8, expr_level: 'low'}",
    "Foxp3: {z_score: -0.5, expr_level: 'low'}",
    "IFNg: {z_score: -0.6, expr_level: 'low'}",
    "  major_lineage='CD4 T cell', minor_lineage='', state_name='Naive', confidence=88",
    "Reasoning: CD4+ T cell with high CD62L and LOW CD44/CD25, classic naive signature. Foxp3/IFNg low rules out Treg/Th1. Use 'Naive' in state_name.",
    "",
    "Example 7 - Ambiguous case (IMPORTANT):",
    "CD11c: {z_score: 1.2, expr_level: 'high'}",
    "CD11b: {z_score: 0.8, expr_level: 'medium'}",
    "F4_80: {z_score: 0.5, expr_level: 'medium'}",
    "CD68: {z_score: 0.3, expr_level: 'medium'}",
    "  major_lineage='Unknown', minor_lineage='', state_name='', confidence=30",
    "Reasoning: CD11c z_score < 1.5, insufficient for DC assignment. No clear lineage markers above threshold",
    "",
    base_context_text,
    context_prompt,
    onto_prompt,
    "",
    "STRICT ANNOTATION PRINCIPLES:",
    "1. major_lineage must stay broad and conservative.",
    "2. minor_lineage may be more specific, but it must remain biologically coherent and compatible with the major_lineage.",
    "3. state_name should capture functional/activity modifiers and should not redefine lineage.",
    "4. Do NOT use a single marker alone to define lineage whenever lineage requires combinational evidence.",
    "5. Use marker combinations, biological compatibility, and exclusion logic.",
    "6. If the cluster belongs to a gated parent population, stay within that biological universe unless strongly contradicted.",
    "7.    CRITICAL - Treg assignment rules:",
    "   - Treg REQUIRES Foxp3 z_score > 1.5 (not just > 0)",
    "   - If Foxp3 z_score < 1.5, DO NOT assign 'Treg' as minor_lineage",
    "   - If IFNg z_score > 1.5, assign 'Th1' or 'Th1-like', NOT Treg",
    "   - If IL17a z_score > 1.5, assign 'Th17', NOT Treg",
    "   - If TNFa z_score > 1.5 with low Foxp3, assign 'Th1-like', NOT Treg",
    "8.    CRITICAL - state_name rules:",
    "   - ONLY include markers in state_name if z_score > 1.0 OR explicitly marked as 'high'/'very_high'",
    "   - If IFNg z_score < 1.0, DO NOT add 'IFNg+' to state_name",
    "   - If TNFa z_score < 1.0, DO NOT add 'TNFa+' to state_name",
    "   - If IL17a z_score < 1.0, DO NOT add 'IL17+' to state_name",
    "   - If a marker is 'low' or 'negative', NEVER include it in state_name",
    "9.    CRITICAL - minor_lineage for CD4 T cells without clear subtype:",
    "   - If no subtype marker (Foxp3, IFNg, IL17a, etc.) has z_score > 1.5, LEAVE minor_lineage EMPTY",
    "   - DO NOT use 'Conventional CD4 T cell' or similar generic terms",
    "   - Empty minor_lineage is acceptable and preferred over vague labels",
    "10. Memory/naive/activated/exhausted usually belong to state_name unless they are essential to subtype identity.",
    "11. 'Unknown' should be used only when the profile is too weak, blank, or biologically irreconcilable.",
    "",
    "FUNCTIONAL PREFIX RULES:",
    if (isTRUE(allow_functional_prefix)) {
      paste(
        "If a functional marker is clearly elevated and biologically meaningful, you MAY represent it in state_name.",
        "Examples: IFNg+, TNFa+, Ki67+, IL4+, Activated, Memory.",
        "Do NOT include trivial markers that add no biological meaning."
      )
    } else {
      "Do NOT output cytokine/proliferation prefixes in state_name."
    },
    "",
    "OUTPUT FORMAT MUST be a simple JSON array containing objects strictly with these keys:",
    " - 'cluster'",
    " - 'major_lineage'",
    " - 'minor_lineage'",
    " - 'state_name'",
    " - 'coarse_name'",
    " - 'fine_name'",
    " - 'confidence'",
    " - 'Biomarker'",
    " - 'reasoning'",
    sep = "\n"
  )

  user_actor <- paste0(
    "Here is the cluster data:\n\n",
    cluster_payload,
    "\n\nPlease output the JSON array only."
  )

  draft_response <- Call_OpenRouter_Chat(
    messages = list(
      list(role = "system", content = sys_actor),
      list(role = "user", content = user_actor)
    ),
    model = model,
    api_key = api_key
  )
  if (is.null(draft_response)) return(NULL)

  message("   [Agent 2: Senior reviewer] Performing hierarchy correction and name normalization...")

  sys_critic <- paste(
    "You are a Senior Principal Investigator in Immunology reviewing flow cytometry annotations.",
    "",
    base_context_text,
    context_prompt,
    onto_prompt,
    "",
    "YOUR JOB:",
    "Review the junior draft and correct lineage errors, subtype errors, state errors, and overconfident hallucinations.",
    "",
    "  VALIDATION CHECKLIST (CHECK EACH ANNOTATION):",
    "",
    "1.   SINGLE-MARKER LINEAGE ERROR:",
    "   - Problem: Lineage assigned based on ONE marker when multiple are required",
    "   - Example: CD3 z_score=2.1 alone   'T cell' (WRONG if CD4/CD8 both low)",
    "   - Fix: Check ALL lineage markers. If only 1 is high, reduce confidence or mark 'Unknown'",
    "",
    "2.   EXCLUSION MARKER VIOLATION:",
    "   - Problem: Assigned lineage contradicts exclusion markers",
    "   - Example: 'CD4 T cell' but CD8 z_score=2.5 (WRONG - CD4/CD8 double positive rare)",
    "   - Fix: Check mutual exclusivity. CD4 high + CD8 high = likely artifact or DP thymocyte",
    "",
    "2b.   MISSING LOW-EXPRESSION REQUIREMENTS:",
    "   - Problem: Subtype assigned without checking that exclusion markers are LOW",
    "   - Examples:",
    "     * 'Treg' but IFNg z=1.8, IL17a z=1.5 (WRONG - should be Th1 or Th17, not Treg)",
    "     * 'M2 Macrophage' but iNOS z=2.0 (WRONG - M2 requires iNOS LOW)",
    "     * 'Naive CD4' but CD44 z=1.5 (WRONG - Naive requires CD44 LOW)",
    "   - Fix: Check that exclusion markers have z_score < 0.5 (or < 0 for strict exclusion)",
    "",
    "3.   FUNCTIONAL MARKER AS LINEAGE:",
    "   - Problem: State/function markers (TNFa, IFNg, IL4, IL17, Ki67, PD1, GranzymeB) used as lineage",
    "   - Example: 'IFNg+ cell' as major_lineage (WRONG)",
    "   - Fix: Move to state_name. Find true lineage from CD3/CD4/CD8/CD19/F4-80/CD11b/etc.",
    "",
    "4.   OVERCONFIDENCE WITH WEAK SIGNAL:",
    "   - Problem: High confidence (>80) but key markers have z_score < 1.5",
    "   - Example: confidence=90 but CD3 z_score=1.2, CD4 z_score=0.8 (WRONG)",
    "   - Fix: Reduce confidence to 40-60 or mark 'Unknown' if z_score < 1.0 for primary markers",
    "",
    "5.   HIERARCHY INCONSISTENCY:",
    "   - Problem: minor_lineage incompatible with major_lineage",
    "   - Example: major='B cell', minor='Th1' (WRONG)",
    "   - Fix: Ensure minor is a valid subtype of major",
    "",
    "6.   OVER-ASSIGNMENT OF SINGLE SUBTYPE:",
    "   - Problem: Too many clusters assigned to the SAME minor_lineage (e.g., all clusters   'Treg')",
    "   - Example: 10 clusters all labeled 'Treg' but only 2 have Foxp3 z_score > 1.5 (WRONG)",
    "   - Fix: Re-check each cluster's defining markers. Assign diverse subtypes based on actual marker patterns:",
    "     * Foxp3 z_score > 1.5   'Treg'",
    "     * IFNg z_score > 1.5 (Foxp3 low)   'Th1' or 'Th1-like'",
    "     * IL17a z_score > 1.5 (Foxp3 low)   'Th17'",
    "     * TNFa z_score > 1.5 (Foxp3 low, IFNg low)   'Activated CD4 T cell' or 'Th1-like'",
    "     * No strong subtype markers   'Conventional CD4 T cell' or leave minor_lineage empty",
    "",
    "  USE Z-SCORE FOR VALIDATION:",
    "- Primary lineage markers (CD3, CD4, CD8, CD19, F4/80, CD68, CD11b) MUST have z_score > 1.5",
    "- If z_score < 1.5 for claimed lineage marker   REDUCE confidence by 20-40 points",
    "- If z_score < 1.0 for claimed lineage marker   Consider 'Unknown' or alternative lineage",
    "",
    "  CRITICAL - Treg-specific validation:",
    "- Treg assignment REQUIRES Foxp3 z_score > 1.5",
    "- If Foxp3 z_score < 1.5, DO NOT assign 'Treg' - use 'Conventional CD4 T cell' or other appropriate subtype",
    "- If reasoning mentions 'Th1', 'Th17', or 'IFNg+/IL17a+' but minor_lineage='Treg', this is WRONG - correct it",
    "",
    "SUBTYPE NORMALIZATION RULES:",
    "1. Foxp3-high (z_score > 1.5) CD4 T-cell clusters should be normalized to minor_lineage = 'Treg' or 'Treg-like'.",
    "2. IFNg, TNFa, IL4, IL17A, Ki67, PD1, GranzymeB are usually state modifiers rather than lineage names.",
    "3. Memory, naive, activated, proliferating, exhausted should usually be encoded in state_name unless they are part of a well-established subtype name.",
    "4. coarse_name must equal major_lineage.",
    "5. fine_name should be constructed as state_name + minor_lineage when minor_lineage is available; otherwise state_name + major_lineage.",
    "6. Prefer standardized immunology labels over literal marker strings whenever possible.",
    "",
    "CONFIDENCE RULES:",
    paste0("If evidence is weak, reduce confidence. If confidence < ", unknown_threshold,
           ", you may output 'Unknown' only when no coherent lineage can be supported."),
    "",
    "Output the FINAL corrected JSON array containing objects strictly with keys:",
    " - 'cluster'",
    " - 'major_lineage'",
    " - 'minor_lineage'",
    " - 'state_name'",
    " - 'coarse_name'",
    " - 'fine_name'",
    " - 'confidence'",
    " - 'Biomarker'",
    " - 'reasoning'",
    sep = "\n"
  )

  user_critic <- paste0(
    "RAW CLUSTER DATA:\n", cluster_payload,
    "\n\nJUNIOR ASSISTANT'S DRAFT:\n", draft_response,
    "\n\nPlease review, correct biological errors, normalize the naming hierarchy, and output ONLY the FINAL JSON array."
  )

  final_response <- Call_OpenRouter_Chat(
    messages = list(
      list(role = "system", content = sys_critic),
      list(role = "user", content = user_critic)
    ),
    model = model,
    api_key = api_key
  )
  if (is.null(final_response)) return(NULL)

  json_str <- .extract_json_array(final_response)
  if (is.null(json_str)) {
    warning("   Unable to extract a valid JSON array from the LLM response.")
    return(NULL)
  }

  parsed_final <- tryCatch(jsonlite::fromJSON(json_str), error = function(e) NULL)
  if (is.null(parsed_final)) {
    warning("   JSON parsing failed.")
    return(NULL)
  }

  needed_cols <- c(
    "cluster", "major_lineage", "minor_lineage", "state_name",
    "coarse_name", "fine_name", "confidence", "Biomarker", "reasoning"
  )
  for (cc in needed_cols) {
    if (!cc %in% colnames(parsed_final)) parsed_final[[cc]] <- NA
  }

  if (!nrow(parsed_final)) return(NULL)

  normalized_rows <- lapply(seq_len(nrow(parsed_final)), function(k) {
    norm_res <- Normalize_Annotation_Label(
      major_lineage = parsed_final$major_lineage[k],
      minor_lineage = parsed_final$minor_lineage[k],
      state_name    = parsed_final$state_name[k],
      biomarker     = parsed_final$Biomarker[k],
      reasoning     = parsed_final$reasoning[k],
      base_population_context = base_population_context,
      allow_functional_prefix = allow_functional_prefix,
      ontology = ontology,
      unknown_threshold = unknown_threshold,
      confidence = parsed_final$confidence[k]
    )

    data.frame(
      cluster       = as.character(parsed_final$cluster[k]),
      major_lineage = norm_res$major_lineage,
      minor_lineage = norm_res$minor_lineage,
      state_name    = norm_res$state_name,
      coarse_name   = norm_res$coarse_name,
      fine_name     = norm_res$fine_name,
      confidence    = suppressWarnings(as.numeric(parsed_final$confidence[k])),
      Biomarker     = .clean_text_scalar(parsed_final$Biomarker[k], ""),
      reasoning     = .clean_text_scalar(parsed_final$reasoning[k], ""),
      stringsAsFactors = FALSE
    )
  })

  parsed_final <- do.call(rbind, normalized_rows)

  message("\n   [Review summary]:")
  for (k in seq_len(nrow(parsed_final))) {
    conf_val <- suppressWarnings(as.numeric(parsed_final$confidence[k]))
    final_nm <- parsed_final$fine_name[k]
    status_icon <- if (identical(parsed_final$major_lineage[k], "Unknown")) {
      "?"
    } else if (!is.na(conf_val) && conf_val < 60) {
      "~"
    } else {
      "+"
    }

    conf_str <- if (!is.na(conf_val)) paste0(" [Confidence: ", conf_val, "%]") else ""
    subtype_str <- if (nzchar(parsed_final$minor_lineage[k])) paste0(" | subtype=", parsed_final$minor_lineage[k]) else ""
    state_str   <- if (nzchar(parsed_final$state_name[k])) paste0(" | state=", parsed_final$state_name[k]) else ""

    message(sprintf(
      "     - [%s] Cluster %s -> coarse=%s | fine=%s%s%s%s: %s",
      status_icon,
      parsed_final$cluster[k],
      parsed_final$coarse_name[k],
      final_nm,
      conf_str,
      subtype_str,
      state_str,
      parsed_final$reasoning[k]
    ))
  }

  parsed_final
}

Annotate_Clusters_Dual_Track <- function(
    res_obj,
    rules = NULL,
    model = "openai/gpt-4o",
    use_llm_only = FALSE,
    ontology = NULL,
    max_clusters_per_batch = 20,
    base_population_context = NULL,
    allow_functional_prefix = TRUE,
    unknown_threshold = 35,
    annotation_granularity = c("fine", "coarse", "both")
) {
  annotation_granularity <- match.arg(annotation_granularity)
  df <- res_obj$Data

  ignore_cols <- c(
    "tSNE1", "tSNE2", "UMAP1", "UMAP2", "Cluster", "SampleID", "Group",
    "Time", "TIME", "CellType", "id", "ClusterLabel", "Pseudotime", "Gate_Status",
    "Cluster_Name", "Lineage_Name", "State_Name", "Major_Lineage", "Minor_Lineage",
    "Coarse_Name", "Fine_Name"
  )

  numeric_cols <- colnames(df)[vapply(df, is.numeric, logical(1))]
  marker_cols  <- setdiff(numeric_cols, ignore_cols)

  summary_list <- Build_Cluster_Summary_Values(df, marker_cols, group_col = "Cluster")

  defs <- data.frame(
    cluster         = names(summary_list),
    major_lineage   = NA_character_,
    minor_lineage   = NA_character_,
    state_name      = NA_character_,
    coarse_name     = NA_character_,
    fine_name       = NA_character_,
    confidence      = NA_real_,
    method          = NA_character_,
    Biomarker       = NA_character_,
    reasoning       = NA_character_,
    stringsAsFactors = FALSE
  )

  #              -       
  if (!use_llm_only && length(rules) > 0) {
    message("\n>>> [Track 1] Running rule-based matching with negative marker support...")

    for (i in seq_along(summary_list)) {
      cl_data <- summary_list[[i]]$expression_profile
      cl_name <- summary_list[[i]]$cluster

      #    Match_Rules_With_Negatives   
      matched_rule <- Match_Rules_With_Negatives(cl_data, rules)

      if (!is.null(matched_rule)) {
        #            
        rule <- rules[[matched_rule]]
        if (is.character(rule)) {
          biomarker_str <- paste(rule, collapse = ", ")
        } else {
          pos_markers <- rule$positive %||% character(0)
          neg_markers <- rule$negative %||% character(0)
          biomarker_str <- paste0(
            "Positive: ", paste(pos_markers, collapse = ", "),
            if (length(neg_markers) > 0) paste0(" | Negative: ", paste(neg_markers, collapse = ", ")) else ""
          )
        }

        norm_res <- Normalize_Annotation_Label(
          major_lineage = matched_rule,
          minor_lineage = "",
          state_name    = "",
          biomarker     = biomarker_str,
          reasoning     = "User defined rules matched (with negative markers)",
          base_population_context = base_population_context,
          allow_functional_prefix = allow_functional_prefix,
          ontology = ontology,
          unknown_threshold = unknown_threshold,
          confidence = 100
        )

        defs$major_lineage[defs$cluster == cl_name] <- norm_res$major_lineage
        defs$minor_lineage[defs$cluster == cl_name] <- norm_res$minor_lineage
        defs$state_name[defs$cluster == cl_name]    <- norm_res$state_name
        defs$coarse_name[defs$cluster == cl_name]   <- norm_res$coarse_name
        defs$fine_name[defs$cluster == cl_name]     <- norm_res$fine_name
        defs$method[defs$cluster == cl_name]        <- "Rule-based"
        defs$confidence[defs$cluster == cl_name]    <- 100
        defs$Biomarker[defs$cluster == cl_name]     <- biomarker_str
        defs$reasoning[defs$cluster == cl_name]     <- "User defined rules matched (with negative markers)"
      }
    }
  } else if (use_llm_only) {
    message("\n>>> [Mode switch] Pure Agentic-LLM mode selected; skipping rule-based hard matching...")
  }

  unannotated_clusters <- defs$cluster[is.na(defs$fine_name)]

  if (length(unannotated_clusters) > 0) {
    message(paste0(
      "\n>>> [Track 2] Detected ", length(unannotated_clusters),
      " clusters requiring intelligent analysis; launching the Agentic LLM workflow..."
    ))

    if (!is.null(ontology)) {
      message("   Mounted a controlled ontology dictionary to constrain AI output.")
    }

    if (!is.null(base_population_context) && nzchar(base_population_context)) {
      message(sprintf("   Injected global base context: [%s]", base_population_context))
    }

    un_list <- summary_list[unannotated_clusters]
    payload <- jsonlite::toJSON(un_list, auto_unbox = TRUE, pretty = TRUE)
    api_key <- Sys.getenv("OPENROUTER_API_KEY")

    parsed_final <- Run_Agentic_LLM_Loop(
      cluster_payload         = payload,
      model                   = model,
      api_key                 = api_key,
      ontology                = ontology,
      base_population_context = base_population_context,
      max_clusters_per_batch  = max_clusters_per_batch,
      allow_functional_prefix = allow_functional_prefix,
      unknown_threshold       = unknown_threshold
    )

    if (!is.null(parsed_final)) {
      for (j in seq_len(nrow(parsed_final))) {
        idx <- defs$cluster == as.character(parsed_final$cluster[j])
        if (!any(idx)) next

        proposed_major <- if ("major_lineage" %in% colnames(parsed_final)) parsed_final$major_lineage[j] else NA
        proposed_minor <- if ("minor_lineage" %in% colnames(parsed_final)) parsed_final$minor_lineage[j] else ""
        proposed_state <- if ("state_name" %in% colnames(parsed_final)) parsed_final$state_name[j] else ""
        proposed_biomk <- if ("Biomarker" %in% colnames(parsed_final)) parsed_final$Biomarker[j] else ""
        proposed_reason <- if ("reasoning" %in% colnames(parsed_final)) parsed_final$reasoning[j] else ""
        proposed_conf <- if ("confidence" %in% colnames(parsed_final)) parsed_final$confidence[j] else NA

        norm_res <- Normalize_Annotation_Label(
          major_lineage = proposed_major,
          minor_lineage = proposed_minor,
          state_name    = proposed_state,
          biomarker     = proposed_biomk,
          reasoning     = proposed_reason,
          base_population_context = base_population_context,
          allow_functional_prefix = allow_functional_prefix,
          ontology = ontology,
          unknown_threshold = unknown_threshold,
          confidence = proposed_conf
        )

        defs$major_lineage[idx] <- norm_res$major_lineage
        defs$minor_lineage[idx] <- norm_res$minor_lineage
        defs$state_name[idx]    <- norm_res$state_name
        defs$coarse_name[idx]   <- norm_res$coarse_name
        defs$fine_name[idx]     <- norm_res$fine_name
        defs$method[idx]        <- "Agentic-LLM"
        defs$reasoning[idx]     <- proposed_reason
        defs$confidence[idx]    <- suppressWarnings(as.numeric(proposed_conf))
        defs$Biomarker[idx]     <- proposed_biomk
      }
    } else {
      warning("   Agentic LLM workflow interrupted; remaining clusters were marked as Unknown.")
      defs$major_lineage[is.na(defs$major_lineage)] <- "Unknown"
      defs$minor_lineage[is.na(defs$minor_lineage)] <- ""
      defs$state_name[is.na(defs$state_name)]       <- ""
      defs$coarse_name[is.na(defs$coarse_name)]     <- "Unknown"
      defs$fine_name[is.na(defs$fine_name)]         <- "Unknown (API Error)"
      defs$confidence[is.na(defs$confidence)]       <- 0
    }
  }

  defs$coarse_name[is.na(defs$coarse_name) | !nzchar(defs$coarse_name)] <- "Unknown"
  defs$fine_name[is.na(defs$fine_name) | !nzchar(defs$fine_name)]       <- defs$coarse_name[is.na(defs$fine_name) | !nzchar(defs$fine_name)]
  defs$major_lineage[is.na(defs$major_lineage) | !nzchar(defs$major_lineage)] <- defs$coarse_name[is.na(defs$major_lineage) | !nzchar(defs$major_lineage)]
  defs$minor_lineage[is.na(defs$minor_lineage)] <- ""
  defs$state_name[is.na(defs$state_name)]       <- ""

  report_file <- "Cluster_Annotation_Report.csv"
  utils::write.csv(defs, report_file, row.names = FALSE)
  message("\n>>> Detailed annotation report exported to: [", report_file, "]")

  coarse_map <- setNames(defs$coarse_name, defs$cluster)
  fine_map   <- setNames(defs$fine_name, defs$cluster)
  major_map  <- setNames(defs$major_lineage, defs$cluster)
  minor_map  <- setNames(defs$minor_lineage, defs$cluster)
  state_map  <- setNames(defs$state_name, defs$cluster)

  res_obj$Data$Coarse_Name   <- as.character(coarse_map[as.character(res_obj$Data$Cluster)])
  res_obj$Data$Fine_Name     <- as.character(fine_map[as.character(res_obj$Data$Cluster)])
  res_obj$Data$Major_Lineage <- as.character(major_map[as.character(res_obj$Data$Cluster)])
  res_obj$Data$Minor_Lineage <- as.character(minor_map[as.character(res_obj$Data$Cluster)])
  res_obj$Data$State_Name    <- as.character(state_map[as.character(res_obj$Data$Cluster)])

  if (annotation_granularity == "coarse") {
    res_obj$Data$Cluster_Name <- res_obj$Data$Coarse_Name
  } else if (annotation_granularity == "fine") {
    res_obj$Data$Cluster_Name <- res_obj$Data$Fine_Name
  } else {
    res_obj$Data$Cluster_Name <- paste0(res_obj$Data$Coarse_Name, " | ", res_obj$Data$Fine_Name)
  }

  res_obj$Data$Lineage_Name <- res_obj$Data$Major_Lineage
  res_obj$Cluster_Definitions <- defs

  # Post-processing validation
  message("\n>>> [Post-processing] Validating annotation quality...")

  low_conf_idx <- which(defs$confidence < unknown_threshold & defs$method == "Agentic-LLM")
  if (length(low_conf_idx) > 0) {
    message(sprintf("       Found %d low-confidence clusters (confidence < %d), marking as 'Unknown'",
                    length(low_conf_idx), unknown_threshold))
    defs$major_lineage[low_conf_idx] <- "Unknown"
    defs$minor_lineage[low_conf_idx] <- ""
    defs$state_name[low_conf_idx] <- ""
    defs$coarse_name[low_conf_idx] <- "Unknown"
    defs$fine_name[low_conf_idx] <- "Unknown (Low Confidence)"

    for (idx in low_conf_idx) {
      cl <- defs$cluster[idx]
      res_obj$Data$Major_Lineage[res_obj$Data$Cluster == cl] <- "Unknown"
      res_obj$Data$Minor_Lineage[res_obj$Data$Cluster == cl] <- ""
      res_obj$Data$State_Name[res_obj$Data$Cluster == cl] <- ""
      res_obj$Data$Coarse_Name[res_obj$Data$Cluster == cl] <- "Unknown"
      res_obj$Data$Fine_Name[res_obj$Data$Cluster == cl] <- "Unknown (Low Confidence)"
      if (annotation_granularity == "coarse") {
        res_obj$Data$Cluster_Name[res_obj$Data$Cluster == cl] <- "Unknown"
      } else if (annotation_granularity == "fine") {
        res_obj$Data$Cluster_Name[res_obj$Data$Cluster == cl] <- "Unknown (Low Confidence)"
      } else {
        res_obj$Data$Cluster_Name[res_obj$Data$Cluster == cl] <- "Unknown | Unknown (Low Confidence)"
      }
    }
  }

  message("\n>>> [Annotation Summary]")
  major_counts <- table(defs$major_lineage)
  message(sprintf("   Total clusters annotated: %d", nrow(defs)))
  message("   Cell type distribution:")
  for (ct in names(sort(major_counts, decreasing = TRUE))) {
    pct <- round(100 * major_counts[ct] / nrow(defs), 1)
    message(sprintf("     - %s: %d clusters (%.1f%%)", ct, major_counts[ct], pct))
  }

  unknown_count <- sum(defs$major_lineage == "Unknown")
  unknown_pct <- round(100 * unknown_count / nrow(defs), 1)
  if (unknown_pct > 30) {
    message(sprintf("       WARNING: %.1f%% of clusters are 'Unknown' (threshold: 30%%)", unknown_pct))
    message("   Suggestions:")
    message("     - Check if marker panel is sufficient for cell type identification")
    message("     - Consider lowering unknown_threshold parameter")
    message("     - Review cluster quality (may need re-clustering with different parameters)")
  } else {
    message(sprintf("     Quality check passed: %.1f%% 'Unknown' clusters (threshold: 30%%)", unknown_pct))
  }

  conf_stats <- summary(defs$confidence[defs$method == "Agentic-LLM"])
  message("\n   Confidence score distribution (LLM-annotated clusters):")
  message(sprintf("     Min: %.1f | Median: %.1f | Mean: %.1f | Max: %.1f",
                  conf_stats["Min."], conf_stats["Median"], conf_stats["Mean"], conf_stats["Max."]))

  message("\n>>> Annotation workflow completed!")
  res_obj
}

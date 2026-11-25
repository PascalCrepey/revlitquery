# R/08_shiny_helpers.R
# Pure, testable helper functions for the Shiny Query Optimizer application.

#' Split a textarea string of terms on commas, semicolons, or newlines
#'
#' @param txt Character scalar.
#' @return Cleaned lowercase unique character vector.
parse_term_list <- function(txt) {
  if (is.null(txt) || length(txt) == 0L || is.na(txt)) return(character(0))
  raw <- as.character(txt[1L])
  if (!nzchar(stringr::str_trim(raw))) return(character(0))
  parts <- unlist(stringr::str_split(raw, "[,;\r\n]+"))
  parts <- tolower(stringr::str_squish(parts))
  unique(parts[nzchar(parts)])
}

#' Parse optional custom concept blocks from YAML text
#'
#' @param txt YAML text representing a named mapping of term lists.
#' @return List with `ok` (logical), `blocks` (named list or NULL), `error` (character or NULL).
parse_custom_blocks_yaml <- function(txt) {
  if (is.null(txt) || length(txt) == 0L || is.na(txt)) {
    return(list(ok = TRUE, blocks = NULL, error = NULL))
  }
  raw <- stringr::str_trim(as.character(txt[1L]))
  if (!nzchar(raw)) {
    return(list(ok = TRUE, blocks = NULL, error = NULL))
  }

  parsed <- tryCatch(
    yaml::yaml.load(raw),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "yaml_err")
  )
  if (inherits(parsed, "yaml_err")) {
    return(list(ok = FALSE, blocks = NULL, error = paste("Invalid YAML in custom blocks:", parsed$msg)))
  }
  if (!is.list(parsed) || is.null(names(parsed)) || any(!nzchar(stringr::str_trim(names(parsed))))) {
    return(list(
      ok = FALSE,
      blocks = NULL,
      error = "Custom blocks YAML must be a named mapping (e.g. `pathogen:\n  - h5n1`)."
    ))
  }

  cleaned <- list()
  for (nm in names(parsed)) {
    vals <- parsed[[nm]]
    if (is.null(vals) || is.list(vals) && length(vals) == 0L) {
      return(list(ok = FALSE, blocks = NULL, error = paste0("Block '", nm, "' cannot be empty.")))
    }
    vec <- tolower(stringr::str_squish(as.character(unlist(vals, use.names = FALSE))))
    vec <- unique(vec[nzchar(vec)])
    if (length(vec) == 0L) {
      return(list(ok = FALSE, blocks = NULL, error = paste0("Block '", nm, "' has no valid non-empty terms.")))
    }
    cleaned[[stringr::str_squish(nm)]] <- vec
  }

  list(ok = TRUE, blocks = cleaned, error = NULL)
}

#' Convert discovered concept blocks list into editable YAML text
#'
#' @param blocks List of `list(name = ..., terms = ...)`.
#' @return Character scalar YAML string.
blocks_to_yaml <- function(blocks) {
  if (is.null(blocks) || length(blocks) == 0L) return("")
  named <- stats::setNames(
    lapply(blocks, function(b) as.list(as.character(b$terms))),
    vapply(seq_along(blocks), function(i) blocks[[i]]$name %||% paste0("concept_", i), character(1))
  )
  yaml::as.yaml(named)
}

#' Validate and assemble optimizer config list from Shiny sidebar inputs
#'
#' @param required_txt Comma/newline-separated required terms string.
#' @param banned_txt Comma/newline-separated banned terms string.
#' @param exclude_txt Comma/newline-separated exclusion terms string.
#' @param custom_blocks_yaml Optional YAML text defining custom concept blocks.
#' @param max_blocks Maximum number of concept blocks.
#' @param max_terms_per_block Maximum number of terms per concept block.
#' @param max_gs_chars Maximum Google Scholar query character limit.
#' @return List with `ok`, `config`, and `error`.
build_config_from_inputs <- function(required_txt = "",
                                     banned_txt = "",
                                     exclude_txt = "",
                                     custom_blocks_yaml = "",
                                     max_blocks = 3L,
                                     max_terms_per_block = 5L,
                                     max_gs_chars = 256L) {
  cb_res <- parse_custom_blocks_yaml(custom_blocks_yaml)
  if (!isTRUE(cb_res$ok)) {
    return(list(ok = FALSE, config = NULL, error = cb_res$error))
  }

  cfg <- list(
    required_terms = parse_term_list(required_txt),
    banned_terms = parse_term_list(banned_txt),
    exclude_terms = parse_term_list(exclude_txt),
    max_blocks = as.integer(max_blocks %||% 3L),
    max_terms_per_block = as.integer(max_terms_per_block %||% 5L),
    max_gs_chars = as.integer(max_gs_chars %||% 256L)
  )
  if (!is.null(cb_res$blocks)) {
    cfg$custom_blocks <- cb_res$blocks
  }

  list(ok = TRUE, config = cfg, error = NULL)
}

#' Recompute recommended query on an existing Pareto frontier when sliders move
#'
#' @param pareto_df Pareto frontier tibble.
#' @param min_recall Minimum seed recall in `0..1`.
#' @param max_hits Maximum PubMed hit threshold (`Inf` or numeric).
#' @return List with updated `pareto_df`, `recommended_sig`, and `recommended_status`.
recompute_recommendation <- function(pareto_df, min_recall = 0.60, max_hits = Inf) {
  if (is.null(pareto_df) || nrow(pareto_df) == 0L) {
    return(list(
      pareto_df = pareto_df,
      recommended_sig = NA_character_,
      recommended_status = "empty"
    ))
  }

  if (is.null(max_hits) || length(max_hits) == 0L || is.na(max_hits) || max_hits <= 0) {
    max_hits <- Inf
  }

  meets_both <- pareto_df |>
    dplyr::filter(.data$combined_recall >= min_recall & .data$pm_hits <= max_hits)

  if (nrow(meets_both) > 0L) {
    rec_sig <- meets_both$signature[1L]
    rec_status <- "meets_constraints"
  } else {
    meets_rec <- pareto_df |> dplyr::filter(.data$combined_recall >= min_recall)
    if (nrow(meets_rec) > 0L) {
      rec_sig <- meets_rec$signature[which.min(meets_rec$pm_hits)]
      rec_status <- "exceeds_max_hits"
    } else {
      rec_sig <- pareto_df$signature[1L]
      rec_status <- "below_min_recall"
    }
  }

  pareto_df <- pareto_df |>
    dplyr::mutate(is_recommended = .data$signature == rec_sig)

  list(
    pareto_df = pareto_df,
    recommended_sig = rec_sig,
    recommended_status = rec_status
  )
}

#' Build interactive Plotly Pareto frontier chart
#'
#' @param eval_df All evaluated queries tibble.
#' @param pareto_df Pareto-optimal queries tibble.
#' @param selected_rank Currently inspected Pareto rank (`Q1`, `Q2`, ...).
#' @param min_recall Target minimum recall (`0..1`).
#' @param max_hits Optional max PubMed hit threshold (`Inf` if none).
#' @param selected_eval_id Currently inspected `eval_id`.
#' @return A `plotly` htmlwidget.
build_pareto_plotly <- function(eval_df,
                                pareto_df,
                                selected_rank = NULL,
                                min_recall = 0.60,
                                max_hits = Inf,
                                selected_eval_id = NULL) {
  valid_eval <- if (!is.null(eval_df) && nrow(eval_df) > 0L) {
    if (!("is_custom" %in% names(eval_df))) eval_df$is_custom <- FALSE
    if (!("custom_label" %in% names(eval_df))) eval_df$custom_label <- ""
    eval_df |> dplyr::filter(!is.na(.data$pm_hits), !is.na(.data$combined_recall))
  } else {
    tibble::tibble()
  }

  if (nrow(valid_eval) == 0L) {
    p_empty <- plotly::plot_ly(
      x = numeric(0), y = numeric(0),
      type = "scatter", mode = "markers", source = "pareto"
    ) |>
      plotly::layout(
        xaxis = list(visible = FALSE),
        yaxis = list(visible = FALSE),
        annotations = list(list(
          text = "Run the optimizer or add a custom query to visualize the sensitivity vs. specificity Pareto frontier.",
          xref = "paper", yref = "paper", x = 0.5, y = 0.5, showarrow = FALSE,
          font = list(size = 14, color = "#6c757d")
        ))
      ) |>
      plotly::config(displaylogo = FALSE)
    return(plotly::event_register(p_empty, "plotly_click"))
  }

  dom_df <- valid_eval |>
    dplyr::filter(!.data$is_pareto & !.data$is_custom) |>
    dplyr::mutate(
      x_plot = pmax(as.numeric(.data$pm_hits), 1),
      y_pct = 100 * .data$combined_recall,
      hover = sprintf(
        "<b>Evaluated #%d</b> (%s)<br>Recall: %.1f%% (%d/%d)<br>PubMed hits: %s<br>GS chars: %d<br><code>%s</code>",
        .data$eval_id, .data$stage, .data$y_pct, .data$seeds_found, .data$seeds_resolvable,
        format(.data$pm_hits, big.mark = ","), .data$gs_chars,
        htmltools::htmlEscape(stringr::str_trunc(.data$google_scholar_query, 90))
      )
    )

  cust_df <- valid_eval |>
    dplyr::filter(.data$is_custom) |>
    dplyr::mutate(
      x_plot = pmax(as.numeric(.data$pm_hits), 1),
      y_pct = 100 * .data$combined_recall,
      is_sel = if (!is.null(selected_eval_id) && !is.na(selected_eval_id)) .data$eval_id == selected_eval_id else FALSE,
      lbl = ifelse(nzchar(.data$custom_label), .data$custom_label, paste0("Custom #", .data$eval_id)),
      hover = sprintf(
        "<b>%s</b> (Custom%s)%s<br>Recall: %.1f%% (%d/%d)<br>PubMed hits: %s<br>GS chars: %d<br><code>%s</code>",
        htmltools::htmlEscape(.data$lbl),
        ifelse(.data$is_pareto, ", Pareto", ""),
        ifelse(.data$is_sel, " [Selected]", ""),
        .data$y_pct, .data$seeds_found, .data$seeds_resolvable,
        format(.data$pm_hits, big.mark = ","), .data$gs_chars,
        htmltools::htmlEscape(stringr::str_trunc(.data$google_scholar_query, 90))
      )
    )

  p_df <- if (!is.null(pareto_df) && nrow(pareto_df) > 0L) {
    pareto_df |>
      dplyr::filter(!is.na(.data$pm_hits), !is.na(.data$combined_recall)) |>
      dplyr::arrange(.data$pm_hits, dplyr::desc(.data$combined_recall)) |>
      dplyr::mutate(
        x_plot = pmax(as.numeric(.data$pm_hits), 1),
        y_pct = 100 * .data$combined_recall,
        is_sel = if (!is.null(selected_eval_id) && !is.na(selected_eval_id)) {
          .data$eval_id == selected_eval_id
        } else if (!is.null(selected_rank) && !is.na(selected_rank)) {
          .data$rank == selected_rank
        } else {
          FALSE
        },
        pt_color = dplyr::case_when(
          .data$is_recommended ~ "#d93025",
          TRUE ~ "#1f6feb"
        ),
        pt_size = ifelse(.data$is_sel, 15, ifelse(.data$is_recommended, 12, 10)),
        pt_line_color = ifelse(.data$is_sel, "#111827", "#ffffff"),
        pt_line_width = ifelse(.data$is_sel, 2.5, 1.5),
        lbl = paste0("Q", .data$rank, ifelse(.data$is_recommended, "*", "")),
        hover = sprintf(
          "<b>Q%d%s</b>%s<br>Recall: %.1f%% (%d/%d)<br>PubMed hits: %s<br>GS chars: %d<br><code>%s</code>",
          .data$rank,
          ifelse(.data$is_recommended, " (Recommended)", ""),
          ifelse(.data$is_sel, " [Selected]", ""),
          .data$y_pct, .data$seeds_found, .data$seeds_resolvable,
          format(.data$pm_hits, big.mark = ","), .data$gs_chars,
          htmltools::htmlEscape(stringr::str_trunc(.data$google_scholar_query, 90))
        )
      )
  } else {
    tibble::tibble()
  }

  p <- plotly::plot_ly(source = "pareto")

  if (nrow(dom_df) > 0L) {
    p <- p |>
      plotly::add_markers(
        data = dom_df,
        x = ~x_plot, y = ~y_pct,
        customdata = ~eval_id,
        text = ~hover, hoverinfo = "text",
        name = "Evaluated (dominated)",
        marker = list(color = "rgba(150, 160, 175, 0.55)", size = 7)
      )
  }

  if (nrow(p_df) > 0L) {
    p <- p |>
      plotly::add_trace(
        data = p_df,
        x = ~x_plot, y = ~y_pct,
        type = "scatter", mode = "lines+markers+text",
        customdata = ~eval_id,
        text = ~lbl, textposition = "top center",
        textfont = list(size = 11, color = "#1f2937"),
        hovertext = ~hover, hoverinfo = "text",
        name = "Pareto optimal",
        line = list(color = "#1f6feb", width = 2),
        marker = list(
          color = p_df$pt_color,
          size = p_df$pt_size,
          line = list(color = p_df$pt_line_color, width = p_df$pt_line_width)
        )
      )
  }

  if (nrow(cust_df) > 0L) {
    p <- p |>
      plotly::add_trace(
        data = cust_df,
        x = ~x_plot, y = ~y_pct,
        type = "scatter", mode = "markers+text",
        customdata = ~eval_id,
        text = ~lbl, textposition = "bottom center",
        textfont = list(size = 11, color = "#6b21a8"),
        hovertext = ~hover, hoverinfo = "text",
        name = "Custom query",
        marker = list(
          symbol = "diamond",
          color = "#9333ea",
          size = ifelse(cust_df$is_sel, 15, 12),
          line = list(color = ifelse(cust_df$is_sel, "#111827", "#ffffff"), width = ifelse(cust_df$is_sel, 2.5, 1.5))
        )
      )
  }

  shapes <- list(
    list(
      type = "line", xref = "paper", yref = "y",
      x0 = 0, x1 = 1, y0 = 100 * min_recall, y1 = 100 * min_recall,
      line = list(color = "#6b7280", width = 1.2, dash = "dash")
    )
  )
  if (!is.null(max_hits) && !is.na(max_hits) && is.finite(max_hits) && max_hits > 0) {
    shapes[[length(shapes) + 1L]] <- list(
      type = "line", xref = "x", yref = "paper",
      x0 = max(as.numeric(max_hits), 1), x1 = max(as.numeric(max_hits), 1),
      y0 = 0, y1 = 1,
      line = list(color = "#d93025", width = 1.2, dash = "dot")
    )
  }

  p <- p |>
    plotly::layout(
      xaxis = list(title = "Total PubMed Hits (log10 scale)", type = "log"),
      yaxis = list(title = "Seed Recall (%)", range = c(-3, 108)),
      shapes = shapes,
      legend = list(orientation = "h", x = 0, y = 1.12),
      margin = list(t = 40, r = 20, b = 50, l = 55),
      hovermode = "closest"
    ) |>
    plotly::config(displaylogo = FALSE)

  plotly::event_register(p, "plotly_click")
}

#' Build display table for the Seed Coverage Matrix tab
#'
#' @param seeds_df Enriched seeds tibble.
#' @param pareto_df Pareto frontier tibble.
#' @param eval_df Optional all evaluated queries tibble (for custom queries).
#' @return Tibble ready for `DT::datatable()`.
build_coverage_display_df <- function(seeds_df, pareto_df, eval_df = NULL) {
  if (is.null(seeds_df) || nrow(seeds_df) == 0L) return(tibble::tibble())

  resolvable <- as.logical(seeds_df$in_pubmed %||% FALSE)
  out <- tibble::tibble(
    Citekey = seeds_df$citekey,
    Year = seeds_df$year,
    Title = seeds_df$full_title %||% seeds_df$title,
    PubMed = ifelse(seeds_df$in_pubmed, seeds_df$pmid, "Not found")
  )

  if (!is.null(pareto_df) && nrow(pareto_df) > 0L) {
    for (r in seq_len(nrow(pareto_df))) {
      q_col <- paste0("Q", pareto_df$rank[r])
      found_i <- pareto_df$found_idx[[r]]
      out[[q_col]] <- ifelse(
        !resolvable,
        "N/A",
        ifelse(seeds_df$seed_id %in% found_i, "Found", "Missed")
      )
    }
  }

  if (!is.null(eval_df) && nrow(eval_df) > 0L && ("is_custom" %in% names(eval_df))) {
    extra_c <- eval_df |> dplyr::filter(.data$is_custom & !.data$is_pareto)
    for (r in seq_len(nrow(extra_c))) {
      lbl <- if ("custom_label" %in% names(extra_c) && nzchar(extra_c$custom_label[r])) {
        extra_c$custom_label[r]
      } else {
        paste0("#", extra_c$eval_id[r])
      }
      c_col <- make.unique(c(names(out), paste0("Custom: ", lbl)))[length(names(out)) + 1L]
      found_i <- extra_c$found_idx[[r]]
      out[[c_col]] <- ifelse(
        !resolvable,
        "N/A",
        ifelse(seeds_df$seed_id %in% found_i, "Found", "Missed")
      )
    }
  }
  out
}

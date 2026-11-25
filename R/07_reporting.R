# R/07_reporting.R
# Generate CSV outputs, seed coverage matrix, Markdown report, and Pareto frontier plot.

#' Write all optimization artifacts (CSVs, Markdown report, PNG plot) to `output_dir`
#'
#' @param result Optimization result list from `optimize_queries()`.
#' @param output_dir Output directory path.
#' @return Invisibly returns list of generated file paths.
#' @export
write_optimization_report <- function(result, output_dir = "output") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pareto_df <- result$pareto_df
  eval_df <- result$eval_df
  seeds_df <- result$seeds_df
  if (!("is_custom" %in% names(eval_df))) eval_df$is_custom <- FALSE
  if (!("custom_label" %in% names(eval_df))) eval_df$custom_label <- ""
  if (!("is_custom" %in% names(pareto_df))) pareto_df$is_custom <- FALSE
  if (!("custom_label" %in% names(pareto_df))) pareto_df$custom_label <- ""

  # Drop list-columns for flat CSV export
  export_pareto <- pareto_df |>
    dplyr::select(
      "rank", "is_recommended", "stage", "is_custom", "custom_label",
      "combined_recall", "seeds_found", "seeds_resolvable",
      "pm_recall", "pm_hits",
      "gs_chars", "n_blocks", "n_terms",
      "pubmed_query", "google_scholar_query",
      "missed_citekeys", "pubmed_url", "google_scholar_url"
    )

  export_eval <- eval_df |>
    dplyr::arrange(dplyr::desc(.data$is_pareto), dplyr::desc(.data$combined_recall), .data$pm_hits) |>
    dplyr::select(
      "eval_id", "is_pareto", "stage", "is_custom", "custom_label",
      "combined_recall", "seeds_found", "seeds_resolvable",
      "pm_recall", "pm_hits",
      "gs_chars", "n_blocks", "n_terms",
      "pubmed_query", "google_scholar_query",
      "missed_citekeys"
    )

  seed_res <- seeds_df |>
    dplyr::mutate(has_abstract = nzchar(.data$full_abstract)) |>
    dplyr::select(
      "seed_id", "citekey", "year", "full_title",
      "doi", "pmid", "in_pubmed", "has_abstract"
    )

  # Seed coverage matrix: rows = seeds, cols = Q1..Qm + any non-Pareto custom queries
  cov_mat <- tibble::tibble(
    seed_id = seeds_df$seed_id,
    citekey = seeds_df$citekey,
    in_pubmed = seeds_df$in_pubmed,
    title = seeds_df$full_title
  )
  for (r in seq_len(nrow(pareto_df))) {
    q_label <- paste0("Q", pareto_df$rank[r])
    found_i <- pareto_df$found_idx[[r]]
    cov_mat[[q_label]] <- seeds_df$seed_id %in% found_i
  }
  extra_custom <- eval_df |>
    dplyr::filter(.data$is_custom & !.data$is_pareto)
  for (r in seq_len(nrow(extra_custom))) {
    lbl_raw <- if (nzchar(extra_custom$custom_label[r])) extra_custom$custom_label[r] else paste0("#", extra_custom$eval_id[r])
    c_col <- make.unique(c(names(cov_mat), paste0("Custom: ", lbl_raw)))[length(names(cov_mat)) + 1L]
    cov_mat[[c_col]] <- seeds_df$seed_id %in% extra_custom$found_idx[[r]]
  }

  p_pareto_csv <- file.path(output_dir, "pareto_queries.csv")
  p_eval_csv <- file.path(output_dir, "all_evaluated_queries.csv")
  p_seed_csv <- file.path(output_dir, "seed_resolution.csv")
  p_cov_csv <- file.path(output_dir, "seed_coverage_matrix.csv")
  p_excl_csv <- file.path(output_dir, "exclusion_hints.csv")
  p_report_md <- file.path(output_dir, "report.md")
  p_plot_png <- file.path(output_dir, "pareto_frontier.png")

  excl_df <- result$exclusion_hints %||% .empty_exclusion_hints()
  readr::write_csv(export_pareto, p_pareto_csv)
  readr::write_csv(export_eval, p_eval_csv)
  readr::write_csv(seed_res, p_seed_csv)
  readr::write_csv(cov_mat, p_cov_csv)
  readr::write_csv(excl_df, p_excl_csv)

  # Build Pareto plot
  plot_df <- eval_df |>
    dplyr::filter(!is.na(.data$pm_hits), !is.na(.data$combined_recall)) |>
    dplyr::mutate(
      pm_hits_plot = pmax(as.numeric(.data$pm_hits), 1),
      recall_pct = 100 * .data$combined_recall
    )
  p_pts <- pareto_df |>
    dplyr::filter(!is.na(.data$pm_hits), !is.na(.data$combined_recall)) |>
    dplyr::mutate(
      pm_hits_plot = pmax(as.numeric(.data$pm_hits), 1),
      recall_pct = 100 * .data$combined_recall,
      rec_fac = factor(as.character(.data$is_recommended), levels = c("FALSE", "TRUE")),
      label = paste0("Q", .data$rank, ifelse(.data$is_recommended, "*", ""))
    )

  c_pts <- plot_df |>
    dplyr::filter(.data$is_custom) |>
    dplyr::mutate(
      c_label = ifelse(nzchar(.data$custom_label), .data$custom_label, paste0("Custom #", .data$eval_id))
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$pm_hits_plot, y = .data$recall_pct)) +
    ggplot2::geom_point(
      data = dplyr::filter(plot_df, !.data$is_pareto & !.data$is_custom),
      color = "grey75", alpha = 0.7, size = 2
    )
  if (nrow(p_pts) >= 2L) {
    p <- p + ggplot2::geom_line(
      data = dplyr::arrange(p_pts, .data$pm_hits_plot),
      color = "#1f6feb", linewidth = 0.9
    )
  }
  p <- p +
    ggplot2::geom_point(
      data = p_pts,
      ggplot2::aes(color = .data$rec_fac),
      size = 3.5
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "#1f6feb", "TRUE" = "#d93025"),
      labels = c("FALSE" = "Pareto optimal", "TRUE" = "Recommended"),
      drop = FALSE,
      name = NULL
    ) +
    ggrepel::geom_text_repel(
      data = p_pts,
      ggplot2::aes(label = .data$label),
      size = 3.4, fontface = "bold", max.overlaps = 20
    )

  if (nrow(c_pts) > 0L) {
    p <- p +
      ggplot2::geom_point(data = c_pts, color = "#9333ea", shape = 18, size = 4.2) +
      ggrepel::geom_text_repel(
        data = c_pts,
        ggplot2::aes(label = .data$c_label),
        color = "#6b21a8", size = 3.2, fontface = "italic", max.overlaps = 15
      )
  }

  p <- p +
    ggplot2::scale_x_log10(labels = scales::label_comma()) +
    ggplot2::scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
    ggplot2::geom_hline(yintercept = 100 * result$min_recall, linetype = "dashed", color = "grey45") +
    ggplot2::labs(
      title = "PubMed & Google Scholar Query Pareto Frontier",
      subtitle = glue::glue(
        "Evaluated {nrow(eval_df)} queries across {nrow(seeds_df)} BibTeX seed papers (* = recommended",
        "{ifelse(is.finite(result$max_hits), paste0(', max_hits <= ', format(result$max_hits, big.mark = ',')), '')})"
      ),
      x = "Total PubMed Hits (log10 scale)",
      y = "Seed Recall (%)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  if (is.finite(result$max_hits) && result$max_hits > 0) {
    p <- p + ggplot2::geom_vline(xintercept = max(result$max_hits, 1), linetype = "dotted", color = "#d93025")
  }

  ggplot2::ggsave(p_plot_png, plot = p, width = 8.5, height = 5.5, dpi = 150)

  # Build Markdown report
  n_total <- nrow(seeds_df)
  n_pm <- sum(seeds_df$in_pubmed)
  n_abs <- sum(nzchar(seeds_df$full_abstract))

  rec_row <- dplyr::filter(pareto_df, .data$is_recommended)
  if (nrow(rec_row) == 0L) rec_row <- pareto_df[1L, , drop = FALSE]

  status_note <- switch(result$recommended_status,
    meets_constraints = glue::glue("Meets both `min_recall >= {round(100 * result$min_recall)}%` and `max_hits <= {result$max_hits}`."),
    exceeds_max_hits = glue::glue("No Pareto query stayed under `max_hits = {result$max_hits}` while meeting `min_recall >= {round(100 * result$min_recall)}%`; showing lowest-hit query meeting `min_recall`."),
    below_min_recall = glue::glue("No evaluated query reached `min_recall >= {round(100 * result$min_recall)}%`; showing highest-recall Pareto query.")
  )

  # Summary table rows
  tbl_lines <- vapply(seq_len(nrow(pareto_df)), function(i) {
    r <- pareto_df[i, ]
    rec_tag <- if (isTRUE(r$is_recommended)) " **\\***" else ""
    missed <- if (nzchar(r$missed_citekeys)) r$missed_citekeys else "None"
    glue::glue(
      "| Q{r$rank}{rec_tag} | {round(100 * r$combined_recall, 1)}% ({r$seeds_found}/{r$seeds_resolvable}) | {format(r$pm_hits, big.mark = ',')} | {r$gs_chars} | {missed} |"
    )
  }, character(1))

  # Detailed blocks for each Pareto query
  detail_blocks <- vapply(seq_len(nrow(pareto_df)), function(i) {
    r <- pareto_df[i, ]
    rec_badge <- if (isTRUE(r$is_recommended)) " (Recommended)" else ""
    glue::glue(
      "### Q{r$rank}{rec_badge}: {round(100 * r$combined_recall, 1)}% recall | {format(r$pm_hits, big.mark = ',')} PubMed hits | {r$gs_chars} GS chars\n\n",
      "- **Missed seeds**: {ifelse(nzchar(r$missed_citekeys), r$missed_citekeys, 'None')}\n",
      "- **Verify live**: [PubMed]({r$pubmed_url}) | [Google Scholar]({r$google_scholar_url})\n\n",
      "**PubMed query:**\n```text\n{r$pubmed_query}\n```\n\n",
      "**Google Scholar query ({r$gs_chars}/{result$max_gs_chars} chars):**\n```text\n{r$google_scholar_query}\n```\n"
    )
  }, character(1))

  excl_section <- if (nrow(excl_df) > 0L) {
    if (!("concept" %in% names(excl_df))) excl_df$concept <- "general_noise"
    excl_blocks <- attr(excl_df, "exclusion_blocks") %||% list()
    blk_lines <- if (length(excl_blocks) > 0L) {
      c(
        "### Exclusion Concept Blocks (Greedy Set-Cover, 0 Seeds Lost)",
        "",
        vapply(excl_blocks, function(eb) {
          glue::glue("- **{eb$name}** (~{eb$noise_covered_pct}% of sampled noise): `{paste(eb$terms, collapse = '`, `')}`")
        }, character(1)),
        ""
      )
    } else {
      character(0)
    }
    top_eh <- utils::head(excl_df, 25L)
    eh_rows <- vapply(seq_len(nrow(top_eh)), function(i) {
      h <- top_eh[i, ]
      lost_txt <- if (nzchar(h$seeds_lost_citekeys)) paste0(h$seeds_lost_n, " (", h$seeds_lost_citekeys, ")") else "0"
      glue::glue(
        "| {h$concept} | `{h$exclude_spec}` | {h$field} | {format(h$exact_hits_removed, big.mark = ',')} | {h$pct_hits_reduced}% | {lost_txt} | {h$safety_tier} | `{h$pm_clause}` | `{h$gs_clause}` |"
      )
    }, character(1))
    safe_specs <- attr(excl_df, "setcover_ex") %||% excl_df$exclude_spec[excl_df$seeds_lost_n == 0L]
    yml_block <- if (length(safe_specs) > 0L) {
      c("", "Copy-pasteable set-cover `exclude_terms` for `config.yml` (0 seeds lost):", "```yaml", "exclude_terms:", paste0("  - ", utils::head(safe_specs, 12L)), "```")
    } else {
      character(0)
    }
    c(
      "## Specificity Hints: Recommended Terms to Exclude",
      "",
      "Terms and prefix wildcards mined from non-seed PubMed results that reduce total hits when excluded (`NOT` / `-`):",
      "",
      blk_lines,
      "| Concept | Exclude Spec | Scope | Hits Removed | % Reduction | Seeds Lost | Safety | PubMed | Google Scholar |",
      "|---|---|---|---:|---:|---:|---|---|---|",
      eh_rows,
      yml_block,
      ""
    )
  } else {
    c("## Specificity Hints: Recommended Terms to Exclude", "", "No significant exclusion terms discovered.", "")
  }

  custom_rows <- dplyr::filter(eval_df, .data$is_custom)
  custom_section <- if (nrow(custom_rows) > 0L) {
    c_blocks <- vapply(seq_len(nrow(custom_rows)), function(i) {
      r <- custom_rows[i, ]
      lbl <- if (nzchar(r$custom_label)) r$custom_label else paste0("Custom #", r$eval_id)
      pareto_tag <- if (isTRUE(r$is_pareto)) " (Pareto-optimal)" else " (Dominated / Over cap)"
      missed <- if (nzchar(r$missed_citekeys)) r$missed_citekeys else "None"
      glue::glue(
        "### {lbl}{pareto_tag}: {round(100 * r$combined_recall, 1)}% recall ({r$seeds_found}/{r$seeds_resolvable}) | {format(r$pm_hits, big.mark = ',')} PubMed hits\n\n",
        "- **Missed seeds**: {missed}\n",
        "- **Verify live**: [PubMed]({r$pubmed_url}) | [Google Scholar]({r$google_scholar_url})\n\n",
        "**PubMed query:**\n```text\n{r$pubmed_query}\n```\n\n",
        "**Google Scholar (approx, {r$gs_chars} chars):**\n```text\n{r$google_scholar_query}\n```\n"
      )
    }, character(1))
    c("## Custom / Benchmark Queries", "", c_blocks)
  } else {
    character(0)
  }

  md_lines <- c(
    "# Search Query Optimization Report",
    "",
    "## Corpus Summary",
    glue::glue("- **Total BibTeX seeds**: {n_total}"),
    glue::glue("- **Resolved in PubMed**: {n_pm}/{n_total}"),
    glue::glue("- **Seeds with full abstract**: {n_abs}/{n_total}"),
    glue::glue("- **Total queries evaluated**: {nrow(eval_df)} ({nrow(pareto_df)} Pareto-optimal)"),
    "",
    "## Recommended Query Pair (Q" %+% rec_row$rank[1L] %+% ")",
    "",
    status_note,
    "",
    glue::glue("- **Seed recall**: {round(100 * rec_row$combined_recall[1L], 1)}% ({rec_row$seeds_found[1L]}/{rec_row$seeds_resolvable[1L]})"),
    glue::glue("- **PubMed total hits**: {format(rec_row$pm_hits[1L], big.mark = ',')}"),
    glue::glue("- **Google Scholar length**: {rec_row$gs_chars[1L]} / {result$max_gs_chars} chars"),
    glue::glue("- **Links**: [Open in PubMed]({rec_row$pubmed_url[1L]}) | [Open in Google Scholar]({rec_row$google_scholar_url[1L]})"),
    "",
    "**PubMed:**",
    "```text",
    rec_row$pubmed_query[1L],
    "```",
    "",
    "**Google Scholar:**",
    "```text",
    rec_row$google_scholar_query[1L],
    "```",
    "",
    "## Pareto Frontier Overview",
    "",
    "| Query | Seed Recall | PubMed Hits | GS Chars | Missed Seeds |",
    "|---|---:|---:|---:|---|",
    tbl_lines,
    "",
    "![Pareto Frontier](pareto_frontier.png)",
    "",
    custom_section,
    excl_section,
    "## All Pareto-Optimal Query Pairs",
    "",
    detail_blocks,
    "## Tuning Specificity vs. Sensitivity",
    "- Lower `--min-recall` (e.g. `0.80`) or set `--max-hits` (e.g. `500`) to favor higher specificity (fewer extraneous papers).",
    "- Use a YAML config file (`--config config.yml`) to lock `required_terms`, exclude `banned_terms`, add `exclude_terms` (`NOT`), or define `custom_blocks`."
  )

  writeLines(md_lines, p_report_md, useBytes = TRUE)

  invisible(list(
    pareto_csv = p_pareto_csv,
    all_eval_csv = p_eval_csv,
    seed_csv = p_seed_csv,
    coverage_csv = p_cov_csv,
    exclusion_csv = p_excl_csv,
    report_md = p_report_md,
    plot_png = p_plot_png
  ))
}

`%+%` <- function(a, b) paste0(a, b)

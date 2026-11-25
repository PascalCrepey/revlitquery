# R/10_package_api.R
# High-level package entry points for revlitquery:
# - optimize_search_queries(): full end-to-end BibTeX -> PubMed & Google Scholar Pareto optimizer
# - revlitquery_app() / run_revlitquery_app(): interactive Shiny application

#' Write a built-in 6-entry synthetic H5N1 demo BibTeX file to a temporary path
#'
#' Used by the interactive Shiny app's "Load sample H5N1 .bib" button and unit tests
#' so the package does not bundle any external `.bib` files or cached data.
#'
#' @param dest_path Destination file path (defaults to `file.path(tempdir(), "sample_h5n1_cfr.bib")`).
#' @return Character path to the written `.bib` file.
#' @keywords internal
.write_demo_bibtex <- function(dest_path = file.path(tempdir(), "sample_h5n1_cfr.bib")) {
  bib_lines <- c(
    "@article{li2008estimates,",
    "  title = {Estimates of the case fatality rate of human {H5N1} avian influenza infections},",
    "  author = {Li, Fu-Chun and Choi, Bernard C. K. and Sly, Timothy and Pak, Anthony W. P.},",
    "  journal = {Emerging Infectious Diseases},",
    "  year = {2008},",
    "  pmid = {18439345},",
    "  abstract = {Human infections with avian influenza A (H5N1) virus have a high reported case fatality rate. However, mild or asymptomatic cases may be underreported, affecting estimates of the true case fatality risk.}",
    "}",
    "",
    "@article{wang2012serological,",
    "  title = {Serological evidence of human infection with avian influenza {A (H5N1)} virus: a systematic review and meta-analysis},",
    "  author = {Wang, Taia T. and Parides, Michael K. and Palese, Peter},",
    "  journal = {Science},",
    "  year = {2012},",
    "  pmid = {22362880},",
    "  abstract = {The World Health Organization reports a case fatality rate of nearly 60 percent for confirmed human H5N1 infections, yet seroprevalence studies suggest many subclinical infections occur.}",
    "}",
    "",
    "@article{fiebig2011epidemiological,",
    "  title = {Epidemiological and clinical characteristics of human avian influenza {A (H5N1)} cases reported globally},",
    "  author = {Fiebig, Lena and Soyka, Jana and Buda, Silke and Buchholz, Udo and Dehnert, Manuel and Haas, Walter},",
    "  journal = {Eurosurveillance},",
    "  year = {2011},",
    "  pmid = {21679679},",
    "  abstract = {Analysis of confirmed human H5N1 cases shows an overall case-fatality ratio above 50 percent, with substantial variation by age, country, and time to hospitalization.}",
    "}",
    "",
    "@article{kandun2008factors,",
    "  title = {Factors associated with case fatality of human {H5N1} virus infections in {Indonesia}: a case series},",
    "  author = {Kandun, I. Nyoman and Tresnaningsih, Erna and Purba, Wilfried H. and Lee, Vernon and Samaan, Gina and Harun, Syahrial and Soni, Erni and Septiawati, Chita and Setiawati, Tetty and Sariwati, Endang and Wandra, Toni},",
    "  journal = {The Lancet},",
    "  year = {2008},",
    "  pmid = {18959499},",
    "  abstract = {Human influenza A (H5N1) infections in Indonesia have had a high case fatality rate. Early oseltamivir treatment was associated with reduced mortality among confirmed H5N1 patients.}",
    "}",
    "",
    "@article{lai2016global,",
    "  title = {Global epidemiology of avian influenza {A H5N1} virus infection in humans, 1997--2015: a systematic review of individual case data},",
    "  author = {Lai, Shengjie and Qin, Ying and Cowling, Benjamin J. and Ren, Xiang and Wardrop, Nicola A. and Gilbert, Marius and Tsang, Tim K. and Wu, Peng and Feng, Luzhao and Jiang, Hui and Peng, Zhibin and Zheng, Jiandong and Liao, Qiaohong and Li, Sa and Horby, Peter W. and Farrar, Jeremy J. and Gao, George F. and Tatem, Andrew J. and Yu, Hongjie},",
    "  journal = {The Lancet Infectious Diseases},",
    "  year = {2016},",
    "  pmid = {27211899},",
    "  abstract = {Between 1997 and 2015, 907 human H5N1 cases were documented worldwide with an overall case fatality risk of 53.5 percent, varying across clades and demographic groups.}",
    "}",
    "",
    "@article{cowling2013comparative,",
    "  title = {Comparative epidemiology of human infections with avian influenza {A H7N9} and {H5N1} viruses in {China}: a population-based study of laboratory-confirmed cases},",
    "  author = {Cowling, Benjamin J. and Jin, Lianmei and Lau, Eric H. Y. and Liao, Qiaohong and Wu, Peng and Jiang, Hui and Tsang, Tim K. and Zheng, Jiandong and Fang, Vicky J. and Chang, Zhaorui and Ni, Michael Y. and Zhang, Qian and Ip, Dennis K. M. and Yu, Jianxing and Li, Yu and Wang, Liping and Tu, Wenxiao and Meng, Ling and Wu, Joseph T. and Luo, Huiming and Li, Qun and Shu, Yuelong and Li, Zhongjie and Feng, Zijiang and Yang, Weizhong and Wang, Yu and Leung, Gabriel M. and Yu, Hongjie},",
    "  journal = {The Lancet},",
    "  year = {2013},",
    "  pmid = {23803488},",
    "  abstract = {We compared clinical severity and fatality risk of laboratory-confirmed human H7N9 and H5N1 infections in China, finding substantially higher mortality risk among H5N1 patients.}",
    "}"
  )
  writeLines(bib_lines, dest_path, useBytes = TRUE)
  dest_path
}

#' Optimize semantically equivalent PubMed and Google Scholar queries for a BibTeX file
#'
#' End-to-end reverse literature review query optimizer: parses a `.bib` file of seed
#' papers, enriches metadata via PubMed E-utilities, extracts high-lift domain terms
#' and concept blocks, evaluates candidate query combinations along the sensitivity
#' (seed recall) vs. specificity (total PubMed hits) Pareto frontier, mines non-seed
#' PubMed results for concept-grouped exclusion (`NOT` / `-`) terms and prefix wildcards,
#' and writes a report bundle.
#'
#' @param bib_path Path to input `.bib` file of seed papers.
#' @param output_dir Output directory for reports, CSVs, and Pareto plot (or `NULL` to skip disk export).
#' @param config_path Optional YAML config file path.
#' @param min_recall Minimum desired seed recall (0..1), default 0.60.
#' @param max_hits Optional maximum PubMed hit count threshold, default 2000.
#' @param max_evals Maximum live API query evaluations, default 60.
#' @param max_gs_chars Maximum Google Scholar query length, default 256.
#' @param use_cache Logical; cache API responses under `cache_dir`, default TRUE.
#' @param cache_dir Optional cache directory path (defaults to `file.path(output_dir, ".cache")`).
#' @param custom_queries Optional character vector or list of custom PubMed queries.
#' @param email Optional email for NCBI polite pool.
#' @param ncbi_key Optional NCBI API key.
#' @return Optimization result list (invisibly).
#' @export
optimize_search_queries <- function(bib_path,
                                    output_dir = "output",
                                    config_path = NULL,
                                    min_recall = 0.60,
                                    max_hits = 2000,
                                    max_evals = 60L,
                                    max_gs_chars = 256L,
                                    use_cache = TRUE,
                                    cache_dir = NULL,
                                    custom_queries = NULL,
                                    email = Sys.getenv("NCBI_EMAIL", ""),
                                    ncbi_key = Sys.getenv("NCBI_API_KEY", "")) {
  config <- if (!is.null(config_path) && nzchar(config_path)) {
    if (!file.exists(config_path)) stop("Config file not found: ", config_path, call. = FALSE)
    yaml::read_yaml(config_path)
  } else {
    list()
  }

  if (!is.null(custom_queries) && length(custom_queries) > 0L) {
    existing_cq <- .normalize_custom_queries_config(config$custom_queries)
    extra_cq <- .normalize_custom_queries_config(custom_queries)
    config$custom_queries <- c(existing_cq, extra_cq)
  }

  max_gs_chars <- as.integer(config$max_gs_chars %||% max_gs_chars)
  if (is.null(cache_dir) || !nzchar(cache_dir)) {
    cache_dir <- file.path(if (!is.null(output_dir) && nzchar(output_dir)) output_dir else tempdir(), ".cache")
  }

  cli::cli_h1("revlitquery: PubMed & Google Scholar Search Query Optimizer")
  cli::cli_inform(c("i" = "Parsing BibTeX: {.path {bib_path}}"))
  seeds_raw <- parse_bibtex(bib_path)

  cli::cli_inform(c("i" = "Resolving {nrow(seeds_raw)} seed papers against PubMed..."))
  seeds_enriched <- enrich_seeds(
    seeds_df = seeds_raw,
    email = if (nzchar(email)) email else NULL,
    ncbi_key = if (nzchar(ncbi_key)) ncbi_key else NULL,
    cache_dir = cache_dir,
    use_cache = use_cache
  )

  cli::cli_inform(c("i" = "Extracting candidate terms and probing PubMed background frequencies..."))
  terms_df <- extract_candidate_terms(
    seeds_df = seeds_enriched,
    config = config,
    probe_pubmed = TRUE,
    email = if (nzchar(email)) email else NULL,
    ncbi_key = if (nzchar(ncbi_key)) ncbi_key else NULL,
    cache_dir = cache_dir,
    use_cache = use_cache
  )

  blocks <- build_concept_blocks(terms_df = terms_df, seeds_df = seeds_enriched, config = config)
  cli::cli_inform(c("i" = "Discovered {length(blocks)} concept block(s):"))
  for (b in blocks) {
    cli::cli_inform(c("*" = "{b$name}: {paste(b$terms, collapse = ' | ')}"))
  }

  res <- optimize_queries(
    seeds_df = seeds_enriched,
    terms_df = terms_df,
    blocks = blocks,
    config = config,
    min_recall = min_recall,
    max_hits = max_hits,
    max_evals = max_evals,
    max_gs_chars = max_gs_chars,
    email = if (nzchar(email)) email else NULL,
    ncbi_key = if (nzchar(ncbi_key)) ncbi_key else NULL,
    cache_dir = cache_dir,
    use_cache = use_cache
  )

  if (!is.null(output_dir) && nzchar(output_dir)) {
    files <- write_optimization_report(res, output_dir = output_dir)
    res$output_files <- files
  }

  cli::cli_h2("Pareto Frontier Summary")
  p_df <- res$pareto_df
  for (i in seq_len(nrow(p_df))) {
    r <- p_df[i, ]
    star <- if (isTRUE(r$is_recommended)) "*" else " "
    cat(sprintf(
      "%s Q%-2d | Recall: %5.1f%% (%d/%d) | PM hits: %-8s | GS chars: %3d/%d\n",
      star, r$rank, 100 * r$combined_recall, r$seeds_found, r$seeds_resolvable,
      format(r$pm_hits, big.mark = ","), r$gs_chars, max_gs_chars
    ))
    cat(sprintf("     PM: %s\n", r$pubmed_query))
    cat(sprintf("     GS: %s\n", r$google_scholar_query))
  }
  cust_df <- res$eval_df[isTRUE(res$eval_df$is_custom) | (!is.na(res$eval_df$is_custom) & res$eval_df$is_custom), , drop = FALSE]
  if (nrow(cust_df) > 0L) {
    cli::cli_h2("Custom Benchmark Queries")
    for (ci in seq_len(nrow(cust_df))) {
      r <- cust_df[ci, ]
      lbl <- if (nzchar(r$custom_label)) r$custom_label else paste0("Custom #", r$eval_id)
      p_tag <- if (isTRUE(r$is_pareto)) " [Pareto]" else ""
      cat(sprintf(
        "  %-18s%s | Recall: %5.1f%% (%d/%d) | PM hits: %-8s | GS chars: %3d\n",
        lbl, p_tag, 100 * r$combined_recall, r$seeds_found, r$seeds_resolvable,
        format(r$pm_hits, big.mark = ","), r$gs_chars
      ))
      cat(sprintf("     PM: %s\n", r$pubmed_query))
      cat(sprintf("     GS: %s\n", r$google_scholar_query))
    }
  }
  eh_df <- res$exclusion_hints
  if (!is.null(eh_df) && nrow(eh_df) > 0L) {
    cli::cli_h2("Top Exclusion Hints (Grouped by Concept)")
    for (j in seq_len(min(15L, nrow(eh_df)))) {
      h <- eh_df[j, ]
      cat(sprintf(
        "  [%-20s] %-22s | -%s hits (-%.1f%%) | %s | GS: %s\n",
        h$concept %||% "general", h$exclude_spec, format(h$exact_hits_removed, big.mark = ","),
        h$pct_hits_reduced, h$safety_tier, h$gs_clause
      ))
    }
  }
  if (!is.null(output_dir) && nzchar(output_dir)) {
    cli::cli_inform(c("v" = "Report and artifacts saved to {.path {output_dir}}"))
  }
  invisible(res)
}

.revlitquery_ui <- function() {
  tags <- shiny::tags
  HTML <- shiny::HTML
  bslib::page_sidebar(
    title = tags$span(
      fontawesome::fa_i("magnifying-glass-chart", class = "me-2"),
      "revlitquery: PubMed & Google Scholar Query Optimizer"
    ),
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#1f6feb"),
    fillable = FALSE,
    tags$head(
      tags$style(HTML("
        .query-box {
          background: #f8fafc;
          border: 1px solid #e2e8f0;
          border-radius: 6px;
          padding: 10px 12px;
          font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
          font-size: 0.86rem;
          white-space: pre-wrap;
          word-break: break-word;
          margin-bottom: 0;
        }
        .seed-status-pill {
          font-size: 0.82rem;
          padding: 6px 10px;
          border-radius: 6px;
          background: #f1f5f9;
          border: 1px solid #e2e8f0;
        }
        .block-pill {
          display: inline-block;
          background: #eef2ff;
          color: #1e40af;
          border: 1px solid #c7d2fe;
          border-radius: 999px;
          padding: 2px 9px;
          font-size: 0.8rem;
          margin: 2px 4px 2px 0;
        }
        .excl-pill {
          display: inline-block;
          background: #fef2f2;
          color: #991b1b;
          border: 1px solid #fecaca;
          border-radius: 999px;
          padding: 2px 9px;
          font-size: 0.8rem;
          margin: 2px 4px 2px 0;
        }
      ")),
      tags$script(HTML("
        window.copyQueryText = function(elId, btnEl) {
          var el = document.getElementById(elId);
          if (!el) return;
          var text = el.innerText || el.textContent || '';
          var finish = function() {
            if (!btnEl) return;
            var orig = btnEl.dataset.origLabel || btnEl.innerText;
            btnEl.dataset.origLabel = orig;
            btnEl.innerText = 'Copied!';
            setTimeout(function() { btnEl.innerText = orig; }, 1400);
          };
          if (navigator.clipboard && navigator.clipboard.writeText) {
            navigator.clipboard.writeText(text).then(finish).catch(function() {
              var ta = document.createElement('textarea');
              ta.value = text;
              document.body.appendChild(ta);
              ta.select();
              document.execCommand('copy');
              document.body.removeChild(ta);
              finish();
            });
          } else {
            var ta = document.createElement('textarea');
            ta.value = text;
            document.body.appendChild(ta);
            ta.select();
            document.execCommand('copy');
            document.body.removeChild(ta);
            finish();
          }
        };
      "))
    ),
    sidebar = bslib::sidebar(
      width = 345,
      tags$h6("1. Seed Papers (BibTeX)", class = "fw-bold mb-2"),
      shiny::fileInput("bib_file", "Upload .bib file", accept = c(".bib", ".txt"), width = "100%"),
      shiny::actionButton(
        "load_sample", "Load sample H5N1 .bib",
        icon = fontawesome::fa_i("flask"),
        class = "btn-outline-secondary btn-sm w-100 mb-2"
      ),
      shiny::uiOutput("bib_status_ui"),
      tags$hr(class = "my-3"),
      tags$h6("2. Sensitivity vs. Specificity", class = "fw-bold mb-2"),
      shiny::sliderInput(
        "min_recall", "Minimum seed recall (%)",
        min = 20, max = 100, value = 75, step = 5, post = "%"
      ),
      shiny::numericInput(
        "max_hits", "Max PubMed hits cap (<= 2000 recommended)",
        value = 2000, min = 1, step = 100
      ),
      shiny::sliderInput(
        "max_gs_chars", "Google Scholar char limit",
        min = 100, max = 256, value = 256, step = 8
      ),
      bslib::accordion(
        id = "adv_acc",
        open = FALSE,
        class = "mt-2",
        bslib::accordion_panel(
          "Search & API Settings",
          shiny::sliderInput("max_evals", "Max live query evaluations", min = 10, max = 120, value = 60, step = 5),
          shiny::sliderInput("max_blocks", "Max concept blocks", min = 1, max = 4, value = 3, step = 1),
          shiny::sliderInput("max_terms_per_block", "Max terms per block", min = 1, max = 8, value = 5, step = 1),
          bslib::input_switch("use_cache", "Use disk cache (.cache/)", value = TRUE),
          shiny::textInput("email", "API email (NCBI polite pool)", value = Sys.getenv("NCBI_EMAIL", "")),
          shiny::passwordInput("ncbi_key", "NCBI API key (optional)", value = Sys.getenv("NCBI_API_KEY", ""))
        ),
        bslib::accordion_panel(
          "Term Overrides & Custom Blocks",
          shiny::textAreaInput(
            "required_terms", "Required terms (comma/newline)",
            placeholder = "e.g. h5n1", rows = 2
          ),
          shiny::textAreaInput(
            "banned_terms", "Banned candidate terms",
            placeholder = "e.g. world health organization", rows = 2
          ),
          shiny::textAreaInput(
            "exclude_terms", "Exclude terms (NOT / -)",
            placeholder = "e.g. vaccin*, poultry, h7n9", rows = 2
          ),
          shiny::textAreaInput(
            "custom_blocks_yaml", "Custom concept blocks (YAML)",
            placeholder = "pathogen:\n  - h5n1\n  - h5\nseverity:\n  - fatal*\n  - mortality",
            rows = 5
          ),
          shiny::actionButton(
            "copy_blocks_to_editor", "Copy discovered blocks into editor",
            icon = fontawesome::fa_i("wand-magic-sparkles"),
            class = "btn-outline-secondary btn-sm w-100 mt-1"
          )
        )
      ),
      shiny::actionButton(
        "run_opt", "Run Query Optimizer",
        icon = fontawesome::fa_i("magnifying-glass"),
        class = "btn-primary w-100 fw-semibold mt-3"
      )
    ),
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      fill = FALSE,
      bslib::value_box(
        title = "Seed Papers",
        value = shiny::textOutput("vb_seeds", inline = TRUE),
        showcase = fontawesome::fa_i("book-open"),
        theme = "light",
        tags$p(shiny::textOutput("vb_seeds_sub", inline = TRUE), class = "mb-0 small text-muted")
      ),
      bslib::value_box(
        title = "Selected Recall",
        value = shiny::textOutput("vb_recall", inline = TRUE),
        showcase = fontawesome::fa_i("bullseye"),
        theme = "primary",
        tags$p(shiny::textOutput("vb_recall_sub", inline = TRUE), class = "mb-0 small")
      ),
      bslib::value_box(
        title = "PubMed Total Hits",
        value = shiny::textOutput("vb_pm_hits", inline = TRUE),
        showcase = fontawesome::fa_i("database"),
        theme = "info",
        tags$p(shiny::textOutput("vb_pm_hits_sub", inline = TRUE), class = "mb-0 small")
      ),
      bslib::value_box(
        title = "Google Scholar Chars",
        value = shiny::textOutput("vb_gs_chars", inline = TRUE),
        showcase = fontawesome::fa_i("graduation-cap"),
        theme = "secondary",
        tags$p(shiny::textOutput("vb_gs_sub", inline = TRUE), class = "mb-0 small")
      )
    ),
    bslib::navset_card_underline(
      id = "main_tabs",
      bslib::nav_panel(
        "Pareto Explorer",
        icon = fontawesome::fa_i("chart-line"),
        bslib::card_body(
          plotly::plotlyOutput("pareto_plot", height = "340px"),
          bslib::accordion(
            id = "custom_q_acc",
            open = FALSE,
            class = "mt-3",
            bslib::accordion_panel(
              "Add & Analyze a Custom Query",
              icon = fontawesome::fa_i("plus-circle"),
              shiny::textAreaInput(
                "custom_pm_query", "PubMed query to evaluate",
                rows = 2,
                placeholder = 'e.g. (("H5N1"[Title/Abstract] OR "H5"[Title/Abstract]) AND "Influenza"[Title/Abstract] AND "fatal*"[Title/Abstract]) NOT ("H7N9"[Title/Abstract] OR "H1N1"[Title/Abstract])'
              ),
              bslib::layout_columns(
                col_widths = c(4, 8),
                shiny::textInput("custom_query_label", "Label (optional)", placeholder = "e.g. Manuscript original"),
                shiny::textInput("custom_gs_query", "Google Scholar query (optional; auto-translated if blank)")
              ),
              shiny::actionButton(
                "add_custom_query", "Evaluate & Add Query",
                icon = fontawesome::fa_i("plus"),
                class = "btn-outline-primary btn-sm fw-semibold"
              )
            )
          ),
          tags$div(class = "mt-3", shiny::uiOutput("selected_query_card")),
          tags$h6("Pareto-Optimal Query Pairs", class = "fw-bold mt-4 mb-2"),
          DT::DTOutput("pareto_table")
        )
      ),
      bslib::nav_panel(
        "Seed Coverage Matrix",
        icon = fontawesome::fa_i("table-cells"),
        bslib::card_body(
          tags$p(
            "Inspect which BibTeX seed papers are resolved in PubMed and retrieved by each Pareto query.",
            class = "text-muted small mb-2"
          ),
          DT::DTOutput("coverage_table")
        )
      ),
      bslib::nav_panel(
        "Terms & Concept Blocks",
        icon = fontawesome::fa_i("tags"),
        bslib::card_body(
          shiny::uiOutput("blocks_summary_ui"),
          tags$h6("Candidate Terms Ranked by Seed Lift", class = "fw-bold mt-3 mb-2"),
          DT::DTOutput("terms_table")
        )
      ),
      bslib::nav_panel(
        "Exclusion Hints",
        icon = fontawesome::fa_i("filter-circle-xmark"),
        bslib::card_body(
          tags$p(
            "Terms and prefix wildcards frequently appearing in non-seed PubMed results, organized by exclusion concept. Select rows or add the greedy set-cover exclusions to the sidebar.",
            class = "text-muted small mb-2"
          ),
          shiny::uiOutput("exclusion_blocks_ui"),
          tags$div(
            class = "d-flex flex-wrap gap-2 mb-3 mt-2",
            shiny::actionButton("add_selected_excl", "Add selected rows to sidebar", class = "btn-primary btn-sm"),
            shiny::actionButton("add_setcover_excl", "Add top set-cover safe exclusions", class = "btn-outline-primary btn-sm"),
            shiny::actionButton("add_safe_excl", "Add all safe (0-loss) exclusions", class = "btn-outline-secondary btn-sm")
          ),
          DT::DTOutput("exclusion_hints_table")
        )
      ),
      bslib::nav_panel(
        "All Evaluated Queries",
        icon = fontawesome::fa_i("list-check"),
        bslib::card_body(
          DT::DTOutput("all_eval_table")
        )
      ),
      bslib::nav_panel(
        "Export & Report",
        icon = fontawesome::fa_i("file-export"),
        bslib::card_body(
          tags$div(
            class = "d-flex flex-wrap gap-2 mb-3",
            shiny::downloadButton("dl_zip", "Download Full Report Bundle (.zip)", class = "btn-primary btn-sm"),
            shiny::downloadButton("dl_pareto_csv", "Download Pareto Queries (.csv)", class = "btn-outline-secondary btn-sm"),
            shiny::downloadButton("dl_config_yml", "Download CLI Config (.yml)", class = "btn-outline-secondary btn-sm")
          ),
          shiny::uiOutput("report_preview_ui")
        )
      )
    )
  )
}

.revlitquery_make_server <- function(initial_result = NULL, initial_bib_label = NULL) {
  function(input, output, session) {
  tags <- shiny::tags
  HTML <- shiny::HTML
  req <- shiny::req
  init_seeds <- if (!is.null(initial_result)) initial_result$seeds_df else NULL
  init_lbl <- if (!is.null(initial_bib_label) && nzchar(initial_bib_label)) {
    initial_bib_label
  } else if (!is.null(initial_result)) {
    "Preloaded optimization result"
  } else {
    "No BibTeX loaded yet"
  }
  init_sel <- if (!is.null(initial_result) && !is.null(initial_result$pareto_df) && nrow(initial_result$pareto_df) > 0L) {
    rec_eids <- initial_result$pareto_df$eval_id[initial_result$pareto_df$is_recommended]
    if (length(rec_eids) > 0L) rec_eids[1L] else initial_result$pareto_df$eval_id[1L]
  } else {
    NULL
  }
  seeds_raw <- shiny::reactiveVal(init_seeds)
  bib_label <- shiny::reactiveVal(init_lbl)
  opt_result <- shiny::reactiveVal(initial_result)
  selected_eval_id <- shiny::reactiveVal(init_sel)

  load_bib_path <- function(path, label) {
    tryCatch({
      df <- parse_bibtex(path)
      seeds_raw(df)
      bib_label(label)
      opt_result(NULL)
      selected_eval_id(NULL)
      shiny::showNotification(
        sprintf("Loaded %d BibTeX entries from %s.", nrow(df), label),
        type = "message", duration = 3
      )
    }, error = function(e) {
      shiny::showNotification(paste("Failed to parse BibTeX:", conditionMessage(e)), type = "error", duration = 6)
    })
  }

  shiny::observeEvent(input$bib_file, {
    req(input$bib_file$datapath)
    load_bib_path(input$bib_file$datapath, input$bib_file$name)
  })

  shiny::observeEvent(input$load_sample, {
    sample_path <- .write_demo_bibtex()
    load_bib_path(sample_path, "sample_h5n1_cfr.bib")
  })

  output$bib_status_ui <- shiny::renderUI({
    df <- seeds_raw()
    if (is.null(df)) {
      tags$div(class = "seed-status-pill text-muted", "Upload a .bib file or load the sample dataset.")
    } else {
      n_abs <- sum(nzchar(df$abstract))
      tags$div(
        class = "seed-status-pill",
        tags$strong(bib_label()),
        tags$div(sprintf("%d entries (%d with BibTeX abstract)", nrow(df), n_abs), class = "text-muted")
      )
    }
  })

  shiny::observeEvent(input$run_opt, {
    df_raw <- seeds_raw()
    if (is.null(df_raw) || nrow(df_raw) == 0L) {
      shiny::showNotification("Please load a BibTeX file first.", type = "warning")
      return()
    }

    cfg_res <- build_config_from_inputs(
      required_txt = input$required_terms,
      banned_txt = input$banned_terms,
      exclude_txt = input$exclude_terms,
      custom_blocks_yaml = input$custom_blocks_yaml,
      max_blocks = input$max_blocks,
      max_terms_per_block = input$max_terms_per_block,
      max_gs_chars = input$max_gs_chars
    )
    if (!isTRUE(cfg_res$ok)) {
      shiny::showNotification(cfg_res$error, type = "error", duration = 7)
      return()
    }
    cfg <- cfg_res$config

    min_rec <- as.numeric(input$min_recall) / 100
    max_h <- if (is.null(input$max_hits) || is.na(input$max_hits) || input$max_hits <= 0) Inf else as.numeric(input$max_hits)
    email_v <- if (nzchar(stringr::str_trim(input$email))) stringr::str_trim(input$email) else NULL
    ncbi_v <- if (nzchar(stringr::str_trim(input$ncbi_key))) stringr::str_trim(input$ncbi_key) else NULL
    cache_dir <- ".cache"

    res <- tryCatch({
      shiny::withProgress(message = "Optimizing search queries...", value = 0, {
        shiny::incProgress(0.05, detail = "Enriching BibTeX seeds via PubMed")
        enriched <- enrich_seeds(
          seeds_df = df_raw,
          email = email_v,
          ncbi_key = ncbi_v,
          cache_dir = cache_dir,
          use_cache = isTRUE(input$use_cache),
          progress_cb = function(msg) shiny::incProgress(0.15 / max(nrow(df_raw), 1L), detail = msg)
        )

        shiny::incProgress(0.05, detail = "Extracting candidate terms & probing frequencies")
        terms_df <- extract_candidate_terms(
          seeds_df = enriched,
          config = cfg,
          probe_pubmed = TRUE,
          email = email_v,
          ncbi_key = ncbi_v,
          cache_dir = cache_dir,
          use_cache = isTRUE(input$use_cache),
          progress_cb = function(msg) shiny::incProgress(0.25 / 55, detail = msg)
        )

        blocks <- build_concept_blocks(terms_df = terms_df, seeds_df = enriched, config = cfg)
        if (length(blocks) == 0L) {
          stop("Could not extract any concept blocks from these seed papers.", call. = FALSE)
        }

        shiny::incProgress(0.05, detail = "Evaluating candidate query frontier")
        max_ev <- as.integer(input$max_evals)
        optimize_queries(
          seeds_df = enriched,
          terms_df = terms_df,
          blocks = blocks,
          config = cfg,
          min_recall = min_rec,
          max_hits = max_h,
          max_evals = max_ev,
          max_gs_chars = as.integer(input$max_gs_chars),
          email = email_v,
          ncbi_key = ncbi_v,
          cache_dir = cache_dir,
          use_cache = isTRUE(input$use_cache),
          progress_cb = function(msg) shiny::incProgress(0.45 / max(max_ev, 1L), detail = msg)
        )
      })
    }, error = function(e) {
      shiny::showNotification(paste("Optimization failed:", conditionMessage(e)), type = "error", duration = 8)
      NULL
    })

    if (!is.null(res)) {
      prev_res <- opt_result()
      if (!is.null(prev_res) && !is.null(prev_res$eval_df) && nrow(prev_res$eval_df) > 0L &&
          ("is_custom" %in% names(prev_res$eval_df)) && any(prev_res$eval_df$is_custom) &&
          nrow(prev_res$seeds_df) == nrow(res$seeds_df)) {
        prev_cust <- prev_res$eval_df[prev_res$eval_df$is_custom, , drop = FALSE]
        for (ci in seq_len(nrow(prev_cust))) {
          res <- add_custom_query_to_result(res, prev_cust[ci, , drop = FALSE])
        }
      }
      opt_result(res)
      rec_eids <- res$pareto_df$eval_id[res$pareto_df$is_recommended]
      selected_eval_id(if (length(rec_eids) > 0L) rec_eids[1L] else if (nrow(res$pareto_df) > 0L) res$pareto_df$eval_id[1L] else NULL)
      shiny::showNotification(
        sprintf("Done! Evaluated %d queries (%d Pareto-optimal).", nrow(res$eval_df), nrow(res$pareto_df)),
        type = "message", duration = 4
      )
    }
  })

  shiny::observeEvent(input$add_custom_query, {
    df_raw <- seeds_raw()
    if (is.null(df_raw) || nrow(df_raw) == 0L) {
      shiny::showNotification("Please load a BibTeX file first.", type = "warning")
      return()
    }
    pm_txt <- stringr::str_squish(as.character(input$custom_pm_query %||% ""))
    if (!nzchar(pm_txt)) {
      shiny::showNotification("Enter a PubMed query to evaluate.", type = "warning")
      return()
    }

    min_rec <- as.numeric(input$min_recall) / 100
    max_h <- if (is.null(input$max_hits) || is.na(input$max_hits) || input$max_hits <= 0) Inf else as.numeric(input$max_hits)
    email_v <- if (nzchar(stringr::str_trim(input$email))) stringr::str_trim(input$email) else NULL
    ncbi_v <- if (nzchar(stringr::str_trim(input$ncbi_key))) stringr::str_trim(input$ncbi_key) else NULL
    cache_dir <- ".cache"

    updated <- tryCatch({
      shiny::withProgress(message = "Evaluating custom query...", value = 0.1, {
        cur <- opt_result()
        if (is.null(cur)) {
          shiny::incProgress(0.25, detail = "Resolving BibTeX seed papers in PubMed...")
          enriched <- enrich_seeds(
            seeds_df = df_raw,
            email = email_v, ncbi_key = ncbi_v,
            cache_dir = cache_dir, use_cache = isTRUE(input$use_cache)
          )
          cur <- list(
            eval_df = tibble::tibble(),
            pareto_df = tibble::tibble(),
            exclusion_hints = .empty_exclusion_hints(),
            recommended_sig = NA_character_,
            recommended_status = "empty",
            min_recall = min_rec,
            max_hits = max_h,
            max_gs_chars = as.integer(input$max_gs_chars),
            seeds_df = enriched,
            terms_df = tibble::tibble(),
            blocks = list()
          )
        }
        cur$min_recall <- min_rec
        cur$max_hits <- max_h

        shiny::incProgress(0.55, detail = "Querying PubMed hit count and seed coverage...")
        c_row <- evaluate_custom_query(
          pm_query = pm_txt,
          seeds_df = cur$seeds_df,
          label = input$custom_query_label,
          gs_query = input$custom_gs_query,
          email = email_v, ncbi_key = ncbi_v,
          cache_dir = cache_dir, use_cache = isTRUE(input$use_cache)
        )
        add_custom_query_to_result(cur, c_row)
      })
    }, error = function(e) {
      shiny::showNotification(paste("Custom query failed:", conditionMessage(e)), type = "error", duration = 7)
      NULL
    })

    if (!is.null(updated)) {
      opt_result(updated)
      selected_eval_id(updated$last_added_eval_id)
      added_row <- dplyr::filter(updated$eval_df, .data$eval_id == updated$last_added_eval_id)
      shiny::showNotification(
        sprintf(
          "Added custom query: %.1f%% recall (%d/%d), %s PubMed hits.",
          100 * added_row$combined_recall[1L],
          added_row$seeds_found[1L], added_row$seeds_resolvable[1L],
          format(added_row$pm_hits[1L], big.mark = ",")
        ),
        type = "message", duration = 5
      )
    }
  })

  current_result <- shiny::reactive({
    res <- opt_result()
    if (is.null(res)) return(NULL)
    min_rec <- as.numeric(input$min_recall) / 100
    max_h <- if (is.null(input$max_hits) || is.na(input$max_hits) || input$max_hits <= 0) Inf else as.numeric(input$max_hits)
    eval_df <- res$eval_df
    if (nrow(eval_df) == 0L) return(res)
    eval_df$within_max_hits <- is.na(eval_df$pm_hits) | (eval_df$pm_hits <= max_h)
    eval_df$is_pareto <- compute_pareto_frontier(eval_df, max_hits = max_h)
    pareto_df <- eval_df |>
      dplyr::filter(.data$is_pareto) |>
      dplyr::arrange(dplyr::desc(.data$combined_recall), .data$pm_hits, .data$gs_chars) |>
      dplyr::mutate(rank = dplyr::row_number())
    rec <- recompute_recommendation(pareto_df, min_recall = min_rec, max_hits = max_h)
    res$eval_df <- eval_df
    res$pareto_df <- rec$pareto_df
    res$recommended_sig <- rec$recommended_sig
    res$recommended_status <- rec$recommended_status
    res$min_recall <- min_rec
    res$max_hits <- max_h
    res
  })

  shiny::observeEvent(list(input$min_recall, input$max_hits), {
    res <- current_result()
    if (is.null(res) || nrow(res$pareto_df) == 0L) return()
    rec_e <- res$pareto_df$eval_id[res$pareto_df$is_recommended]
    if (length(rec_e) > 0L) selected_eval_id(rec_e[1L])
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$pareto_table_rows_selected, {
    res <- current_result()
    idx <- input$pareto_table_rows_selected
    if (!is.null(res) && length(idx) == 1L && idx >= 1L && idx <= nrow(res$pareto_df)) {
      selected_eval_id(res$pareto_df$eval_id[idx])
    }
  })

  shiny::observeEvent(input$all_eval_table_rows_selected, {
    res <- current_result()
    idx <- input$all_eval_table_rows_selected
    if (!is.null(res) && nrow(res$eval_df) > 0L && length(idx) == 1L) {
      ord_df <- res$eval_df |>
        dplyr::arrange(dplyr::desc(.data$is_pareto), dplyr::desc(.data$combined_recall), .data$pm_hits)
      if (idx >= 1L && idx <= nrow(ord_df)) {
        selected_eval_id(ord_df$eval_id[idx])
      }
    }
  })

  shiny::observeEvent(plotly::event_data("plotly_click", source = "pareto"), {
    ed <- plotly::event_data("plotly_click", source = "pareto")
    if (!is.null(ed) && !is.null(ed$customdata)) {
      eid <- suppressWarnings(as.integer(ed$customdata[1L]))
      if (!is.na(eid)) selected_eval_id(eid)
    }
  })

  shiny::observeEvent(input$copy_blocks_to_editor, {
    res <- current_result()
    if (is.null(res) || length(res$blocks) == 0L) {
      shiny::showNotification("Run the optimizer first to discover concept blocks.", type = "warning")
      return()
    }
    yml_txt <- blocks_to_yaml(res$blocks)
    shiny::updateTextAreaInput(session, "custom_blocks_yaml", value = yml_txt)
    shiny::showNotification("Copied discovered blocks into Custom Concept Blocks YAML editor.", type = "message", duration = 3)
  })

  shiny::observeEvent(input$copy_blocks_inline, {
    res <- current_result()
    if (!is.null(res) && length(res$blocks) > 0L) {
      shiny::updateTextAreaInput(session, "custom_blocks_yaml", value = blocks_to_yaml(res$blocks))
      shiny::showNotification("Copied discovered blocks into sidebar YAML editor.", type = "message", duration = 3)
    }
  })

  selected_query_row <- shiny::reactive({
    res <- current_result()
    if (is.null(res) || nrow(res$eval_df) == 0L) return(NULL)
    eid <- selected_eval_id()
    row <- if (!is.null(eid)) dplyr::filter(res$eval_df, .data$eval_id == eid) else tibble::tibble()
    if (nrow(row) == 0L && nrow(res$pareto_df) > 0L) {
      row <- dplyr::filter(res$pareto_df, .data$is_recommended)
      if (nrow(row) == 0L) row <- res$pareto_df[1L, , drop = FALSE]
    }
    if (nrow(row) == 0L) row <- res$eval_df[1L, , drop = FALSE]
    r1 <- row[1L, , drop = FALSE]
    p_match <- dplyr::filter(res$pareto_df, .data$eval_id == r1$eval_id)
    r1$pareto_rank <- if (nrow(p_match) > 0L) p_match$rank[1L] else NA_integer_
    r1$is_recommended <- if (nrow(p_match) > 0L) isTRUE(p_match$is_recommended[1L]) else FALSE
    r1
  })

  output$vb_seeds <- shiny::renderText({
    res <- current_result()
    if (!is.null(res)) return(as.character(nrow(res$seeds_df)))
    df <- seeds_raw()
    if (!is.null(df)) as.character(nrow(df)) else "--"
  })

  output$vb_seeds_sub <- shiny::renderText({
    res <- current_result()
    if (!is.null(res)) {
      return(sprintf("%d in PubMed | %d with abstract", sum(res$seeds_df$in_pubmed), sum(nzchar(res$seeds_df$full_abstract))))
    }
    df <- seeds_raw()
    if (!is.null(df)) "Parsed (click Run to resolve)" else "Load .bib to start"
  })

  output$vb_recall <- shiny::renderText({
    r <- selected_query_row()
    if (is.null(r)) return("--")
    sprintf("%.1f%%", 100 * r$combined_recall)
  })

  output$vb_recall_sub <- shiny::renderText({
    r <- selected_query_row()
    if (is.null(r)) return("Waiting for optimization")
    lbl <- if (!is.na(r$pareto_rank)) {
      paste0("Q", r$pareto_rank)
    } else if (isTRUE(r$is_custom) && nzchar(r$custom_label)) {
      r$custom_label
    } else {
      paste0("#", r$eval_id)
    }
    sprintf("%d / %d PubMed seeds (%s)", r$seeds_found, r$seeds_resolvable, lbl)
  })

  output$vb_pm_hits <- shiny::renderText({
    r <- selected_query_row()
    if (is.null(r) || is.na(r$pm_hits)) return("--")
    format(r$pm_hits, big.mark = ",")
  })

  output$vb_pm_hits_sub <- shiny::renderText({
    r <- selected_query_row()
    res <- current_result()
    if (is.null(r) || is.null(res) || is.na(r$pm_hits)) return("Target: <= max_hits cap")
    if (is.finite(res$max_hits)) {
      if (r$pm_hits <= res$max_hits) {
        sprintf("Within <= %s cap", format(res$max_hits, big.mark = ","))
      } else {
        sprintf("Exceeds %s cap", format(res$max_hits, big.mark = ","))
      }
    } else {
      "No hit cap set"
    }
  })

  output$vb_gs_chars <- shiny::renderText({
    r <- selected_query_row()
    if (is.null(r)) return("--")
    sprintf("%d / %d", r$gs_chars, as.integer(input$max_gs_chars))
  })

  output$vb_gs_sub <- shiny::renderText({
    r <- selected_query_row()
    if (is.null(r)) return("<= 256 char limit")
    if (r$gs_chars <= 256L) "Within Google Scholar limit" else "Exceeds 256 chars"
  })

  output$pareto_plot <- plotly::renderPlotly({
    res <- current_result()
    min_rec <- as.numeric(input$min_recall) / 100
    max_h <- if (is.null(input$max_hits) || is.na(input$max_hits) || input$max_hits <= 0) Inf else as.numeric(input$max_hits)
    if (is.null(res)) {
      return(build_pareto_plotly(tibble::tibble(), tibble::tibble(), NULL, min_rec, max_h))
    }
    build_pareto_plotly(
      eval_df = res$eval_df,
      pareto_df = res$pareto_df,
      selected_rank = NULL,
      min_recall = min_rec,
      max_hits = max_h,
      selected_eval_id = selected_eval_id()
    )
  })

  output$selected_query_card <- shiny::renderUI({
    r <- selected_query_row()
    res <- current_result()
    if (is.null(r) || is.null(res)) {
      return(tags$div(
        class = "alert alert-light border small text-muted mb-0",
        "Load a BibTeX file and click ", tags$strong("Run Query Optimizer"),
        " (or evaluate a custom query above) to inspect PubMed and Google Scholar queries."
      ))
    }

    title_txt <- if (!is.na(r$pareto_rank)) {
      sprintf("Inspector: Q%d", r$pareto_rank)
    } else if (isTRUE(r$is_custom)) {
      lbl <- if (nzchar(r$custom_label)) r$custom_label else paste0("Custom #", r$eval_id)
      paste0("Inspector: ", lbl)
    } else {
      sprintf("Inspector: Evaluated #%d", r$eval_id)
    }
    badge <- if (isTRUE(r$is_recommended)) {
      tags$span(class = "badge bg-danger ms-2", "Recommended")
    } else if (isTRUE(r$is_custom) && isTRUE(r$is_pareto)) {
      tags$span(class = "badge bg-primary ms-2", "Custom (Pareto-optimal)")
    } else if (isTRUE(r$is_custom)) {
      tags$span(class = "badge text-bg-secondary ms-2", "Custom Query")
    } else if (isTRUE(r$is_pareto)) {
      tags$span(class = "badge bg-primary ms-2", "Pareto-optimal")
    } else {
      tags$span(class = "badge text-bg-light border ms-2", "Evaluated")
    }
    missed_txt <- if (nzchar(r$missed_citekeys)) r$missed_citekeys else "None"

    bslib::card(
      class = "border-primary-subtle",
      bslib::card_header(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
        tags$span(tags$strong(title_txt), badge),
        tags$span(
          class = "small text-muted",
          sprintf(
            "Recall: %.1f%% (%d/%d) | PubMed: %s hits | GS: %d/%d chars",
            100 * r$combined_recall, r$seeds_found, r$seeds_resolvable,
            format(r$pm_hits, big.mark = ","), r$gs_chars, res$max_gs_chars
          )
        )
      ),
      bslib::card_body(
        bslib::layout_columns(
          col_widths = c(6, 6),
          tags$div(
            tags$div(
              class = "d-flex justify-content-between align-items-center mb-1",
              tags$strong("PubMed Query"),
              tags$div(
                tags$button(
                  type = "button", class = "btn btn-outline-secondary btn-xs btn-sm py-0 px-2 me-1",
                  onclick = "window.copyQueryText('pm_q_pre', this)", "Copy"
                ),
                tags$a(
                  href = r$pubmed_url, target = "_blank", rel = "noopener",
                  class = "btn btn-outline-primary btn-sm py-0 px-2",
                  "Open in PubMed"
                )
              )
            ),
            tags$pre(id = "pm_q_pre", class = "query-box", r$pubmed_query)
          ),
          tags$div(
            tags$div(
              class = "d-flex justify-content-between align-items-center mb-1",
              tags$strong("Google Scholar Query"),
              tags$div(
                tags$button(
                  type = "button", class = "btn btn-outline-secondary btn-xs btn-sm py-0 px-2 me-1",
                  onclick = "window.copyQueryText('gs_q_pre', this)", "Copy"
                ),
                tags$a(
                  href = r$google_scholar_url, target = "_blank", rel = "noopener",
                  class = "btn btn-outline-primary btn-sm py-0 px-2",
                  "Open in Google Scholar"
                )
              )
            ),
            tags$pre(id = "gs_q_pre", class = "query-box", r$google_scholar_query)
          )
        ),
        tags$div(
          class = "mt-2 small text-muted",
          tags$span(tags$strong("Missed seeds: "), missed_txt)
        )
      )
    )
  })

  output$pareto_table <- DT::renderDT({
    res <- current_result()
    req(res, nrow(res$pareto_df) > 0L)
    p_df <- res$pareto_df
    disp <- tibble::tibble(
      Query = paste0(
        "Q", p_df$rank,
        ifelse(p_df$is_recommended, " *", ""),
        ifelse(isTRUE(p_df$is_custom) | (!is.na(p_df$is_custom) & p_df$is_custom), " [Custom]", "")
      ),
      `Recall %` = round(100 * p_df$combined_recall, 1),
      Found = paste0(p_df$seeds_found, "/", p_df$seeds_resolvable),
      `PM Hits` = p_df$pm_hits,
      `GS Chars` = p_df$gs_chars,
      `PubMed Query` = p_df$pubmed_query,
      `Google Scholar Query` = p_df$google_scholar_query,
      `Missed Seeds` = p_df$missed_citekeys
    )
    sel_idx <- which(p_df$eval_id == (selected_eval_id() %||% -1L))
    DT::datatable(
      disp,
      selection = list(mode = "single", selected = if (length(sel_idx) == 1L) sel_idx else NULL),
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
    ) |>
      DT::formatCurrency("PM Hits", currency = "", digits = 0)
  })

  output$coverage_table <- DT::renderDT({
    res <- current_result()
    if (!is.null(res) && nrow(res$eval_df) > 0L) {
      cov_df <- build_coverage_display_df(res$seeds_df, res$pareto_df, res$eval_df)
      q_cols <- grep("^(Q[0-9]+|Custom: )", names(cov_df), value = TRUE)
      dt <- DT::datatable(cov_df, rownames = FALSE, selection = "none", options = list(pageLength = 20, scrollX = TRUE))
      if (length(q_cols) > 0L) {
        dt <- dt |>
          DT::formatStyle(
            q_cols,
            backgroundColor = DT::styleEqual(c("Found", "Missed", "N/A"), c("#dcfce7", "#fee2e2", "#f1f5f9")),
            color = DT::styleEqual(c("Found", "Missed", "N/A"), c("#166534", "#991b1b", "#64748b")),
            fontWeight = "bold"
          )
      }
      return(dt)
    }

    df <- seeds_raw()
    req(df)
    prev <- df |>
      dplyr::transmute(
        Citekey = .data$citekey,
        Year = .data$year,
        Title = .data$title,
        DOI = .data$doi,
        PMID = .data$pmid,
        `Has Abstract` = ifelse(nzchar(.data$abstract), "Yes", "No")
      )
    DT::datatable(prev, rownames = FALSE, selection = "none", options = list(pageLength = 20, scrollX = TRUE))
  })

  output$blocks_summary_ui <- shiny::renderUI({
    res <- current_result()
    if (is.null(res)) {
      return(tags$div(class = "text-muted small", "Run the optimizer to view discovered concept blocks and candidate term lift scores."))
    }
    blk_cards <- lapply(res$blocks, function(b) {
      tags$div(
        class = "mb-2",
        tags$strong(paste0(b$name, ": ")),
        lapply(b$terms, function(t_str) tags$span(class = "block-pill", t_str))
      )
    })
    dom_q <- attr(res$blocks, "domain_qualifier") %||% character(0)
    if (length(dom_q) > 0L) {
      blk_cards[[length(blk_cards) + 1L]] <- tags$div(
        class = "mb-2",
        tags$strong("domain_context: "),
        lapply(dom_q, function(t_str) tags$span(class = "block-pill", t_str))
      )
    }
    tags$div(
      class = "p-3 bg-light border rounded",
      tags$div(
        class = "d-flex justify-content-between align-items-center mb-2",
        tags$h6("Active Inclusion Concept Blocks", class = "fw-bold mb-0"),
        shiny::actionButton(
          "copy_blocks_inline", "Copy blocks into YAML editor",
          class = "btn-outline-secondary btn-sm"
        )
      ),
      blk_cards
    )
  })

  output$terms_table <- DT::renderDT({
    res <- current_result()
    req(res)
    disp_t <- res$terms_df |>
      dplyr::transmute(
        Term = .data$term,
        Words = .data$n_words,
        `Seed DF (TiAb)` = .data$df_tiab,
        `Seed DF (Title)` = .data$df_title,
        `PubMed BG Hits` = .data$bg_hits,
        Lift = round(.data$lift, 1),
        Score = round(.data$term_score, 3)
      )
    DT::datatable(disp_t, rownames = FALSE, selection = "none", options = list(pageLength = 20, scrollX = TRUE)) |>
      DT::formatCurrency(c("PubMed BG Hits", "Lift"), currency = "", digits = 0)
  })

  output$exclusion_blocks_ui <- shiny::renderUI({
    res <- current_result()
    if (is.null(res)) return(NULL)
    eh <- res$exclusion_hints %||% .empty_exclusion_hints()
    eb_list <- attr(eh, "exclusion_blocks") %||% list()
    if (length(eb_list) == 0L) return(NULL)

    eb_rows <- lapply(eb_list, function(eb) {
      tags$div(
        class = "mb-2",
        tags$strong(sprintf("%s (~%.1f%% of sampled noise): ", eb$name, eb$noise_covered_pct)),
        lapply(eb$terms, function(sp) tags$span(class = "excl-pill", sp))
      )
    })
    tags$div(
      class = "p-3 bg-light border rounded mb-2",
      tags$h6("Exclusion Concept Blocks (Greedy Set-Cover, 0 Seeds Lost)", class = "fw-bold mb-2"),
      eb_rows
    )
  })

  output$exclusion_hints_table <- DT::renderDT({
    res <- current_result()
    req(res)
    eh <- res$exclusion_hints %||% .empty_exclusion_hints()
    if (nrow(eh) == 0L) {
      return(DT::datatable(
        tibble::tibble(Message = "No significant exclusion hints found for this query."),
        rownames = FALSE, selection = "none", options = list(dom = "t")
      ))
    }
    if (!("concept" %in% names(eh))) eh$concept <- "general_noise"
    disp_eh <- eh |>
      dplyr::transmute(
        Concept = .data$concept,
        Term = .data$term,
        Scope = .data$field,
        `Exclude Spec` = .data$exclude_spec,
        `Noise %` = .data$noise_pct,
        `Hits Removed` = .data$exact_hits_removed,
        `% Reduced` = .data$pct_hits_reduced,
        `Seeds Lost` = .data$seeds_lost_n,
        Safety = .data$safety_tier,
        `Missed Citekeys` = .data$seeds_lost_citekeys,
        `PubMed Clause` = .data$pm_clause,
        `GS Clause` = .data$gs_clause
      )
    DT::datatable(
      disp_eh,
      rownames = FALSE,
      selection = "multiple",
      options = list(pageLength = 25, lengthMenu = c(15, 25, 50), scrollX = TRUE)
    ) |>
      DT::formatCurrency("Hits Removed", currency = "", digits = 0) |>
      DT::formatStyle(
        "Safety",
        backgroundColor = DT::styleEqual("Safe (0 seeds lost)", "#dcfce7"),
        color = DT::styleEqual("Safe (0 seeds lost)", "#166534"),
        fontWeight = "bold"
      )
  })

  .append_exclusions_to_sidebar <- function(new_specs) {
    new_specs <- unique(new_specs[nzchar(new_specs)])
    if (length(new_specs) == 0L) {
      shiny::showNotification("No exclusion terms selected.", type = "warning")
      return()
    }
    cur <- parse_term_list(input$exclude_terms)
    merged <- unique(c(cur, new_specs))
    shiny::updateTextAreaInput(session, "exclude_terms", value = paste(merged, collapse = ", "))
    shiny::showNotification(
      sprintf("Added %d exclusion(s) to sidebar. Click 'Run Query Optimizer' to apply.", length(new_specs)),
      type = "message", duration = 4
    )
  }

  shiny::observeEvent(input$add_selected_excl, {
    res <- current_result()
    eh <- if (!is.null(res)) res$exclusion_hints else NULL
    idx <- input$exclusion_hints_table_rows_selected
    if (is.null(eh) || nrow(eh) == 0L || length(idx) == 0L) {
      shiny::showNotification("Select one or more rows in the Exclusion Hints table first.", type = "warning")
      return()
    }
    .append_exclusions_to_sidebar(eh$exclude_spec[idx])
  })

  shiny::observeEvent(input$add_setcover_excl, {
    res <- current_result()
    eh <- if (!is.null(res)) res$exclusion_hints else NULL
    if (is.null(eh) || nrow(eh) == 0L) {
      shiny::showNotification("Run the optimizer first to discover exclusion hints.", type = "warning")
      return()
    }
    sc_specs <- attr(eh, "setcover_ex") %||% utils::head(eh$exclude_spec[eh$seeds_lost_n == 0L], 10L)
    .append_exclusions_to_sidebar(sc_specs)
  })

  shiny::observeEvent(input$add_safe_excl, {
    res <- current_result()
    eh <- if (!is.null(res)) res$exclusion_hints else NULL
    if (is.null(eh) || nrow(eh) == 0L) {
      shiny::showNotification("Run the optimizer first to discover exclusion hints.", type = "warning")
      return()
    }
    safe_specs <- eh$exclude_spec[eh$seeds_lost_n == 0L]
    .append_exclusions_to_sidebar(safe_specs)
  })

  output$all_eval_table <- DT::renderDT({
    res <- current_result()
    req(res, nrow(res$eval_df) > 0L)
    ord_df <- res$eval_df |>
      dplyr::arrange(dplyr::desc(.data$is_pareto), dplyr::desc(.data$combined_recall), .data$pm_hits)
    disp_e <- ord_df |>
      dplyr::transmute(
        ID = .data$eval_id,
        Type = ifelse(.data$is_custom, ifelse(nzchar(.data$custom_label), .data$custom_label, "Custom"), "Auto"),
        Pareto = ifelse(.data$is_pareto, "Yes", "No"),
        Stage = .data$stage,
        `Recall %` = round(100 * .data$combined_recall, 1),
        Found = paste0(.data$seeds_found, "/", .data$seeds_resolvable),
        `PM Hits` = .data$pm_hits,
        `GS Chars` = .data$gs_chars,
        `PubMed Query` = .data$pubmed_query,
        `Google Scholar Query` = .data$google_scholar_query
      )
    sel_idx <- which(ord_df$eval_id == (selected_eval_id() %||% -1L))
    DT::datatable(
      disp_e, rownames = FALSE,
      selection = list(mode = "single", selected = if (length(sel_idx) == 1L) sel_idx else NULL),
      options = list(pageLength = 20, scrollX = TRUE)
    ) |>
      DT::formatCurrency("PM Hits", currency = "", digits = 0)
  })

  output$report_preview_ui <- shiny::renderUI({
    res <- current_result()
    if (is.null(res)) {
      return(tags$div(class = "text-muted small", "Run the optimizer to generate and preview the Markdown report and export bundle."))
    }
    tmp_dir <- tempfile("shiny_rep_prev_")
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    files <- write_optimization_report(res, output_dir = tmp_dir)
    md_txt <- paste(readLines(files$report_md, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    tags$div(
      class = "border rounded p-3 bg-white",
      HTML(markdown::markdownToHTML(text = md_txt, fragment.only = TRUE))
    )
  })

  output$dl_zip <- shiny::downloadHandler(
    filename = function() paste0("query_optimizer_report_", format(Sys.Date(), "%Y%m%d"), ".zip"),
    content = function(file) {
      res <- current_result()
      req(res)
      tmp_dir <- tempfile("shiny_rep_zip_")
      on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
      write_optimization_report(res, output_dir = tmp_dir)
      zip::zipr(zipfile = file, files = list.files(tmp_dir, full.names = TRUE))
    }
  )

  output$dl_pareto_csv <- shiny::downloadHandler(
    filename = function() paste0("pareto_queries_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      res <- current_result()
      req(res)
      out <- res$pareto_df |>
        dplyr::select(
          "rank", "is_recommended", "stage",
          "combined_recall", "seeds_found", "seeds_resolvable",
          "pm_recall", "pm_hits",
          "gs_chars", "n_blocks", "n_terms",
          "pubmed_query", "google_scholar_query",
          "missed_citekeys", "pubmed_url", "google_scholar_url"
        )
      readr::write_csv(out, file)
    }
  )

  output$dl_config_yml <- shiny::downloadHandler(
    filename = function() "optimizer_config.yml",
    content = function(file) {
      cfg_res <- build_config_from_inputs(
        required_txt = input$required_terms,
        banned_txt = input$banned_terms,
        exclude_txt = input$exclude_terms,
        custom_blocks_yaml = input$custom_blocks_yaml,
        max_blocks = input$max_blocks,
        max_terms_per_block = input$max_terms_per_block,
        max_gs_chars = input$max_gs_chars
      )
      cfg_out <- if (isTRUE(cfg_res$ok)) cfg_res$config else list()
      writeLines(yaml::as.yaml(cfg_out), file, useBytes = TRUE)
    }
  )
  }
}

#' Create the revlitquery Shiny application object
#'
#' @param initial_result Optional precomputed result list from [optimize_search_queries()] to display immediately on startup.
#' @param initial_bib_label Optional display label for `initial_result`.
#' @return A `shiny.appobj` object that can be passed to [shiny::runApp()].
#' @export
revlitquery_app <- function(initial_result = NULL, initial_bib_label = NULL) {
  shiny::shinyApp(
    ui = .revlitquery_ui(),
    server = .revlitquery_make_server(initial_result = initial_result, initial_bib_label = initial_bib_label)
  )
}

#' Launch the interactive revlitquery Shiny application
#'
#' @param initial_result Optional precomputed result list from [optimize_search_queries()] to display immediately on startup.
#' @param initial_bib_label Optional display label for `initial_result`.
#' @param launch.browser Logical; open the application in the default browser (default `TRUE`).
#' @param ... Additional arguments passed to [shiny::runApp()].
#' @return Invisibly returns the result of [shiny::runApp()].
#' @export
run_revlitquery_app <- function(initial_result = NULL, initial_bib_label = NULL, launch.browser = TRUE, ...) {
  shiny::runApp(
    revlitquery_app(initial_result = initial_result, initial_bib_label = initial_bib_label),
    launch.browser = launch.browser,
    ...
  )
}

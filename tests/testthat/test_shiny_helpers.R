test_that("Shiny term list and custom block YAML helpers parse and round-trip cleanly", {
  expect_equal(parse_term_list(NULL), character(0))
  expect_equal(parse_term_list("  H5N1, case fatality; Avian Influenza\nh5n1\n "),
               c("h5n1", "case fatality", "avian influenza"))

  empty_cb <- parse_custom_blocks_yaml("   ")
  expect_true(empty_cb$ok)
  expect_null(empty_cb$blocks)

  good_cb <- parse_custom_blocks_yaml("pathogen:\n  - H5N1\n  - avian influenza\noutcome:\n  - case fatality\n")
  expect_true(good_cb$ok)
  expect_equal(names(good_cb$blocks), c("pathogen", "outcome"))
  expect_equal(good_cb$blocks$pathogen, c("h5n1", "avian influenza"))

  bad_cb <- parse_custom_blocks_yaml("pathogen: [unclosed")
  expect_false(bad_cb$ok)
  expect_match(bad_cb$error, "Invalid YAML")

  unnamed_cb <- parse_custom_blocks_yaml("- just_a_list")
  expect_false(unnamed_cb$ok)

  yml_back <- blocks_to_yaml(list(
    list(name = "concept_1", terms = c("h5n1", "avian influenza")),
    list(name = "concept_2", terms = c("case fatality", "mortality"))
  ))
  reparsed <- parse_custom_blocks_yaml(yml_back)
  expect_true(reparsed$ok)
  expect_equal(reparsed$blocks$concept_1, c("h5n1", "avian influenza"))
  expect_equal(reparsed$blocks$concept_2, c("case fatality", "mortality"))

  cfg_res <- build_config_from_inputs(
    required_txt = "h5n1",
    banned_txt = "study, paper",
    exclude_txt = "vaccine",
    custom_blocks_yaml = yml_back,
    max_blocks = 3, max_terms_per_block = 4, max_gs_chars = 220
  )
  expect_true(cfg_res$ok)
  expect_equal(cfg_res$config$required_terms, "h5n1")
  expect_equal(cfg_res$config$banned_terms, c("study", "paper"))
  expect_equal(cfg_res$config$exclude_terms, "vaccine")
  expect_equal(cfg_res$config$max_gs_chars, 220L)
  expect_equal(names(cfg_res$config$custom_blocks), c("concept_1", "concept_2"))
})

test_that("recompute_recommendation handles all constraint regimes and empty tables", {
  empty_res <- recompute_recommendation(tibble::tibble(), 0.8, Inf)
  expect_equal(empty_res$recommended_status, "empty")

  p_df <- tibble::tibble(
    rank = 1:3,
    signature = c("s1", "s2", "s3"),
    combined_recall = c(1.0, 0.8, 0.5),
    pm_hits = c(800L, 120L, 20L),
    is_recommended = c(TRUE, FALSE, FALSE)
  )

  r1 <- recompute_recommendation(p_df, min_recall = 0.75, max_hits = 200)
  expect_equal(r1$recommended_status, "meets_constraints")
  expect_equal(r1$recommended_sig, "s2")
  expect_equal(r1$pareto_df$is_recommended, c(FALSE, TRUE, FALSE))

  r2 <- recompute_recommendation(p_df, min_recall = 0.75, max_hits = 50)
  expect_equal(r2$recommended_status, "exceeds_max_hits")
  expect_equal(r2$recommended_sig, "s2")

  r3 <- recompute_recommendation(p_df, min_recall = 1.0, max_hits = Inf)
  expect_equal(r3$recommended_status, "meets_constraints")
  expect_equal(r3$recommended_sig, "s1")

  p_low <- p_df[3, ]
  r4 <- recompute_recommendation(p_low, min_recall = 0.8, max_hits = Inf)
  expect_equal(r4$recommended_status, "below_min_recall")
  expect_equal(r4$recommended_sig, "s3")
})

test_that("build_pareto_plotly and build_coverage_display_df handle empty and populated data", {
  p_empty <- build_pareto_plotly(tibble::tibble(), tibble::tibble())
  expect_s3_class(p_empty, "plotly")

  eval_df <- tibble::tibble(
    eval_id = 1:4,
    stage = c("stage1", "stage1", "stage2_r1", "custom"),
    is_custom = c(FALSE, FALSE, FALSE, TRUE),
    custom_label = c("", "", "", "Original"),
    is_pareto = c(TRUE, TRUE, FALSE, FALSE),
    rank = c(1L, 2L, NA_integer_, NA_integer_),
    is_recommended = c(TRUE, FALSE, FALSE, FALSE),
    combined_recall = c(1.0, 0.75, 0.75, 0.50),
    seeds_found = c(4L, 3L, 3L, 2L),
    seeds_resolvable = c(4L, 4L, 4L, 4L),
    pm_hits = c(0L, 150L, 400L, 650L),
    gs_chars = c(45L, 62L, 80L, 55L),
    google_scholar_query = c("h5n1 fatality", "h5n1 mortality", "h5n1 human mortality", "h5n1 influenza"),
    found_idx = list(c(1L, 2L, 3L), c(1L, 2L), c(1L, 2L), c(1L))
  )
  pareto_df <- eval_df[eval_df$is_pareto, ]

  p_full <- build_pareto_plotly(eval_df, pareto_df, selected_rank = 1L, min_recall = 0.75, max_hits = 250, selected_eval_id = 4L)
  expect_s3_class(p_full, "plotly")

  seeds_df <- tibble::tibble(
    seed_id = 1:3,
    citekey = c("a1", "b2", "c3"),
    year = c("2008", "2012", "2016"),
    full_title = c("Title 1", "Title 2", "Title 3"),
    in_pubmed = c(TRUE, TRUE, FALSE),
    pmid = c("111", "222", "")
  )
  pareto_cov <- tibble::tibble(
    rank = c(1L, 2L),
    found_idx = list(c(1L, 2L), c(1L))
  )
  cov_tbl <- build_coverage_display_df(seeds_df, pareto_cov, eval_df)
  expect_equal(names(cov_tbl), c("Citekey", "Year", "Title", "PubMed", "Q1", "Q2", "Custom: Original"))
  expect_equal(cov_tbl$Q1, c("Found", "Found", "N/A"))
  expect_equal(cov_tbl$Q2, c("Found", "Missed", "N/A"))
  expect_equal(cov_tbl[["Custom: Original"]], c("Found", "Missed", "N/A"))
})

test_that("revlitquery_app loads cleanly and sample BibTeX action populates seeds_raw in testServer", {
  old_opt <- options(revlitquery.min_interval = 0)
  on.exit(options(old_opt), add = TRUE)
  expect_true(parse_custom_blocks_yaml(NULL)$ok)
  expect_false(parse_custom_blocks_yaml("a: []")$ok)
  expect_false(parse_custom_blocks_yaml("a:\n  - '   '")$ok)
  expect_false(build_config_from_inputs(custom_blocks_yaml = "- invalid")$ok)

  app <- revlitquery_app()
  expect_s3_class(app, "shiny.appobj")

  # End-to-end test of optimize_search_queries() and preloaded revlitquery_app(initial_result = res)
  demo_bib <- .write_demo_bibtex()
  cfg_file <- tempfile(fileext = ".yml")
  out_dir <- tempfile("pkg_api_out_")
  on.exit(unlink(c(cfg_file, out_dir), recursive = TRUE), add = TRUE)
  writeLines("max_gs_chars: 250\n", cfg_file)

  xml_demo <- charToRaw(paste0(
    '<?xml version="1.0"?><PubmedArticleSet>',
    paste(vapply(c("18439345", "22362880", "21679679", "18959499", "27211899", "23803488", paste0("900", 1:8)), function(id) {
      sprintf(
        '<PubmedArticle><MedlineCitation><PMID>%s</PMID><Article><ArticleTitle>Avian influenza H5N1 case fatality and vaccine mice H7N9</ArticleTitle><Abstract><AbstractText>Human H5N1 mortality and case fatality rate in patients.</AbstractText></Abstract></Article></MedlineCitation></PubmedArticle>',
        id
      )
    }, character(1)), collapse = ""),
    '</PubmedArticleSet>'
  ))
  es_demo <- charToRaw(paste0(
    '{"esearchresult":{"count":"120","idlist":["18439345","22362880","21679679","18959499","27211899","23803488",',
    paste(sprintf('"900%d"', 1:8), collapse = ","),
    ']}}'
  ))

  res_api <- httr2::with_mocked_responses(
    function(req) {
      if (grepl("efetch", req$url)) {
        httr2::response(status_code = 200L, headers = list(`Content-Type` = "text/xml"), body = xml_demo)
      } else {
        httr2::response(status_code = 200L, headers = list(`Content-Type` = "application/json"), body = es_demo)
      }
    },
    optimize_search_queries(
      bib_path = demo_bib,
      output_dir = out_dir,
      config_path = cfg_file,
      min_recall = 0.75,
      max_hits = 1000,
      max_evals = 12L,
      use_cache = FALSE,
      custom_queries = "h5n1[tiab] AND fatal*[tiab]"
    )
  )
  expect_true(nrow(res_api$pareto_df) >= 1L)
  expect_true(file.exists(file.path(out_dir, "report.md")))

  # Test Shiny server with initial_result and full interactive flow
  app_pre <- revlitquery_app(initial_result = res_api, initial_bib_label = "demo.bib")
  httr2::with_mocked_responses(
    function(req) {
      if (grepl("efetch", req$url)) {
        httr2::response(status_code = 200L, headers = list(`Content-Type` = "text/xml"), body = xml_demo)
      } else {
        httr2::response(status_code = 200L, headers = list(`Content-Type` = "application/json"), body = es_demo)
      }
    },
    shiny::testServer(app_pre, {
      session$setInputs(
        min_recall = 75,
        max_hits = 1000,
        max_gs_chars = 256,
        max_evals = 12,
        max_blocks = 3,
        max_terms_per_block = 4,
        use_cache = FALSE,
        email = "",
        ncbi_key = "",
        required_terms = "",
        banned_terms = "",
        exclude_terms = "",
        custom_blocks_yaml = ""
      )

      # Verify preloaded state & outputs render cleanly
      expect_equal(bib_label(), "demo.bib")
      expect_equal(output$vb_seeds, "6")
      expect_match(output$vb_seeds_sub, "in PubMed")
      expect_match(output$vb_recall, "%")
      expect_match(output$vb_recall_sub, "PubMed seeds")
      expect_true(nzchar(output$vb_pm_hits))
      expect_match(output$vb_pm_hits_sub, "cap")
      expect_match(output$vb_gs_chars, "256")
      expect_match(output$vb_gs_sub, "Google Scholar")
      expect_false(is.null(output$bib_status_ui))
      expect_false(is.null(output$selected_query_card))
      expect_false(is.null(output$blocks_summary_ui))
      expect_false(is.null(output$exclusion_blocks_ui))
      expect_false(is.null(output$report_preview_ui))
      expect_true(nzchar(output$pareto_plot))
      expect_true(nzchar(output$pareto_table))
      expect_true(nzchar(output$coverage_table))
      expect_true(nzchar(output$terms_table))
      expect_true(nzchar(output$exclusion_hints_table))
      expect_true(nzchar(output$all_eval_table))

      # Trigger table selections, copy blocks, and exclusion quick-add buttons
      session$setInputs(pareto_table_rows_selected = 1L)
      session$setInputs(all_eval_table_rows_selected = nrow(res_api$eval_df))
      expect_false(is.null(output$selected_query_card))
      expect_true(nzchar(output$vb_recall_sub))
      session$setInputs(max_hits = 10)
      expect_match(output$vb_pm_hits_sub, "Exceeds")
      session$setInputs(max_hits = 0)
      expect_match(output$vb_pm_hits_sub, "No hit cap")
      session$setInputs(max_hits = 1000)
      session$setInputs(copy_blocks_to_editor = 1L)
      session$setInputs(copy_blocks_inline = 1L)
      session$setInputs(exclusion_hints_table_rows_selected = 1L, add_selected_excl = 1L)
      session$setInputs(add_setcover_excl = 1L)
      session$setInputs(add_safe_excl = 1L)

      # Trigger downloadHandlers (dl_pareto_csv, dl_config_yml, dl_zip)
      expect_true(file.exists(output$dl_pareto_csv))
      expect_true(file.exists(output$dl_config_yml))
      expect_true(file.exists(output$dl_zip))

      # Trigger load_sample, file upload, custom query evaluation (when opt_result is NULL), and run_opt (preserving custom query)
      session$setInputs(load_sample = 1L)
      expect_equal(bib_label(), "sample_h5n1_cfr.bib")
      expect_equal(nrow(seeds_raw()), 6L)
      # Pre-run coverage_table shows seed preview when opt_result is NULL
      expect_true(nzchar(output$coverage_table))
      expect_equal(output$vb_seeds, "6")
      expect_match(output$vb_seeds_sub, "Parsed")
      expect_equal(output$vb_recall, "--")
      expect_equal(output$vb_pm_hits, "--")
      expect_equal(output$vb_gs_chars, "--")

      session$setInputs(bib_file = list(datapath = demo_bib, name = "uploaded.bib"))
      expect_equal(bib_label(), "uploaded.bib")

      session$setInputs(
        custom_pm_query = "h5n1[tiab] AND mortality[tiab]",
        custom_query_label = "Live custom",
        custom_gs_query = "",
        add_custom_query = 1L
      )
      expect_false(is.null(opt_result()))
      expect_true(any(opt_result()$eval_df$is_custom))
      expect_false(is.null(output$selected_query_card))

      session$setInputs(run_opt = 1L)
      expect_true(nrow(opt_result()$pareto_df) >= 1L)
      expect_true(any(opt_result()$eval_df$is_custom))
    })
  )

  # Also test empty initial Shiny state (before loading any BibTeX)
  shiny::testServer(app, {
    session$setInputs(min_recall = 75, max_hits = 2000, max_gs_chars = 256)
    expect_equal(output$vb_seeds, "--")
    expect_match(output$vb_seeds_sub, "Load .bib")
    expect_false(is.null(output$bib_status_ui))
    expect_false(is.null(output$selected_query_card))
    expect_false(is.null(output$blocks_summary_ui))
    expect_false(is.null(output$report_preview_ui))
    expect_null(output$exclusion_blocks_ui)
    expect_true(nzchar(output$pareto_plot))
    session$setInputs(run_opt = 1L)
    session$setInputs(add_custom_query = 1L)
    session$setInputs(copy_blocks_to_editor = 1L)
    session$setInputs(add_selected_excl = 1L)
    session$setInputs(add_setcover_excl = 1L)
    session$setInputs(add_safe_excl = 1L)
  })
})

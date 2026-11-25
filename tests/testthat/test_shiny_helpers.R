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
  app <- revlitquery_app()
  expect_s3_class(app, "shiny.appobj")

  shiny::testServer(app, {
    expect_null(seeds_raw())
    session$setInputs(load_sample = 1)
    df <- seeds_raw()
    expect_s3_class(df, "tbl_df")
    expect_equal(nrow(df), 6L)
    expect_equal(bib_label(), "sample_h5n1_cfr.bib")
    session$setInputs(custom_pm_query = "", add_custom_query = 1)
    expect_null(opt_result())
  })
})

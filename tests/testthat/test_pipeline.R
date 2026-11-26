test_that("BibTeX parser handles nested braces, multiline abstracts, LaTeX, DOIs, PMIDs, and skips comments", {
  tmp_bib <- tempfile(fileext = ".bib")
  on.exit(unlink(tmp_bib), add = TRUE)

  bib_txt <- paste(
    "@comment{This should be skipped}",
    "@string{EID = {Emerging Infectious Diseases}}",
    "@ARTICLE{smith2020_h5n1,",
    "  title = {Case fatality of {H5N1} in {S}outheastern {A}sia: a \\textit{review} & update},",
    "  author = {Mu{\\\"o}ller, Hans and Garc{\\'i}a, Ana},",
    "  journal = EID,",
    "  year = 2020,",
    "  doi = {https://doi.org/10.1000/XYZ.123},",
    "  url = {https://pubmed.ncbi.nlm.nih.gov/32109876/},",
    "  keywords = {H5N1; Case fatality rate, Avian influenza},",
    "  abstract = \"Line one about case-fatality risk",
    "  and line two with {nested {braces}} and 50\\% mortality.\"",
    "}",
    "@inproceedings{doe2021_cfr,",
    "  title = \"Estimates of {CFR} in human influenza\",",
    "  pmid = {99887766},",
    "  doi = {doi:10.5555/TEST.99},",
    "  abstract = {Human H5N1 mortality estimates across outbreaks.}",
    "}",
    sep = "\n"
  )
  writeLines(bib_txt, tmp_bib, useBytes = TRUE)

  df <- parse_bibtex(tmp_bib)
  expect_equal(nrow(df), 2L)
  expect_equal(df$citekey, c("smith2020_h5n1", "doe2021_cfr"))
  expect_equal(df$entry_type, c("article", "inproceedings"))
  expect_match(df$title[1], "Case fatality of H5N1 in Southeastern Asia: a review & update", fixed = TRUE)
  expect_match(df$abstract[1], "nested braces", fixed = TRUE)
  expect_match(df$abstract[1], "50% mortality", fixed = TRUE)
  expect_equal(df$doi, c("10.1000/xyz.123", "10.5555/test.99"))
  expect_equal(df$pmid, c("32109876", "99887766"))
  expect_match(df$keywords[1], "h5n1; case fatality rate; avian influenza", fixed = TRUE)
  expect_equal(df$year[1], "2020")
})

test_that("Term extraction preserves acronyms, normalizes hyphens, respects config, and builds stem-disjoint blocks", {
  seeds <- tibble::tibble(
    seed_id = 1:4,
    citekey = paste0("p", 1:4),
    full_title = c(
      "Case fatality rate of human H5N1 avian influenza",
      "Mortality and case-fatality risk in H5N1 patients",
      "Estimates of human H5N1 case fatality",
      "Global epidemiology and case fatality of avian influenza H5N1"
    ),
    full_abstract = c(
      "We estimated the case fatality rate among confirmed human H5N1 infections.",
      "High case-fatality ratio was observed in H5N1 avian influenza cases.",
      "Human H5N1 avian influenza continues to cause fatal infections.",
      "Pooled mortality and case fatality risk for H5N1 virus infections."
    ),
    full_text = character(4),
    pmid = as.character(1001:1004),
    oa_id = paste0("W", 1001:1004),
    in_pubmed = TRUE,
    in_openalex = TRUE,
    mesh_terms = rep(list(character(0)), 4),
    all_keywords = rep(list(character(0)), 4)
  )
  seeds$full_text <- paste(seeds$full_title, seeds$full_abstract)

  mock_bg <- function(term) {
    bg_map <- c(
      "h5n1" = 9000L,
      "avian influenza" = 18000L,
      "case fatality" = 25000L,
      "case fatality rate" = 14000L,
      "case fatality risk" = 1800L,
      "mortality" = 1200000L,
      "human" = 15000000L
    )
    if (term %in% names(bg_map)) bg_map[[term]] else 80000L
  }

  terms_df <- extract_candidate_terms(
    seeds,
    config = list(banned_terms = c("infections")),
    probe_pubmed = TRUE,
    eval_bg_fn = mock_bg
  )

  expect_true("h5n1" %in% terms_df$term)
  expect_true("case fatality" %in% terms_df$term)
  expect_false("infections" %in% terms_df$term)
  expect_false("study" %in% terms_df$term)

  cf_row <- terms_df[terms_df$term == "case fatality", ]
  expect_equal(cf_row$df_tiab, 4L) # matches both "case fatality" and "case-fatality"

  blocks <- build_concept_blocks(terms_df, seeds, max_blocks = 3L, max_terms_per_block = 4L)
  expect_gte(length(blocks), 2L)
  # Blocks should be stem-disjoint across anchors
  anchors <- vapply(blocks, function(b) b$terms[1L], character(1))
  expect_false(shares_stem(anchors[1L], anchors[2L]))

  # Custom blocks override
  custom <- build_concept_blocks(
    terms_df, seeds,
    config = list(custom_blocks = list(virus = c("h5n1", "avian influenza"), outcome = c("case fatality", "mortality")))
  )
  expect_equal(length(custom), 2L)
  expect_equal(names(custom[[1]]), c("name", "terms"))
  expect_equal(custom[[1]]$name, "virus")
  expect_equal(custom[[2]]$terms, c("case fatality", "mortality"))
})

test_that("Query translator builds semantically aligned PubMed and Google Scholar (<=256) queries", {
  q <- make_query(
    blocks = list(
      list(name = "pathogen", field = "title", terms = c("H5N1", "avian influenza")),
      list(name = "outcome", field = "tiab", terms = c("case fatality rate", "mortality"))
    ),
    exclude_terms = c("vaccine trial", "poultry")
  )

  pm <- to_pubmed_query(q)
  gs <- to_google_scholar_query(q)

  expect_match(pm, '\\("h5n1"\\[ti\\] OR "avian influenza"\\[ti\\]\\) AND \\("case fatality rate"\\[tiab\\] OR "mortality"\\[tiab\\]\\)')
  expect_match(pm, 'NOT \\("vaccine trial"\\[tiab\\] OR "poultry"\\[tiab\\]\\)')

  # Check Google Scholar syntax: implicit AND, uppercase OR, intitle:, unquoted single token, -exclude
  expect_match(gs, 'intitle:h5n1', fixed = TRUE)
  expect_match(gs, 'intitle:"avian influenza"', fixed = TRUE)
  expect_match(gs, '("case fatality rate" OR mortality)', fixed = TRUE)
  expect_match(gs, '-"vaccine trial" -poultry', fixed = TRUE)
  expect_true(is_valid_gs_query(q, max_gs_chars = 256L))

  # Check URL helpers
  expect_true(startsWith(pubmed_url(pm), "https://pubmed.ncbi.nlm.nih.gov/?term="))
  expect_true(startsWith(google_scholar_url(gs), "https://scholar.google.com/scholar?q="))

  # Overlong Google Scholar query fails validation
  long_q <- make_query(list(
    list(name = "b1", field = "title", terms = paste0("longtermphrase", 1:15))
  ))
  expect_false(is_valid_gs_query(long_q, max_gs_chars = 256L))

  # Redundant superstring inside the same OR block is simplified
  dup_q <- make_query(list(list(name = "b1", field = "tiab", terms = c("human h5n1", "h5n1"))))
  expect_equal(dup_q$blocks[[1]]$terms, "h5n1")
})

test_that("Pareto optimizer refines queries and writes all 6 report files", {
  seeds <- tibble::tibble(
    seed_id = 1:5,
    citekey = paste0("seed", 1:5),
    year = as.character(2008:2012),
    doi = paste0("10.1000/s", 1:5),
    pmid = c("101", "102", "103", "104", ""),
    in_pubmed = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    full_title = c(
      "Case fatality rate of human H5N1 avian influenza",
      "Case fatality risk of human H5N1 infections",
      "Mortality in human H5N1 avian influenza cases",
      "Fatal outcomes in avian influenza A H5N1",
      "Serological case fatality study of H5N1 in Asia"
    ),
    full_abstract = c(
      "Confirmed human H5N1 infections with high case fatality rate.",
      "Estimated case fatality risk for H5N1 influenza.",
      "High mortality in confirmed H5N1 patients.",
      "Avian influenza H5N1 fatal human cases.",
      "Non-PubMed journal study on H5N1 case fatality."
    ),
    full_text = character(5),
    mesh_terms = rep(list(character(0)), 5),
    all_keywords = rep(list(character(0)), 5)
  )
  seeds$full_text <- paste(seeds$full_title, seeds$full_abstract)

  mock_bg <- function(term) {
    m <- c(
      "h5n1" = 9000L, "avian influenza" = 18000L,
      "case fatality" = 25000L, "case fatality rate" = 14000L,
      "case fatality risk" = 1500L, "mortality" = 1200000L,
      "fatal" = 400000L, "human" = 15000000L
    )
    if (term %in% names(m)) m[[term]] else 50000L
  }

  terms_df <- extract_candidate_terms(
    seeds,
    config = list(required_terms = c("h5n1")),
    probe_pubmed = TRUE,
    eval_bg_fn = mock_bg
  )
  blocks <- build_concept_blocks(terms_df, seeds, max_blocks = 3L, max_terms_per_block = 4L)

  eval_pm_mock <- function(pm_query, seed_pmids, ...) {
    has_ti <- grepl("\\[ti\\]", pm_query)
    has_mort <- grepl('"mortality"', pm_query, fixed = TRUE)
    has_fatal <- grepl('"fatal"', pm_query, fixed = TRUE)
    has_cf <- grepl('"case fatality', pm_query, fixed = TRUE)
    n_or <- lengths(regmatches(pm_query, gregexpr(" OR ", pm_query)))

    base_hits <- if (has_mort) 900L else if (has_fatal) 600L else if (has_cf) 180L else 1500L
    if (has_ti) base_hits <- as.integer(base_hits * 0.45)
    base_hits <- as.integer(base_hits + n_or * 35L)

    found <- character(0)
    for (idx in which(seeds$in_pubmed)) {
      and_parts <- strsplit(pm_query, " AND ", fixed = TRUE)[[1]]
      ok_all <- TRUE
      for (ap in and_parts) {
        ms <- stringr::str_match_all(ap, '"([^"]+)"\\[(ti|tiab)\\]')[[1]]
        if (nrow(ms) == 0L) next
        ok_blk <- FALSE
        for (r in seq_len(nrow(ms))) {
          tgt <- if (ms[r, 3] == "ti") seeds$full_title[idx] else seeds$full_text[idx]
          if (term_matches_text(ms[r, 2], tgt)) {
            ok_blk <- TRUE
            break
          }
        }
        if (!ok_blk) { ok_all <- FALSE; break }
      }
      if (ok_all) found <- c(found, seeds$pmid[idx])
    }

    list(
      pm_hits = base_hits,
      pm_found = found,
      pm_recall = length(found) / length(seed_pmids),
      ok = TRUE
    )
  }

  res <- optimize_queries(
    seeds_df = seeds,
    terms_df = terms_df,
    blocks = blocks,
    config = list(required_terms = c("h5n1")),
    min_recall = 0.80,
    max_hits = 500,
    max_evals = 25L,
    max_gs_chars = 256L,
    discover_exclusions = FALSE,
    eval_pubmed_fn = eval_pm_mock
  )

  expect_true(nrow(res$pareto_df) >= 1L)
  expect_true(all(res$pareto_df$gs_chars <= 256L))
  expect_true(all(res$pareto_df$pm_hits <= 500))
  expect_true(all(grepl("h5n1", res$pareto_df$pubmed_query, fixed = TRUE)))
  expect_equal(sum(res$pareto_df$is_recommended), 1L)

  out_dir <- tempfile("qopt_out_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  files <- write_optimization_report(res, output_dir = out_dir)
  for (fp in unlist(files)) {
    expect_true(file.exists(fp), info = paste("Missing output:", fp))
  }
  expect_true(file.exists(file.path(out_dir, "exclusion_hints.csv")))
})

test_that("Block selectivity rejects redundant concepts, prevents generic leakage, and exclusion hints reduce hits safely", {
  seeds <- tibble::tibble(
    seed_id = 1:4,
    citekey = paste0("s", 1:4),
    year = as.character(2008:2011),
    doi = paste0("10.1000/x", 1:4),
    pmid = as.character(501:504),
    oa_id = paste0("W", 501:504),
    in_pubmed = TRUE,
    in_openalex = FALSE,
    full_title = c(
      "Human H5N1 case fatality rate in Indonesia",
      "Case fatality of human H5N1 avian influenza",
      "Mortality of human H5N1 virus infection",
      "Fatal human H5N1 cases in Asia"
    ),
    full_abstract = c(
      "Case fatality rate was 80% among confirmed human H5N1 cases.",
      "High case fatality and mortality in H5N1 patients.",
      "Human H5N1 virus caused fatal pneumonia and mortality.",
      "Human H5N1 case fatality estimates in outbreak."
    ),
    full_text = character(4),
    mesh_terms = rep(list(character(0)), 4),
    all_keywords = rep(list(character(0)), 4)
  )
  seeds$full_text <- paste(seeds$full_title, seeds$full_abstract)

  mock_bg2 <- function(term) {
    m <- c(
      "h5n1" = 9000L,
      "avian influenza" = 15000L,
      "influenza" = 135000L,
      "case fatality" = 15000L,
      "fatality" = 30000L,
      "fatal" = 172000L,
      "mortality" = 1200000L,
      "died" = 305000L,
      "respiratory" = 690000L,
      "case" = 2800000L,
      "human" = 3400000L,
      "human cases" = 6300L,
      "patients" = 4000000L
    )
    if (term %in% names(m)) m[[term]] else 50000L
  }
  terms_df <- extract_candidate_terms(seeds, probe_pubmed = TRUE, eval_bg_fn = mock_bg2)

  # Joint selectivity mock: "h5n1 AND influenza" barely reduces hits (8600 / 9000 = 95%),
  # whereas outcome terms drop to 450 / 9000 (5%)
  mock_joint <- function(t1, t2) {
    if (t2 %in% c("influenza", "avian influenza")) return(8600L)
    if (t2 %in% c("case fatality", "fatality", "fatal", "mortality", "died")) return(450L)
    if (t2 %in% c("human", "human cases", "patients")) return(3400L)
    7000L
  }
  blocks <- build_concept_blocks(terms_df, seeds, max_blocks = 3L, probe_selectivity = TRUE, eval_joint_fn = mock_joint)
  expect_true(length(blocks) >= 2L)
  expect_false("influenza" %in% unlist(lapply(blocks, `[[`, "terms")))
  expect_false("avian influenza" %in% unlist(lapply(blocks[-1L], `[[`, "terms")))
  expect_false("case" %in% unlist(lapply(blocks, `[[`, "terms")))

  # Unsupervised term clustering and orthogonal block coherence checks
  tc_map <- attr(blocks, "term_clusters")
  expect_true(is.character(tc_map) && length(tc_map) == nrow(terms_df))
  blk_names <- vapply(blocks, `[[`, character(1), "name")
  expect_equal(blk_names[1L], "primary_topic")
  expect_true(all(startsWith(blk_names[-1L], "concept_")))

  # Stage 2b rescue uses discovered clusters & high-lift outlier rescue cleanly
  sec_blk <- blocks[[2L]]
  q_sec_only <- make_query(list(
    list(name = "primary_topic", field = "tiab", terms = "h5n1"),
    list(name = sec_blk$name, field = "tiab", terms = sec_blk$terms[1L])
  ))
  nbrs <- .generate_neighbors(
    pareto_queries = list(q_sec_only),
    pareto_missed_idx = list(c(3L, 4L)),
    blocks = blocks,
    terms_df = terms_df,
    seeds_df = seeds,
    config = list(),
    seen_sigs = query_signature(q_sec_only)
  )
  expect_true(length(nbrs) >= 1L)

  # OpenAlex completely removed from project signatures
  expect_false("use_openalex" %in% names(formals(enrich_seeds)))
  expect_false("use_openalex" %in% names(formals(optimize_queries)))
  expect_false(exists("to_openalex_filter", mode = "function"))

  # Field-aware exclusion translation
  q_ex <- make_query(
    list(list(name = "b1", field = "tiab", terms = "h5n1")),
    exclude_terms = c("h7n9", "title:poultry")
  )
  expect_match(to_pubmed_query(q_ex), '"h7n9"[tiab] OR "poultry"[ti]', fixed = TRUE)
  expect_match(to_google_scholar_query(q_ex), 'h5n1 -h7n9 -intitle:poultry', fixed = TRUE)

  # Mocked exclusion hint discovery (with concept classification and wildcard support)
  mock_esearch <- function(term, retmax = 0L, ...) {
    if (retmax > 0L) {
      return(list(count = 3000L, idlist = c("501", as.character(9001:9010)), ok = TRUE))
    }
    if (grepl('"vaccin', term, fixed = TRUE)) return(list(count = 1200L, idlist = character(0), ok = TRUE))
    if (grepl('"poultry"', term, fixed = TRUE)) return(list(count = 900L, idlist = character(0), ok = TRUE))
    list(count = 300L, idlist = character(0), ok = TRUE)
  }
  mock_efetch <- function(pmids, ...) {
    tibble::tibble(
      pmid = pmids,
      pm_title = paste("Poultry vaccine and vaccination trial study", pmids),
      pm_abstract = paste("H5N1 poultry vaccine and vaccinated chickens in H7N9 mice cells", pmids),
      pm_mesh = rep(list(character(0)), length(pmids)),
      pm_keywords = rep(list(character(0)), length(pmids))
    )
  }

  hints <- discover_exclusion_hints(
    base_query = make_query(list(list(name = "b1", field = "tiab", terms = "h5n1"))),
    seeds_df = seeds,
    base_pm_hits = 3000L,
    sample_n = 10L,
    top_verify = 5L,
    pm_esearch_fn = mock_esearch,
    pm_efetch_fn = mock_efetch
  )
  expect_true(nrow(hints) >= 2L)
  expect_true("concept" %in% names(hints))
  expect_true(any(hints$term %in% c("vaccine", "vaccin*")))
  expect_equal(hints$seeds_lost_n[hints$term %in% c("vaccine", "vaccin*")][1L], 0L)
  expect_equal(hints$exact_hits_removed[hints$term %in% c("vaccine", "vaccin*")][1L], 1200L)
  expect_true(length(attr(hints, "exclusion_blocks")) >= 1L)
  expect_true(all(startsWith(hints$concept, "topic_")))
})

test_that("Custom queries handle multiple quoted and bare [ti]/[Title] tags", {
  multi_ti <- '("h5n1"[ti] OR "avian flu"[Title]) AND fatal*[tiab] NOT ("wild"[ti] OR birds[ti] OR cattle[Title])'
  gs_multi <- pubmed_to_google_scholar_approx(multi_ti)
  expect_match(gs_multi, "intitle:h5n1", fixed = TRUE)
  expect_match(gs_multi, 'intitle:"avian flu"', fixed = TRUE)
  expect_match(gs_multi, "-intitle:wild", fixed = TRUE)
  expect_match(gs_multi, "-intitle:birds", fixed = TRUE)
  expect_match(gs_multi, "-intitle:cattle", fixed = TRUE)

  seeds <- tibble::tibble(
    seed_id = 1:2,
    citekey = c("k1", "k2"),
    year = c("2008", "2012"),
    doi = paste0("10.1000/k", 1:2),
    pmid = c("11", "12"),
    in_pubmed = c(TRUE, TRUE),
    full_title = c("H5N1 fatality", "H5N1 mortality"),
    full_abstract = c("H5N1 fatal", "H5N1 cases"),
    full_text = c("H5N1 fatality", "H5N1 mortality"),
    mesh_terms = rep(list(character(0)), 2),
    all_keywords = rep(list(character(0)), 2)
  )
  mock_pm <- function(pm_query, seed_pmids, ...) {
    list(pm_hits = 150L, pm_found = c("11", "12"), pm_recall = 1.0, ok = TRUE)
  }
  res_multi <- evaluate_custom_query(
    pm_query = multi_ti,
    seeds_df = seeds,
    label = "Multi title tag query",
    eval_pubmed_fn = mock_pm
  )
  expect_equal(res_multi$pm_hits, 150L)
  expect_equal(res_multi$seeds_found, 2L)
  expect_match(res_multi$google_scholar_query, 'intitle:"avian flu"', fixed = TRUE)
})

test_that("Custom queries can be translated, evaluated, added to results, and reported", {
  raw_pm <- '((H5N1[Title/Abstract] OR H5[Title/Abstract]) AND (Influenza[Title/Abstract]) AND (fatal*[Title/Abstract])) NOT (H7N9[Title/Abstract] OR H1N1[Title/Abstract] OR "poultry"[ti])'
  gs_approx <- pubmed_to_google_scholar_approx(raw_pm)
  expect_match(gs_approx, "(H5N1 OR H5)", fixed = TRUE)
  expect_match(gs_approx, "fatal", fixed = TRUE)
  expect_match(gs_approx, "-H7N9 -H1N1 -intitle:poultry", fixed = TRUE)
  expect_false(grepl("\\[Title/Abstract\\]|\\*", gs_approx))

  seeds <- tibble::tibble(
    seed_id = 1:3,
    citekey = c("k1", "k2", "k3"),
    year = c("2008", "2012", "2016"),
    doi = paste0("10.1000/k", 1:3),
    pmid = c("11", "12", "13"),
    in_pubmed = c(TRUE, TRUE, TRUE),
    full_title = c("H5N1 fatality", "H5N1 mortality", "H5N1 serology"),
    full_abstract = c("H5N1 fatal", "H5N1 cases", "H5N1 serosurvey"),
    full_text = c("H5N1 fatality", "H5N1 mortality", "H5N1 serology"),
    mesh_terms = rep(list(character(0)), 3),
    all_keywords = rep(list(character(0)), 3)
  )

  mock_pm_cust <- function(pm_query, seed_pmids, ...) {
    if (grepl("fatal", pm_query)) {
      return(list(pm_hits = 120L, pm_found = c("11", "12"), pm_recall = 2 / 3, ok = TRUE))
    }
    list(pm_hits = 300L, pm_found = c("11"), pm_recall = 1 / 3, ok = TRUE)
  }

  c_row <- evaluate_custom_query(
    pm_query = 'H5N1[tiab] AND fatal*[tiab]',
    seeds_df = seeds,
    label = "Manuscript original",
    eval_pubmed_fn = mock_pm_cust
  )
  expect_true(c_row$is_custom)
  expect_equal(c_row$custom_label, "Manuscript original")
  expect_equal(c_row$seeds_found, 2L)
  expect_equal(c_row$pm_hits, 120L)

  # Merge into minimal result (before full optimization)
  min_res <- list(
    eval_df = tibble::tibble(), pareto_df = tibble::tibble(),
    exclusion_hints = .empty_exclusion_hints(),
    min_recall = 0.6, max_hits = 2000, max_gs_chars = 256L,
    seeds_df = seeds, terms_df = tibble::tibble(), blocks = list()
  )
  merged <- add_custom_query_to_result(min_res, c_row)
  expect_equal(nrow(merged$eval_df), 1L)
  expect_equal(nrow(merged$pareto_df), 1L)
  expect_true(merged$pareto_df$is_custom[1L])

  # Merge a dominated second custom query and ensure report writes cleanly
  c_dom <- evaluate_custom_query(
    pm_query = 'H5N1[tiab] AND serology[tiab]',
    seeds_df = seeds,
    label = "Narrower serology",
    eval_pubmed_fn = mock_pm_cust
  )
  merged2 <- add_custom_query_to_result(merged, c_dom)
  expect_equal(nrow(merged2$eval_df), 2L)
  expect_equal(nrow(merged2$pareto_df), 1L) # second is dominated

  out_dir <- tempfile("qopt_cust_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  files <- write_optimization_report(merged2, output_dir = out_dir)
  cov_csv <- readr::read_csv(files$coverage_csv, show_col_types = FALSE)
  expect_true("Custom: Narrower serology" %in% names(cov_csv))
  rep_txt <- paste(readLines(files$report_md), collapse = "\n")
  expect_match(rep_txt, "Custom / Benchmark Queries", fixed = TRUE)
})

test_that("Prefix wildcards, subtype code-family synonyms, domain qualifiers, and sibling subtype exclusions work together", {
  expect_true(all(term_matches_text("fatal*", c("case fatality rate", "fatal outcome", "ended fatally"))))
  expect_true(.is_code_family("h5n1", "h5"))
  expect_false(.is_code_family("h5n1", "h7n9"))
  expect_false(.is_code_family("h5n1", "h5n6"))

  # Wildcard subsumption inside make_query and clean GS translation
  q_w <- make_query(
    list(
      list(name = "primary_topic", field = "tiab", terms = c("h5n1", "h5")),
      list(name = "domain_context", field = "tiab", terms = "influenza"),
      list(name = "concept_outcome", field = "tiab", terms = c("fatal*", "fatality", "case fatality"))
    ),
    exclude_terms = c("h7n9", "h1n1", "h1n1pdm09")
  )
  expect_equal(q_w$blocks[[3L]]$terms, "fatal*")
  expect_equal(
    to_pubmed_query(q_w),
    '(("h5n1"[tiab] OR "h5"[tiab]) AND ("influenza"[tiab]) AND ("fatal*"[tiab])) NOT ("h7n9"[tiab] OR "h1n1"[tiab] OR "h1n1pdm09"[tiab])'
  )
  expect_equal(
    to_google_scholar_query(q_w),
    "(h5n1 OR h5) influenza fatal -h7n9 -h1n1 -h1n1pdm09"
  )
})

test_that("BibTeX parser and translator edge cases (concatenation, malformed fields, empty queries) are handled", {
  expect_equal(extract_pmid_value("", "", "Note with PMID: 98765432 in text"), "98765432")
  expect_equal(extract_pmid_value("", "", "no pmid here"), "")
  expect_error(parse_bibtex(file.path(tempdir(), "does_not_exist_12345.bib")), "not found")

  empty_bib <- tempfile(fileext = ".bib")
  on.exit(unlink(empty_bib), add = TRUE)
  writeLines("Just plain text without any bib entries", empty_bib)
  expect_error(parse_bibtex(empty_bib), "No valid BibTeX entries")

  edge_bib <- tempfile(fileext = ".bib")
  on.exit(unlink(edge_bib), add = TRUE)
  writeLines(c(
    "@comment{ignored comment}",
    "@preamble{\"some preamble\"}",
    "@article{edge2024,,",
    "  $badtoken,",
    "  malformed_no_equals,",
    "  title = \"Part One \" # {Part Two} # \" Part Three\",",
    "  year = 2024,",
    "  author = {Doe, Jane},",
    "}"
  ), edge_bib)
  df_edge <- parse_bibtex(edge_bib)
  expect_equal(nrow(df_edge), 1L)
  expect_equal(df_edge$title[1L], "Part One Part Two Part Three")

  q_empty <- make_query(list())
  expect_equal(to_pubmed_query(q_empty), "")
  expect_equal(to_google_scholar_query(q_empty), "")
  expect_false(is_valid_gs_query(q_empty))
  expect_equal(pubmed_to_google_scholar_approx("   "), "")

  q_ex_tags <- make_query(
    list(list(name = "c1", field = "tiab", terms = "h5n1")),
    exclude_terms = c("ti:ferret", "tiab:mouse")
  )
  expect_match(to_pubmed_query(q_ex_tags), '"ferret"[ti]', fixed = TRUE)
  expect_match(to_pubmed_query(q_ex_tags), '"mouse"[tiab]', fixed = TRUE)
  expect_match(to_google_scholar_query(q_ex_tags), "-intitle:ferret -mouse", fixed = TRUE)
})

test_that("PubMed API clients (cache, esearch, efetch, eval_query) work with httr2 mocked responses", {
  old_opt <- options(revlitquery.min_interval = 0)
  on.exit({
    options(old_opt)
    .rate_state$pubmed_interval <- 0.45
  }, add = TRUE)
  tmp_cache <- tempfile("pm_cache_")
  on.exit(unlink(tmp_cache, recursive = TRUE), add = TRUE)

  ck <- cache_key("test_p", list(a = 1))
  expect_null(cache_get(tmp_cache, ck, use_cache = TRUE))
  expect_null(cache_get(tmp_cache, ck, use_cache = FALSE))
  cache_set(tmp_cache, ck, list(val = 99), use_cache = TRUE)
  expect_equal(cache_get(tmp_cache, ck, use_cache = TRUE)$val, 99)

  resp_429 <- httr2::response(status_code = 429L, headers = list(`Retry-After` = "1"))
  resp_503 <- httr2::response(status_code = 503L)
  expect_true(.is_short_transient(resp_429))
  expect_true(.is_short_transient(resp_503))

  # 1. pubmed_esearch success + cache reuse + retmax=0 reuse of retmax=2500 cache
  es_json_ok <- charToRaw(paste0(
    '{"esearchresult":{"count":"2","idlist":["111","222"],',
    '"warninglist":{"phrasesignored":[],"quotedphrasesnotfound":[]},',
    '"errorlist":{"phrasesnotfound":[]}}}'
  ))
  httr2::with_mocked_responses(
    function(req) {
      httr2::response(
        status_code = 200L,
        headers = list(`Content-Type` = "application/json"),
        body = es_json_ok
      )
    },
    {
      r1 <- pubmed_esearch("h5n1[tiab]", retmax = 2500L, email = "a@b.c", api_key = "k1", cache_dir = tmp_cache, use_cache = TRUE)
      expect_true(r1$ok)
      expect_equal(r1$count, 2L)
      expect_equal(r1$idlist, c("111", "222"))
      # Second call hits exact cache
      r1_cached <- pubmed_esearch("h5n1[tiab]", retmax = 2500L, cache_dir = tmp_cache, use_cache = TRUE)
      expect_equal(r1_cached$count, 2L)
      # retmax = 0L reuses cached retmax = 2500L entry
      r0_reused <- pubmed_esearch("h5n1[tiab]", retmax = 0L, cache_dir = tmp_cache, use_cache = TRUE)
      expect_true(r0_reused$ok)
      expect_equal(r0_reused$count, 2L)
    }
  )

  # pubmed_esearch API ERROR and non-numeric count and network error
  httr2::with_mocked_responses(
    function(req) {
      httr2::response(
        status_code = 200L,
        headers = list(`Content-Type` = "application/json"),
        body = charToRaw('{"esearchresult":{"ERROR":"Bad query"}}')
      )
    },
    {
      r_err <- pubmed_esearch("bad_q", use_cache = FALSE)
      expect_false(r_err$ok)
    }
  )
  httr2::with_mocked_responses(
    function(req) {
      httr2::response(
        status_code = 200L,
        headers = list(`Content-Type` = "application/json"),
        body = charToRaw('{"esearchresult":{"count":"not_a_number"}}')
      )
    },
    {
      r_nan <- pubmed_esearch("nan_q", use_cache = FALSE)
      expect_false(r_nan$ok)
    }
  )
  httr2::with_mocked_responses(
    function(req) stop("Simulated connection error"),
    {
      r_exc <- pubmed_esearch("exc_q", use_cache = FALSE)
      expect_false(r_exc$ok)
    }
  )

  # 2. pubmed_efetch XML parsing, batch cache, per-PMID cache, and batch-scan harvest
  expect_equal(nrow(pubmed_efetch(character(0), cache_dir = tmp_cache)), 0L)
  xml_body <- charToRaw(paste0(
    '<?xml version="1.0"?>',
    '<PubmedArticleSet>',
    '<PubmedArticle><MedlineCitation><PMID>111</PMID>',
    '<Article><ArticleTitle>H5N1 fatality study</ArticleTitle>',
    '<Abstract><AbstractText>Severe human H5N1 mortality.</AbstractText></Abstract></Article>',
    '<MeshHeadingList><MeshHeading><DescriptorName>Influenza A Virus, H5N1 Subtype</DescriptorName></MeshHeading></MeshHeadingList>',
    '<KeywordList><Keyword>H5N1</Keyword><Keyword>CFR</Keyword></KeywordList>',
    '</MedlineCitation></PubmedArticle>',
    '<PubmedArticle><MedlineCitation><PMID>222</PMID>',
    '<Article><ArticleTitle>Second H5N1 article</ArticleTitle>',
    '<Abstract><AbstractText>Case fatality risk in patients.</AbstractText></Abstract></Article>',
    '</MedlineCitation></PubmedArticle>',
    '</PubmedArticleSet>'
  ))
  httr2::with_mocked_responses(
    function(req) {
      httr2::response(status_code = 200L, headers = list(`Content-Type` = "text/xml"), body = xml_body)
    },
    {
      ef1 <- pubmed_efetch(c("111", "222"), email = "a@b.c", api_key = "k1", cache_dir = tmp_cache, use_cache = TRUE)
      expect_equal(nrow(ef1), 2L)
      expect_equal(ef1$pmid, c("111", "222"))
      expect_true("cfr" %in% ef1$pm_keywords[[1L]])
      # Exact batch cache hit
      ef_batch <- pubmed_efetch(c("111", "222"), cache_dir = tmp_cache, use_cache = TRUE)
      expect_equal(nrow(ef_batch), 2L)
      # Per-PMID cache hit for subset ("111")
      ef_sub <- pubmed_efetch("111", cache_dir = tmp_cache, use_cache = TRUE)
      expect_equal(nrow(ef_sub), 1L)
      # Delete per-PMID cache for "222" and verify harvest from batch cache file
      unlink(file.path(tmp_cache, "pm_art_v1_222.rds"))
      ef_harv <- pubmed_efetch("222", cache_dir = tmp_cache, use_cache = TRUE)
      expect_equal(nrow(ef_harv), 1L)
      expect_equal(ef_harv$pmid, "222")
    }
  )

  # Empty XML and failed efetch fallback
  httr2::with_mocked_responses(
    function(req) httr2::response(status_code = 200L, body = charToRaw("<PubmedArticleSet></PubmedArticleSet>")),
    {
      ef_empty <- pubmed_efetch("999", use_cache = FALSE)
      expect_equal(nrow(ef_empty), 0L)
    }
  )
  httr2::with_mocked_responses(
    function(req) stop("Network failure"),
    {
      expect_warning(ef_fail <- pubmed_efetch("998", use_cache = FALSE), "PubMed efetch failed")
      expect_equal(nrow(ef_fail), 0L)
    }
  )

  # 3. pubmed_eval_query: fast path (all IDs in first call) and multi-request batch UID path (total_hits > length(idlist))
  ev_fast <- pubmed_eval_query("h5n1[tiab]", seed_pmids = c("111", "222"), cache_dir = tmp_cache, use_cache = TRUE)
  expect_true(ev_fast$ok)
  expect_equal(ev_fast$pm_hits, 2L)
  expect_equal(ev_fast$pm_recall, 1.0)

  call_no <- 0L
  httr2::with_mocked_responses(
    function(req) {
      call_no <<- call_no + 1L
      if (call_no == 1L) {
        httr2::response(
          status_code = 200L,
          headers = list(`Content-Type` = "application/json"),
          body = charToRaw('{"esearchresult":{"count":"5000","idlist":["111"]}}')
        )
      } else {
        httr2::response(
          status_code = 200L,
          headers = list(`Content-Type` = "application/json"),
          body = charToRaw('{"esearchresult":{"count":"2","idlist":["111","222"]}}')
        )
      }
    },
    {
      ev_multi <- pubmed_eval_query("broad_q[tiab]", seed_pmids = c("111", "222"), use_cache = FALSE)
      expect_true(ev_multi$ok)
      expect_equal(ev_multi$pm_hits, 5000L)
      expect_equal(sort(ev_multi$pm_found), c("111", "222"))
      expect_equal(ev_multi$pm_recall, 1.0)
    }
  )
})

test_that("enrich_seeds resolves PMIDs via DOI, exact title, year-constrained title words, and fallback title words", {
  expect_equal(title_jaccard("Avian influenza H5N1", "Avian influenza H5N1"), 1.0)
  expect_equal(title_jaccard("", "Avian influenza"), 0)
  expect_equal(clean_title_for_query('Title with "quotes" and [brackets] (parens)'), "Title with quotes and brackets parens")
  expect_match(title_sig_words("The case fatality of avian influenza in the human population", 4L), "case fatality avian influenza")

  raw_seeds <- tibble::tibble(
    seed_id = 1:5,
    citekey = paste0("s", 1:5),
    entry_type = "article",
    title = c(
      "Already has PMID",
      "Resolved by DOI match",
      "Resolved by exact title query",
      "Resolved by significant words and year fallback",
      "Unresolved gray literature report"
    ),
    abstract = c("Short abs", "", "Existing longer abstract text here", "", ""),
    doi = c("", "10.1000/doi2", "", "", ""),
    pmid = c("101", "", "", "", ""),
    keywords = c("h5n1; cfr", "", "", "", ""),
    year = c("2008", "2010", "2012", "2015", ""),
    author = "Author",
    journal = "Journal"
  )

  mock_es <- function(term, retmax = 1L, ...) {
    if (grepl("10.1000/doi2", term, fixed = TRUE)) {
      return(list(count = 1L, idlist = "102", warnings = character(0), ok = TRUE))
    }
    if (grepl('"Resolved by exact title query"[title]', term, fixed = TRUE)) {
      return(list(count = 1L, idlist = "103", warnings = character(0), ok = TRUE))
    }
    if (grepl("2015[dp]", term, fixed = TRUE)) {
      return(list(count = 0L, idlist = character(0), warnings = character(0), ok = TRUE))
    }
    if (grepl("significant[title]", term, fixed = TRUE)) {
      return(list(count = 1L, idlist = "104", warnings = character(0), ok = TRUE))
    }
    list(count = 0L, idlist = character(0), warnings = character(0), ok = TRUE)
  }

  mock_ef <- function(pmids, ...) {
    tibble::tibble(
      pmid = c("101", "102", "103", "104"),
      pm_title = c(
        "Already has PMID",
        "Resolved by DOI match",
        "Resolved by exact title query",
        "Resolved by significant words and year fallback"
      ),
      pm_abstract = c(
        "Longer PubMed abstract replacing short BibTeX abstract",
        "Fetched abstract for DOI paper",
        "Short",
        "Fetched abstract for sig words paper"
      ),
      pm_mesh = list(c("influenza, human"), character(0), character(0), character(0)),
      pm_keywords = list(c("avian flu"), c("h5n1"), character(0), character(0))
    ) |>
      dplyr::filter(.data$pmid %in% pmids)
  }

  msgs <- character(0)
  enriched <- enrich_seeds(
    seeds_df = raw_seeds,
    pm_esearch_fn = mock_es,
    pm_efetch_fn = mock_ef,
    progress_cb = function(m) msgs <<- c(msgs, m)
  )
  expect_equal(enriched$pmid, c("101", "102", "103", "104", ""))
  expect_equal(enriched$in_pubmed, c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_true(length(msgs) == 5L)
  expect_match(enriched$full_abstract[1L], "Longer PubMed abstract")
  expect_match(enriched$full_abstract[3L], "Existing longer abstract")
})

test_that("Prefix wildcard synthesis, live background probing, igraph other-clustering, and Stage 2a/2b edge cases work", {
  old_opt <- options(revlitquery.min_interval = 0)
  on.exit(options(old_opt), add = TRUE)
  # 8 seeds with morphological variants (fatal/fatality/fatalities, mortality/mortalities, human/humans)
  # plus 2 non-PICO clusterable terms ("oseltamivir", "zanamivir")
  seeds8 <- tibble::tibble(
    seed_id = 1:8,
    citekey = paste0("p", 1:8),
    year = as.character(2008:2015),
    doi = "",
    pmid = as.character(201:208),
    in_pubmed = rep(TRUE, 8),
    full_title = c(
      "H5N1 avian influenza fatality in humans",
      "H5N1 avian influenza fatalities in human patients",
      "H5N1 avian influenza fatal outcomes in humans",
      "H5N1 avian influenza mortality in human cases",
      "H5N1 avian influenza mortalities and oseltamivir treatment",
      "H5N1 avian influenza fatal disease and oseltamivir therapy",
      "H5 avian influenza fatality and zanamivir oseltamivir",
      "H5 avian influenza zoonotic transmission in humans"
    ),
    full_abstract = full_title,
    full_text = full_title,
    mesh_terms = rep(list(character(0)), 8),
    all_keywords = rep(list(character(0)), 8)
  )

  # Test live probe callback & NA fallback when eval_bg_fn = NULL inside httr2::with_mocked_responses
  httr2::with_mocked_responses(
    function(req) {
      httr2::response(
        status_code = 200L,
        headers = list(`Content-Type` = "application/json"),
        body = charToRaw('{"esearchresult":{"count":"1200","idlist":[]}}')
      )
    },
    {
      p_msgs <- character(0)
      t_live <- extract_candidate_terms(
        seeds_df = seeds8,
        config = list(),
        max_candidates = 25L,
        probe_pubmed = TRUE,
        eval_bg_fn = NULL,
        use_cache = FALSE,
        progress_cb = function(m) p_msgs <<- c(p_msgs, m)
      )
      expect_true(nrow(t_live) > 0L)
      expect_true(length(p_msgs) > 0L)
      expect_true("fatal*" %in% t_live$term)

      # max_blocks = 1L and 3L branches (exercises Block 1 code-family synonym expansion h5n1 + h5)
      b_one <- build_concept_blocks(t_live, seeds8, max_blocks = 1L, probe_selectivity = TRUE)
      expect_equal(length(b_one), 1L)
      t_syn <- extract_candidate_terms(
        seeds_df = seeds8,
        probe_pubmed = TRUE,
        eval_bg_fn = function(t) if (t == "h5n1") 8000L else if (t == "h5") 15000L else 250000L
      )
      b_three <- build_concept_blocks(t_syn, seeds8, max_blocks = 3L, probe_selectivity = TRUE, eval_joint_fn = function(t1, t2) 120L)
      expect_true(length(b_three) >= 2L)
      expect_true(all(c("h5n1", "h5") %in% b_three[[1L]]$terms))
    }
  )

  # Test igraph fallback clustering on "other" terms when no PICO outcome/population/epi terms exist
  seeds_other <- tibble::tibble(
    seed_id = 1:6,
    citekey = paste0("o", 1:6),
    year = "2020",
    pmid = as.character(301:306),
    in_pubmed = rep(TRUE, 6),
    full_title = c(
      "H5N1 oseltamivir bindpocket",
      "H5N1 oseltamivir bindpocket structure",
      "H5N1 oseltamivir crystal",
      "H5N1 zanamivir neuraminidase",
      "H5N1 zanamivir neuraminidase inhibitor",
      "H5N1 zanamivir cleft"
    ),
    full_abstract = full_title,
    full_text = full_title,
    mesh_terms = rep(list(character(0)), 6),
    all_keywords = rep(list(character(0)), 6)
  )
  t_other <- extract_candidate_terms(
    seeds_df = seeds_other,
    probe_pubmed = TRUE,
    eval_bg_fn = function(t) if (t == "h5n1") 8000L else 400L
  )
  b_other <- build_concept_blocks(
    terms_df = t_other,
    seeds_df = seeds_other,
    max_blocks = 3L,
    probe_selectivity = TRUE,
    eval_joint_fn = function(t1, t2) 50L
  )
  expect_true(length(b_other) >= 2L)
  expect_match(b_other[[2L]]$name, "^concept_")

  # Test estimate_local_query edge cases, 1-block generate_initial_candidates, and Stage 2a/2b optimizer with exclusions & custom_queries
  expect_equal(estimate_local_query(make_query(list()), t_other, seeds_other)$local_found, 0L)
  q_unk <- make_query(
    list(
      list(name = "primary_topic", field = "tiab", terms = c("h5", "unknown phrase")),
      list(name = "domain_context", field = "tiab", terms = "influenza"),
      list(name = "concept_host", field = "title", terms = "unknown_single")
    ),
    exclude_terms = c("h7n9", "title:ferret")
  )
  loc_unk <- estimate_local_query(q_unk, t_other, seeds_other)
  expect_equal(length(loc_unk$local_cov), 6L)

  expect_equal(length(generate_initial_candidates(list(), t_other, seeds_other)), 0L)
  cands_1b <- generate_initial_candidates(b_other[1L], t_other, seeds_other, config = list(required_terms = "nonexistent_term_xyz"))
  expect_equal(length(cands_1b), 1L)
  expect_equal(length(compute_pareto_frontier(tibble::tibble())), 0L)
})

test_that("Full optimize_queries Stage 2a exclusion refinement, custom_queries config, and reporting blocks work", {
  seeds <- tibble::tibble(
    seed_id = 1:6,
    citekey = paste0("s", 1:6),
    year = as.character(2010:2015),
    doi = "",
    pmid = as.character(1:6),
    in_pubmed = rep(TRUE, 6),
    full_title = c(
      "H5N1 influenza fatality in humans",
      "H5N1 influenza fatal outcome in humans",
      "H5N1 influenza fatalities in human cases",
      "H5N1 influenza fatal pneumonia in humans",
      "H5N1 influenza deadly infection in humans",
      "H5N1 influenza zoonotic infection in humans"
    ),
    full_abstract = full_title,
    full_text = full_title,
    mesh_terms = rep(list(character(0)), 6),
    all_keywords = rep(list(character(0)), 6)
  )

  terms_df <- extract_candidate_terms(
    seeds_df = seeds,
    probe_pubmed = TRUE,
    eval_bg_fn = function(t) {
      switch(t, "h5n1" = 8000L, "h5" = 15000L, "influenza" = 120000L, "fatal*" = 150000L,
             "deadly" = 30000L, "zoonotic" = 40000L, "humans" = 900000L, "human" = 1200000L, 25000L)
    }
  )
  blocks <- list(
    list(name = "primary_topic", terms = c("h5n1")),
    list(name = "concept_outcome", terms = c("fatal*"))
  )
  attr(blocks, "domain_qualifier") <- "influenza"

  mock_eval_pm <- function(pm_query, seed_pmids, ...) {
    has_not <- grepl("NOT", pm_query, fixed = TRUE)
    has_zoon <- grepl("zoonotic", pm_query, fixed = TRUE)
    has_deadly <- grepl("deadly", pm_query, fixed = TRUE)
    found <- as.character(1:4)
    if (has_deadly) found <- union(found, "5")
    if (has_zoon) found <- union(found, "6")
    hits <- if (has_not) 180L else 650L
    list(pm_hits = hits, pm_found = found, pm_recall = length(found) / 6, ok = TRUE)
  }
  mock_es <- function(term, retmax = 0L, ...) {
    if (retmax > 0L) {
      return(list(count = 650L, idlist = c(as.character(1:6), paste0("n", 1:15)), ok = TRUE))
    }
    list(count = 140L, idlist = character(0), ok = TRUE)
  }
  mock_ef <- function(pmids, ...) {
    tibble::tibble(
      pmid = pmids,
      pm_title = rep("Recombinant H7N9 and H1N1 vaccine mice pathogenicity in cells", length(pmids)),
      pm_abstract = rep("Experimental vaccine immunization in mice and ferrets against H7N9 and H1N1(pdm09) virus.", length(pmids)),
      pm_mesh = rep(list(character(0)), length(pmids)),
      pm_keywords = rep(list(c("vaccine", "mice")), length(pmids))
    )
  }

  cfg <- list(
    custom_queries = list(
      "h5n1[tiab] AND fatal*[tiab]",
      list(label = "Named custom", pubmed_query = "h5n1[tiab] AND influenza[tiab]")
    )
  )

  res <- optimize_queries(
    seeds_df = seeds,
    terms_df = terms_df,
    blocks = blocks,
    config = cfg,
    min_recall = 0.80,
    max_hits = 500,
    max_evals = 32L,
    discover_exclusions = TRUE,
    eval_pubmed_fn = mock_eval_pm,
    pm_esearch_fn = mock_es,
    pm_efetch_fn = mock_ef,
    progress_cb = function(m) NULL
  )

  expect_true(nrow(res$exclusion_hints) > 0L)
  expect_true(any(res$eval_df$stage == "stage2_excl"))
  expect_equal(sum(res$eval_df$is_custom), 2L)

  # Re-adding duplicate custom query replaces existing row cleanly
  dup_c <- evaluate_custom_query("h5n1[tiab] AND fatal*[tiab]", seeds, label = "Updated lbl", gs_query = "h5n1 fatal", eval_pubmed_fn = mock_eval_pm)
  res_dup <- add_custom_query_to_result(res, dup_c)
  expect_equal(sum(res_dup$eval_df$is_custom), 2L)

  # Test optimize_queries when all seeds have in_pubmed = FALSE and max_hits is tight (triggers warning & fallback status)
  seeds_no_pm <- seeds
  seeds_no_pm$in_pubmed <- rep(FALSE, 6)
  seeds_no_pm$pmid <- rep("", 6)
  expect_warning(
    res_nopm <- optimize_queries(
      seeds_df = seeds_no_pm,
      terms_df = terms_df,
      blocks = blocks,
      config = list(custom_queries = "h5n1[tiab]"),
      min_recall = 0.99,
      max_hits = 10,
      max_evals = 6L,
      discover_exclusions = FALSE,
      eval_pubmed_fn = mock_eval_pm
    ),
    "satisfied max_hits"
  )
  expect_equal(res_nopm$recommended_status, "exceeds_max_hits")

  out_dir <- tempfile("rep_full_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  files <- write_optimization_report(res, output_dir = out_dir)
  md_txt <- paste(readLines(files$report_md), collapse = "\n")
  expect_match(md_txt, "Exclusion Concept Blocks", fixed = TRUE)
  expect_match(md_txt, "exclude_terms:", fixed = TRUE)
})

test_that("Unsupervised Hybrid Contrastive Graph + Set-Cover works across any scientific domain (Oncology, Ecology, Machine Learning)", {
  # 1. Test unsupervised exclusion co-occurrence clustering (cluster_exclusion_terms)
  expect_equal(cluster_exclusion_terms(character(0), list()), character(0))
  expect_equal(cluster_exclusion_terms("xenograft*", list(c(TRUE, FALSE))), "topic_xenograft")

  # Simulate 12 noise papers splitting into 2 distinct co-occurring noise subtopics:
  # Docs 1..6: preclinical mouse xenograft / cell culture noise
  # Docs 7..12: surgical resection / laparoscopy noise
  cov_noise <- list(
    c(rep(TRUE, 6), rep(FALSE, 6)),   # xenograft
    c(rep(TRUE, 5), FALSE, rep(FALSE, 6)), # mice
    c(FALSE, rep(TRUE, 5), rep(FALSE, 6)), # vitro
    c(rep(FALSE, 6), rep(TRUE, 6)),   # resection
    c(rep(FALSE, 6), rep(TRUE, 5), FALSE), # laparoscopic
    c(rep(FALSE, 6), FALSE, rep(TRUE, 5))  # biopsy
  )
  ex_terms <- c("xenograft", "mice", "vitro", "resection", "laparoscopic", "biopsy")
  ex_mesh <- c(
    rep(list(c("mice", "tumor cells, cultured")), 6),
    rep(list(c("surgical procedures, operative", "biopsy")), 6)
  )
  ex_lbls <- cluster_exclusion_terms(
    terms = ex_terms,
    noise_cov_list = cov_noise,
    rank_scores = c(100, 90, 80, 95, 85, 75),
    noise_mesh = ex_mesh
  )
  expect_equal(length(unique(ex_lbls)), 2L)
  expect_equal(ex_lbls[1L], ex_lbls[2L])
  expect_equal(ex_lbls[4L], ex_lbls[5L])
  expect_false(ex_lbls[1L] == ex_lbls[4L])

  # 2. Oncology corpus: unsupervised concept discovery + prefix wildcard synthesis
  onc_seeds <- tibble::tibble(
    seed_id = 1:6,
    citekey = paste0("onc", 1:6),
    year = as.character(2019:2024),
    doi = "",
    pmid = as.character(501:506),
    in_pubmed = rep(TRUE, 6),
    full_title = c(
      "Pembrolizumab immunotherapy in melanoma metastasis randomized trial",
      "Pembrolizumab therapy in melanoma metastases randomized clinical trial",
      "Pembrolizumab treatment in metastatic melanoma patients remission",
      "Nivolumab immunotherapy in melanoma metastasis progression survival",
      "Nivolumab treatment in metastatic melanoma randomized cohort",
      "Checkpoint inhibitor immunotherapy in melanoma metastases remission"
    ),
    full_abstract = full_title,
    full_text = full_title,
    mesh_terms = rep(list(c("melanoma", "immunotherapy")), 6),
    all_keywords = rep(list(c("checkpoint inhibitor")), 6)
  )

  onc_terms <- extract_candidate_terms(
    seeds_df = onc_seeds,
    probe_pubmed = TRUE,
    eval_bg_fn = function(t) if (t == "melanoma") 140000L else 40000L
  )
  expect_true("metasta*" %in% onc_terms$term)
  onc_blocks <- build_concept_blocks(
    terms_df = onc_terms,
    seeds_df = onc_seeds,
    max_blocks = 3L,
    probe_selectivity = TRUE,
    eval_joint_fn = function(t1, t2) 2500L
  )
  expect_true(length(onc_blocks) >= 2L)
  q_onc <- make_query(list(
    list(name = "topic", field = "tiab", terms = c("melanoma", "metasta*")),
    list(name = "intervention", field = "tiab", terms = "immunotherap*")
  ))
  expect_match(to_google_scholar_query(q_onc), "metastases|metastasis|metastatic")
  expect_match(to_google_scholar_query(q_onc), "immunotherapy", fixed = TRUE)

  # 3. Non-medical Ecology / Biodiversity corpus (completely outside biomedicine/PICO)
  eco_seeds <- tibble::tibble(
    seed_id = 1:6,
    citekey = paste0("eco", 1:6),
    year = as.character(2018:2023),
    doi = "",
    pmid = as.character(601:606),
    in_pubmed = rep(TRUE, 6),
    full_title = c(
      "Pollinator biodiversity decline in fragmented agricultural landscapes",
      "Wildbee pollinator richness decline across fragmented farmland landscapes",
      "Bumblebee pollinator abundance loss in fragmented agricultural habitats",
      "Pollinator diversity loss in fragmented agroecosystem landscapes",
      "Solitary bee pollinator decline across fragmented cropland habitats",
      "Pollinator richness loss in fragmented agricultural ecosystems"
    ),
    full_abstract = full_title,
    full_text = full_title,
    mesh_terms = rep(list(character(0)), 6),
    all_keywords = rep(list(character(0)), 6)
  )
  eco_terms <- extract_candidate_terms(
    seeds_df = eco_seeds,
    probe_pubmed = FALSE
  )
  eco_blocks <- build_concept_blocks(
    terms_df = eco_terms,
    seeds_df = eco_seeds,
    max_blocks = 3L,
    probe_selectivity = FALSE
  )
  expect_true(length(eco_blocks) >= 2L)
  all_eco_terms <- unlist(lapply(eco_blocks, `[[`, "terms"))
  expect_true(any(grepl("pollinat|fragment", all_eco_terms)))
})

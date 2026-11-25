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

  # Semantic family coherence checks
  expect_equal(assign_term_family("case fatality"), "outcome_severity")
  expect_equal(assign_term_family("mortality"), "outcome_severity")
  expect_equal(assign_term_family("died"), "outcome_severity")
  expect_equal(assign_term_family("human cases"), "population_host")
  expect_equal(assign_term_family("patients"), "population_host")
  expect_equal(assign_term_family("seroprevalence"), "epi_methods")
  expect_equal(assign_term_family("respiratory"), "other")

  blk_names <- vapply(blocks, `[[`, character(1), "name")
  expect_equal(blk_names[1L], "primary_topic")
  expect_true("outcome_severity" %in% blk_names)
  out_blk <- blocks[[which(blk_names == "outcome_severity")[1L]]]
  expect_true(all(vapply(out_blk$terms, assign_term_family, character(1)) == "outcome_severity"))
  if ("population_host" %in% blk_names) {
    pop_blk <- blocks[[which(blk_names == "population_host")[1L]]]
    expect_true(all(vapply(pop_blk$terms, assign_term_family, character(1)) == "population_host"))
  }

  # Stage 2b rescue keeps outcome_severity pure (does not inject population_host terms)
  q_out_only <- make_query(list(
    list(name = "primary_topic", field = "tiab", terms = "h5n1"),
    list(name = "outcome_severity", field = "tiab", terms = "case fatality")
  ))
  nbrs <- .generate_neighbors(
    pareto_queries = list(q_out_only),
    pareto_missed_idx = list(c(3L, 4L)),
    blocks = blocks,
    terms_df = terms_df,
    seeds_df = seeds,
    config = list(),
    seen_sigs = query_signature(q_out_only)
  )
  for (nq in nbrs) {
    for (nb_blk in nq$blocks) {
      if (nb_blk$name == "outcome_severity") {
        expect_true(all(vapply(nb_blk$terms, assign_term_family, character(1)) == "outcome_severity"))
      }
    }
  }

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
  expect_equal(assign_exclusion_concept("h7n9"), "competing_subtypes")
  expect_equal(assign_exclusion_concept("vaccin*"), "vaccines_immunology")
  expect_equal(assign_exclusion_concept("mice"), "experimental_lab")
  expect_equal(assign_exclusion_concept("poultry"), "veterinary_ecology")
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
      list(name = "outcome_severity", field = "tiab", terms = c("fatal*", "fatality", "case fatality"))
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

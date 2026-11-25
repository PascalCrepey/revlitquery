# R/09_exclusion_hints.R
# Mine non-seed PubMed hits ("noise" articles) to discover concept-organized,
# wildcard-enabled exclusion terms (NOT / -) and greedy set-cover exclusion
# blocks that reduce hit count while preserving seed recall.

.empty_exclusion_hints <- function() {
  tbl <- tibble::tibble(
    concept = character(0),
    term = character(0),
    field = character(0),
    exclude_spec = character(0),
    pm_clause = character(0),
    gs_clause = character(0),
    noise_pct = numeric(0),
    exact_hits_removed = integer(0),
    pct_hits_reduced = numeric(0),
    seeds_lost_n = integer(0),
    seeds_lost_citekeys = character(0),
    safety_tier = character(0)
  )
  attr(tbl, "exclusion_blocks") <- list()
  attr(tbl, "setcover_bundles") <- list()
  attr(tbl, "code_bundle") <- character(0)
  attr(tbl, "all_code_ex") <- character(0)
  tbl
}

#' Classify a candidate exclusion term into a scientific exclusion concept
#'
#' @param term Single candidate exclusion term (may end with `*`).
#' @return One of `"competing_subtypes"`, `"experimental_lab"`,
#'   `"vaccines_immunology"`, `"veterinary_ecology"`, or `"general_noise"`.
assign_exclusion_concept <- function(term) {
  t_raw <- tolower(stringr::str_squish(as.character(term %||% "")))
  if (!nzchar(t_raw)) return("general_noise")
  t_low <- gsub("\\*", "", t_raw)

  # 1. Competing virus subtypes, strains, or other pathogens
  if (grepl("^[a-z]+[0-9]+[a-z0-9]*$", t_low) ||
      grepl("\\b(h[0-9]+n[0-9]+|pdm09|sars[a-z0-9]*|covid[a-z0-9]*|mers[a-z0-9]*|ndv|newcastle)\\b", t_low, perl = TRUE)) {
    return("competing_subtypes")
  }
  # 2. Vaccines, immunology & prophylaxis trials
  if (grepl("\\b(vaccin[a-z]*|immun[a-z]*|antigen[a-z]*|antibod[a-z]*|adjuvant[a-z]*|neutraliz[a-z]*|prophylax[a-z]*|titer[a-z]*|titres|dose[a-z]*|challenge[a-z]*|protect[a-z]*|efficacy)\\b", t_low, perl = TRUE)) {
    return("vaccines_immunology")
  }
  # 3. Cellular, molecular & experimental animal lab models
  if (grepl("\\b(cell[a-z]*|mice|mouse|ferret[a-z]*|guinea|vitro|vivo|receptor[a-z]*|binding|gene|genes|express[a-z]*|cytokin[a-z]*|epithel[a-z]*|polymerase|pb2|pb1|haemagglutinin|hemagglutinin|neuraminidase|nucle[a-z]*|necro[a-z]*|interfer[a-z]*|mutat[a-z]*|amino|protein[a-z]*|replicat[a-z]*|assay[a-z]*|plasmid[a-z]*|recombin[a-z]*|crystal[a-z]*|structur[a-z]*|strain[a-z]*|isolat[a-z]*|pathogenicity|lethal[a-z]*|lungs|lung|brain|organs|tissue[a-z]*|systemic)\\b", t_low, perl = TRUE)) {
    return("experimental_lab")
  }
  # 4. Veterinary, wildlife, poultry, livestock & phylogenetic ecology
  if (grepl("\\b(duck[a-z]*|mallard[a-z]*|chicken[a-z]*|poultry|bird[a-z]*|waterfowl|migratory|swan[a-z]*|geese|goose|flock[a-z]*|farm[a-z]*|cattle|bovin[a-z]*|dairy|milk|swine|pig|pigs|cat|cats|seal[a-z]*|mink[a-z]*|aiv|aivs|hpaiv[a-z]*|hpai|lpaiv[a-z]*|lpai|clade[a-z]*|lineage[a-z]*|phylo[a-z]*|phyloge[a-z]*|reassort[a-z]*|ecology|ecological|reservoir[a-z]*|wild)\\b", t_low, perl = TRUE)) {
    return("veterinary_ecology")
  }
  "general_noise"
}

# Synthesize prefix wildcards (e.g. vaccin*, cell*, mutat*, phylo*, reassort*, receptor*, duck*, ferret*)
# from noise unigrams when >= 2 noise words share a 5..8-char prefix and have 0 seed loss in tiab or title
.synthesize_exclusion_wildcards <- function(doc_grams, noise_texts, noise_titles,
                                            seed_texts, seed_titles, min_noise_cnt, stops) {
  all_u <- unique(unlist(doc_grams))
  u <- all_u[grepl("^[a-z]{4,}$", all_u) & !(all_u %in% stops)]
  if (length(u) < 2L) return(character(0))
  # Restrict prefix wildcard synthesis to scientific exclusion concepts (avoid wildcarding generic English verbs/adverbs)
  u_cpts <- vapply(u, assign_exclusion_concept, character(1))
  u <- u[u_cpts != "general_noise"]
  if (length(u) < 2L) return(character(0))

  u_counts <- table(unlist(lapply(doc_grams, function(g) intersect(g, u))))
  cand_prefixes <- unique(unlist(lapply(5:8, function(len) {
    long_u <- u[nchar(u) >= len]
    if (length(long_u) == 0L) character(0) else substr(long_u, 1L, len)
  })))
  cand_prefixes <- cand_prefixes[!(cand_prefixes %in% stops)]
  if (length(cand_prefixes) == 0L) return(character(0))

  matched_sets <- lapply(cand_prefixes, function(p) sort(u[startsWith(u, p)]))
  keep_p <- vapply(matched_sets, length, integer(1)) >= 2L
  cand_prefixes <- cand_prefixes[keep_p]
  matched_sets <- matched_sets[keep_p]
  if (length(cand_prefixes) == 0L) return(character(0))

  set_sigs <- vapply(matched_sets, paste, character(1), collapse = "|")
  ord_len <- order(nchar(cand_prefixes), decreasing = TRUE)
  cand_prefixes <- cand_prefixes[ord_len]
  matched_sets <- matched_sets[ord_len]
  set_sigs <- set_sigs[ord_len]
  uniq_mask <- !duplicated(set_sigs)
  cand_prefixes <- cand_prefixes[uniq_mask]
  matched_sets <- matched_sets[uniq_mask]

  kept <- character(0)
  kept_roots <- character(0)
  kept_gains <- integer(0)

  for (i in seq_along(cand_prefixes)) {
    p <- cand_prefixes[i]
    m_u <- matched_sets[[i]]
    max_single <- max(as.integer(u_counts[m_u]), na.rm = TRUE)
    w_term <- paste0(p, "*")

    # Check if w_term is safe (0 seed loss) in tiab or title and removes noise
    loss_tiab <- sum(term_matches_text(w_term, seed_texts))
    loss_ti <- sum(term_matches_text(w_term, seed_titles))
    if (loss_tiab > 0L && loss_ti > 0L) next

    hits_tiab <- sum(term_matches_text(w_term, noise_texts))
    hits_ti <- sum(term_matches_text(w_term, noise_titles))
    eff_hits <- if (loss_tiab == 0L) hits_tiab else hits_ti

    # Keep if it removes >= min_noise_cnt noise papers and unifies >= 2 forms cleanly
    if ((loss_tiab == 0L && hits_tiab >= min_noise_cnt && (hits_tiab > max_single || length(m_u) >= 3L)) ||
        (loss_tiab > 0L && loss_ti == 0L && hits_ti >= max(3L, floor(min_noise_cnt * 0.75)))) {
      root5 <- substr(p, 1L, 5L)
      ex_idx <- match(root5, kept_roots)
      if (is.na(ex_idx)) {
        kept <- c(kept, w_term)
        kept_roots <- c(kept_roots, root5)
        kept_gains <- c(kept_gains, eff_hits)
      } else if (eff_hits > kept_gains[ex_idx]) {
        kept[ex_idx] <- w_term
        kept_gains[ex_idx] <- eff_hits
      }
    }
  }
  kept
}

#' Discover concept-grouped candidate exclusion (NOT / -) terms from non-seed PubMed results
#'
#' Samples articles returned by `base_query` that are NOT in `seeds_df` across the full
#' publication timeline, mines frequent 1- and 2-grams plus morphological prefix wildcards
#' (`*`), checks seed loss locally against PubMed-resolved seeds (including PubMed Author
#' Keywords), groups exclusions by concept, builds greedy set-cover exclusion blocks,
#' and returns up to `max_hints` ranked exclusion hints.
#'
#' @param base_query Canonical `search_query` object to analyze.
#' @param seeds_df Enriched seeds tibble.
#' @param base_pm_hits Integer total PubMed hits of `base_query`.
#' @param sample_n Maximum number of non-seed PubMed articles to sample.
#' @param top_verify Number of top exclusion candidates to verify live via `esearch`.
#' @param max_hints Maximum number of exclusion hints to return and display (default 35).
#' @param email Optional email for NCBI polite pool.
#' @param ncbi_key Optional NCBI API key.
#' @param cache_dir Cache directory path.
#' @param use_cache Logical; whether to read/write cache.
#' @param pm_esearch_fn Optional dependency injection for [pubmed_esearch()].
#' @param pm_efetch_fn Optional dependency injection for [pubmed_efetch()].
#' @param progress_cb Optional progress callback function.
#' @return Tibble of ranked exclusion hints with `concept` column and `exclusion_blocks` attribute.
discover_exclusion_hints <- function(base_query,
                                     seeds_df,
                                     base_pm_hits = NA_integer_,
                                     sample_n = 160L,
                                     top_verify = 12L,
                                     max_hints = 35L,
                                     email = NULL,
                                     ncbi_key = NULL,
                                     cache_dir = ".cache",
                                     use_cache = TRUE,
                                     pm_esearch_fn = pubmed_esearch,
                                     pm_efetch_fn = pubmed_efetch,
                                     progress_cb = NULL) {
  if (is.null(base_query) || length(base_query$blocks) == 0L) return(.empty_exclusion_hints())

  pm_base <- to_pubmed_query(base_query)
  if (!nzchar(pm_base)) return(.empty_exclusion_hints())

  seed_pmids <- unique(as.character(seeds_df$pmid[nzchar(seeds_df$pmid)]))
  retmax <- min(max(as.integer(sample_n) * 4L + length(seed_pmids), 600L), 2500L)

  if (is.function(progress_cb)) progress_cb("Sampling non-seed PubMed results for exclusion hints...")
  res_ids <- pm_esearch_fn(
    term = pm_base, retmax = retmax,
    email = email, api_key = ncbi_key, cache_dir = cache_dir, use_cache = use_cache
  )
  if (!isTRUE(res_ids$ok) || length(res_ids$idlist) == 0L) return(.empty_exclusion_hints())

  if (is.na(base_pm_hits)) base_pm_hits <- as.integer(res_ids$count %||% length(res_ids$idlist))
  noise_all <- setdiff(res_ids$idlist, seed_pmids)
  if (length(noise_all) < 5L) return(.empty_exclusion_hints())

  # Stratified sample across the entire returned ID list to cover all publication years
  if (length(noise_all) > as.integer(sample_n)) {
    strat_idx <- unique(round(seq(1, length(noise_all), length.out = as.integer(sample_n))))
    noise_pmids <- noise_all[strat_idx]
  } else {
    noise_pmids <- noise_all
  }

  noise_meta <- pm_efetch_fn(
    pmids = noise_pmids, email = email, api_key = ncbi_key,
    cache_dir = cache_dir, use_cache = use_cache
  )
  if (is.null(noise_meta) || nrow(noise_meta) < 5L) return(.empty_exclusion_hints())

  n_noise <- nrow(noise_meta)
  noise_titles <- as.character(noise_meta$pm_title %||% "")
  noise_kw_str <- vapply(noise_meta$pm_keywords %||% rep(list(character(0)), n_noise), function(kw) {
    paste(kw, collapse = " ")
  }, character(1))
  noise_texts <- stringr::str_squish(paste(noise_titles, as.character(noise_meta$pm_abstract %||% ""), noise_kw_str))
  norm_noise <- normalize_search_text(noise_texts)

  stops <- unique(c(
    get_stopwords(),
    "demonstrate", "demonstrated", "ability", "various", "long", "efficiently",
    "signs", "shown", "additional", "similar", "different", "recent", "previously",
    "present", "suggest", "suggested", "indicate", "indicated", "provide", "provided", "provides",
    "subsequent", "subsequently", "exhibit", "exhibited", "exhibiting", "activate", "activated"
  ))
  pos_terms <- unique(tolower(gsub("\\*", "", unlist(lapply(base_query$blocks, `[[`, "terms")))))
  pos_stems <- unique(unlist(lapply(pos_terms, term_stems)))

  # Evaluate seed loss against PubMed-resolved seeds when available (so non-PubMed records don't block valid PubMed exclusions)
  eval_seeds_df <- if (any(seeds_df$in_pubmed %||% FALSE)) {
    seeds_df[which(seeds_df$in_pubmed), , drop = FALSE]
  } else {
    seeds_df
  }
  seed_titles <- as.character(eval_seeds_df$full_title %||% eval_seeds_df$title)
  seed_texts <- as.character(eval_seeds_df$full_text %||% paste(seed_titles, eval_seeds_df$abstract))
  n_seeds <- nrow(eval_seeds_df)
  max_allow_loss <- max(2L, floor(0.08 * n_seeds))

  doc_grams <- lapply(norm_noise, .extract_doc_ngrams, max_n = 2L, stops = stops)
  gram_df <- table(unlist(doc_grams))
  min_noise_cnt <- max(3L, ceiling(0.03 * n_noise))
  freq_cands <- names(gram_df[gram_df >= min_noise_cnt])

  # Include parenthesized compound subtype codes like "A(H1N1)pdm09" -> "h1n1pdm09"
  raw_code_2grams <- names(gram_df)[grepl("^[a-z]+[0-9]+[a-z0-9]* [a-z]+[0-9]+[a-z0-9]*$", names(gram_df))]
  raw_code_2grams <- Filter(function(g2) {
    parts <- strsplit(g2, " ", fixed = TRUE)[[1]]
    any(grepl(paste0("(?i)(?<![a-z0-9])", parts[1L], "[()]+", parts[2L], "(?![a-z0-9])"), noise_texts, perl = TRUE))
  }, raw_code_2grams)
  merged_code_2grams <- gsub(" ", "", raw_code_2grams, fixed = TRUE)
  code_tokens <- unique(c(
    names(gram_df)[grepl("^[a-z]+[0-9]+[a-z0-9]*$", names(gram_df))],
    merged_code_2grams
  ))

  # Synthesize prefix wildcards (e.g. vaccin*, cell*, mutat*, phylo*, reassort*, receptor*, duck*, ferret*)
  wildcard_cands <- .synthesize_exclusion_wildcards(
    doc_grams = doc_grams,
    noise_texts = noise_texts,
    noise_titles = noise_titles,
    seed_texts = seed_texts,
    seed_titles = seed_titles,
    min_noise_cnt = min_noise_cnt,
    stops = stops
  )

  cands <- unique(c(wildcard_cands, freq_cands, code_tokens))
  cands <- cands[nchar(cands) >= 3L & !grepl("^[0-9\\s]+$", cands)]

  # Drop terms that overlap positive query terms, their code family, or their stems
  cands <- Filter(function(t_str) {
    t_bare <- gsub("\\*", "", t_str)
    if (t_bare %in% pos_terms) return(FALSE)
    if (any(vapply(pos_terms, .is_code_family, logical(1), t2 = t_bare))) return(FALSE)
    if (any(vapply(pos_terms, function(pt) grepl(paste0("(?<![a-z0-9])", escape_regex(pt), "(?![a-z0-9])"), t_bare, perl = TRUE), logical(1)))) return(FALSE)
    if (length(intersect(term_stems(t_bare), pos_stems)) > 0L) return(FALSE)
    TRUE
  }, cands)
  if (length(cands) == 0L) return(.empty_exclusion_hints())

  in_gram <- intersect(cands, names(gram_df))
  cands_ord <- in_gram[order(gram_df[in_gram], decreasing = TRUE)]
  cands <- unique(c(wildcard_cands, utils::head(cands_ord, 240L), intersect(cands, code_tokens)))

  zero_loss_codes <- Filter(function(ct) {
    sum(term_matches_text(ct, seed_texts)) == 0L
  }, intersect(cands, code_tokens))

  scored <- list()
  for (t_str in cands) {
    t_bare <- gsub("\\*", "", t_str)
    is_code <- t_bare %in% code_tokens
    is_wild <- endsWith(t_str, "*")
    cov_n_tiab <- term_matches_text(t_str, noise_texts)
    cov_n_ti <- term_matches_text(t_str, noise_titles)
    n_cov_tiab <- sum(cov_n_tiab) / n_noise
    n_cov_ti <- sum(cov_n_ti) / n_noise

    s_loss_tiab_idx <- which(term_matches_text(t_str, seed_texts))
    s_loss_ti_idx <- which(term_matches_text(t_str, seed_titles))
    n_loss_tiab <- length(s_loss_tiab_idx)
    n_loss_ti <- length(s_loss_ti_idx)

    # Choose best scope: tiab if 0 seeds lost, else title if 0 seeds lost in title
    use_scope <- NULL
    min_tiab_rate <- if (is_code) 0.025 else 0.035
    if (n_loss_tiab == 0L && n_cov_tiab >= min_tiab_rate) {
      use_scope <- list(field = "tiab", noise_rate = n_cov_tiab, noise_cov = cov_n_tiab, loss_idx = s_loss_tiab_idx)
    } else if (n_loss_ti == 0L && n_cov_ti >= 0.03) {
      use_scope <- list(field = "title", noise_rate = n_cov_ti, noise_cov = cov_n_ti, loss_idx = s_loss_ti_idx)
    } else if (n_loss_tiab <= max_allow_loss && n_cov_tiab >= 0.06 && n_loss_tiab <= n_loss_ti) {
      use_scope <- list(field = "tiab", noise_rate = n_cov_tiab, noise_cov = cov_n_tiab, loss_idx = s_loss_tiab_idx)
    } else if (n_loss_ti <= max_allow_loss && n_cov_ti >= 0.045) {
      use_scope <- list(field = "title", noise_rate = n_cov_ti, noise_cov = cov_n_ti, loss_idx = s_loss_ti_idx)
    }

    if (is.null(use_scope)) next
    est_rem <- as.integer(round(base_pm_hits * use_scope$noise_rate))
    n_lost <- length(use_scope$loss_idx)
    concept <- assign_exclusion_concept(t_str)
    scope_bonus <- if (use_scope$field == "tiab" && n_lost == 0L) 1.35 else 1.0
    code_bonus <- if (is_code && use_scope$field == "tiab" && n_lost == 0L) 1.50 else 1.0
    wild_bonus <- if (is_wild && n_lost == 0L) 1.35 else 1.0
    unigram_bonus <- if (!grepl(" ", t_str, fixed = TRUE)) 1.15 else 1.0
    concept_bonus <- if (concept != "general_noise") 1.20 else 0.90

    scored[[length(scored) + 1L]] <- list(
      concept = concept,
      term = t_str,
      is_code = is_code,
      is_wild = is_wild,
      field = use_scope$field,
      exclude_spec = if (use_scope$field == "title") paste0("title:", t_str) else t_str,
      noise_pct = round(100 * use_scope$noise_rate, 1),
      noise_cov = use_scope$noise_cov,
      est_removed = est_rem,
      seeds_lost_n = n_lost,
      seeds_lost_citekeys = paste(eval_seeds_df$citekey[use_scope$loss_idx], collapse = "; "),
      rank_score = (est_rem * scope_bonus * code_bonus * wild_bonus * unigram_bonus * concept_bonus) / (1 + 6 * n_lost)
    )
  }

  if (length(scored) == 0L) return(.empty_exclusion_hints())

  # Deduplicate redundant literal words that are already subsumed by a safe wildcard in the same scope
  scored_ord <- scored[order(vapply(scored, `[[`, numeric(1), "rank_score"), decreasing = TRUE)]
  keep_sc <- rep(TRUE, length(scored_ord))
  for (i in seq_along(scored_ord)) {
    if (!keep_sc[i]) next
    it_i <- scored_ord[[i]]
    if (isTRUE(it_i$is_wild) && it_i$seeds_lost_n == 0L) {
      for (j in seq_along(scored_ord)) {
        if (i == j || !keep_sc[j]) next
        it_j <- scored_ord[[j]]
        if (identical(it_i$field, it_j$field) && it_j$seeds_lost_n >= it_i$seeds_lost_n &&
            .term_subsumes(it_i$term, it_j$term)) {
          keep_sc[j] <- FALSE
        }
      }
    }
  }
  scored_ord <- scored_ord[keep_sc]

  # Select diverse candidates across all exclusion concepts (up to max_hints, e.g. 35)
  code_idx <- which(vapply(scored_ord, function(x) isTRUE(x$is_code) && x$seeds_lost_n == 0L, logical(1)))
  wild_idx <- which(vapply(scored_ord, function(x) isTRUE(x$is_wild) && x$seeds_lost_n == 0L, logical(1)))
  concept_picks <- integer(0)
  for (cpt in c("competing_subtypes", "experimental_lab", "vaccines_immunology", "veterinary_ecology", "general_noise")) {
    c_rows <- which(vapply(scored_ord, `[[`, character(1), "concept") == cpt)
    if (length(c_rows) > 0L) concept_picks <- c(concept_picks, utils::head(c_rows, 8L))
  }
  all_picks <- sort(unique(c(
    seq_len(min(15L, length(scored_ord))),
    utils::head(code_idx, 5L),
    utils::head(wild_idx, 6L),
    concept_picks
  )))
  sel_cands <- scored_ord[utils::head(all_picks, as.integer(max_hints))]

  n_verify <- min(as.integer(top_verify), length(sel_cands))
  rows <- vector("list", length(sel_cands))
  for (i in seq_along(sel_cands)) {
    item <- sel_cands[[i]]
    tag <- if (item$field == "title") "[ti]" else "[tiab]"
    pm_cl <- sprintf('"%s"%s', item$term, tag)
    gs_cl <- paste0("-", .gs_format_term(item$term, in_title = identical(item$field, "title")))

    if (i <= n_verify) {
      if (is.function(progress_cb)) {
        progress_cb(sprintf("Verifying exclusion hint %d/%d: %s", i, n_verify, item$exclude_spec))
      }
      check_q <- sprintf("(%s) AND %s", pm_base, pm_cl)
      res_v <- pm_esearch_fn(
        term = check_q, retmax = 0L,
        email = email, api_key = ncbi_key, cache_dir = cache_dir, use_cache = use_cache
      )
      exact_rem <- if (isTRUE(res_v$ok) && !is.na(res_v$count)) as.integer(res_v$count) else item$est_removed
    } else {
      exact_rem <- as.integer(item$est_removed)
    }

    pct_red <- round(100 * exact_rem / max(as.numeric(base_pm_hits), 1), 1)
    tier <- if (item$seeds_lost_n == 0L) {
      "Safe (0 seeds lost)"
    } else if (item$seeds_lost_n <= 2L) {
      paste0("Low-risk (", item$seeds_lost_n, " lost)")
    } else {
      paste0("Trade-off (", item$seeds_lost_n, " lost)")
    }

    rows[[i]] <- tibble::tibble(
      concept = item$concept,
      term = item$term,
      field = item$field,
      exclude_spec = item$exclude_spec,
      pm_clause = pm_cl,
      gs_clause = gs_cl,
      noise_pct = item$noise_pct,
      exact_hits_removed = exact_rem,
      pct_hits_reduced = pct_red,
      seeds_lost_n = as.integer(item$seeds_lost_n),
      seeds_lost_citekeys = item$seeds_lost_citekeys,
      safety_tier = tier,
      .is_verified = (i <= n_verify)
    )
  }

  out_df <- dplyr::bind_rows(rows) |>
    dplyr::arrange(.data$seeds_lost_n, dplyr::desc(.data$.is_verified), dplyr::desc(.data$exact_hits_removed)) |>
    dplyr::select(-".is_verified")

  # Re-align noise_cov list with out_df order for greedy set-cover bundling
  spec_map <- stats::setNames(lapply(sel_cands, `[[`, "noise_cov"), vapply(sel_cands, `[[`, character(1), "exclude_spec"))
  safe_rows <- which(out_df$seeds_lost_n == 0L & out_df$exact_hits_removed > 0L)

  # 1. Build Concept-Grouped Exclusion Blocks (each block is a greedy set-cover of safe exclusions in that concept)
  exclusion_blocks <- list()
  for (cpt in c("competing_subtypes", "experimental_lab", "vaccines_immunology", "veterinary_ecology", "general_noise")) {
    c_idx <- safe_rows[out_df$concept[safe_rows] == cpt]
    if (length(c_idx) == 0L) next
    cov_cur <- rep(FALSE, n_noise)
    picked_specs <- character(0)
    while (length(picked_specs) < 5L) {
      rem <- setdiff(out_df$exclude_spec[c_idx], picked_specs)
      if (length(rem) == 0L) break
      gains <- vapply(rem, function(sp) {
        cv <- spec_map[[sp]]
        if (is.null(cv)) return(0L)
        sum(cv & !cov_cur)
      }, integer(1))
      if (max(gains) <= 0L) break
      best_sp <- rem[which.max(gains)]
      picked_specs <- c(picked_specs, best_sp)
      cov_cur <- cov_cur | spec_map[[best_sp]]
    }
    if (length(picked_specs) > 0L) {
      exclusion_blocks[[length(exclusion_blocks) + 1L]] <- list(
        name = cpt,
        terms = picked_specs,
        noise_covered_n = sum(cov_cur),
        noise_covered_pct = round(100 * sum(cov_cur) / n_noise, 1)
      )
    }
  }

  # 2. Build global greedy set-cover exclusion sequence across safe scientific concept exclusions
  global_setcover <- character(0)
  if (length(safe_rows) > 0L) {
    sci_safe_rows <- safe_rows[out_df$concept[safe_rows] != "general_noise"]
    use_safe_rows <- if (length(sci_safe_rows) >= 3L) sci_safe_rows else safe_rows
    cov_g <- rep(FALSE, n_noise)
    avail_sp <- out_df$exclude_spec[use_safe_rows]
    while (length(global_setcover) < min(14L, length(avail_sp))) {
      rem <- setdiff(avail_sp, global_setcover)
      if (length(rem) == 0L) break
      gains <- vapply(rem, function(sp) {
        cv <- spec_map[[sp]]
        if (is.null(cv)) return(-1)
        new_n <- sum(cv & !cov_g)
        if (new_n == 0L) return(-1)
        # Penalize long Google Scholar exclusion tokens slightly so set-cover stays compact
        gs_len <- nchar(out_df$gs_clause[match(sp, out_df$exclude_spec)])
        new_n / (gs_len ^ 0.35)
      }, numeric(1))
      if (max(gains) <= 0) break
      best_sp <- rem[which.max(gains)]
      global_setcover <- c(global_setcover, best_sp)
      cov_g <- cov_g | spec_map[[best_sp]]
    }
  }

  # 3. Build sibling subtype/code exclusion bundle (top 2 verified zero-loss tiab codes + any zero-loss prefix extensions)
  verified_safe_codes <- out_df$term[
    out_df$seeds_lost_n == 0L &
      out_df$field == "tiab" &
      grepl("^[a-z]+[0-9]+[a-z0-9]*$", out_df$term) &
      out_df$exact_hits_removed > 0L
  ]
  if (length(verified_safe_codes) > 0L) {
    top2_codes <- utils::head(verified_safe_codes, 2L)
    if (all(c("h7n9", "h1n1") %in% top2_codes)) {
      top2_codes <- c("h7n9", "h1n1")
    }
    ext_codes <- Filter(function(ct) {
      !(ct %in% top2_codes) && any(vapply(top2_codes, function(base_c) {
        nchar(ct) > nchar(base_c) && startsWith(ct, base_c)
      }, logical(1)))
    }, zero_loss_codes)
    attr(out_df, "code_bundle") <- unique(c(top2_codes, sort(ext_codes)))
    attr(out_df, "all_code_ex") <- unique(c(verified_safe_codes, sort(ext_codes)))
  } else {
    attr(out_df, "code_bundle") <- character(0)
    attr(out_df, "all_code_ex") <- character(0)
  }

  attr(out_df, "exclusion_blocks") <- exclusion_blocks
  attr(out_df, "setcover_ex") <- global_setcover
  out_df
}

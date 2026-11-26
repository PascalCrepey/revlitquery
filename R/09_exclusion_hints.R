# R/09_exclusion_hints.R
# Unsupervised exclusion concept discovery and wildcard mining on non-seed PubMed
# articles ("noise" sample): NPMI + Jaccard + MeSH/keyword co-occurrence graph
# community detection (`igraph::cluster_louvain()`) and greedy set-cover blocks.

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

#' Unsupervised co-occurrence graph clustering of candidate exclusion terms
#'
#' Partitions candidate exclusion terms into data-driven noise topic clusters
#' based on their empirical co-occurrence across sampled non-seed articles
#' (combining Normalized Pointwise Mutual Information, Jaccard overlap, shared
#' MeSH/Author Keyword profiles, stem similarity, and alphanumeric code pattern
#' affinity), using Louvain modularity maximization ([igraph::cluster_louvain()]).
#'
#' @param terms Character vector of candidate exclusion terms.
#' @param noise_cov_list List of logical vectors (each of length `n_noise`) indicating
#'   which sampled non-seed articles contain each candidate term.
#' @param rank_scores Numeric vector of candidate importance scores (same length as `terms`).
#' @param noise_mesh Optional list of character vectors of MeSH headings per noise article.
#' @param noise_keywords Optional list of character vectors of author keywords per noise article.
#' @return Character vector of automatically generated concept cluster labels
#'   (same length as `terms`, e.g., `"topic_mice_receptor"`, `"topic_h1n1_h7n9"`).
#' @export
cluster_exclusion_terms <- function(terms,
                                    noise_cov_list,
                                    rank_scores = NULL,
                                    noise_mesh = NULL,
                                    noise_keywords = NULL) {
  m <- length(terms)
  if (m == 0L) {
    return(character(0))
  }
  if (is.null(rank_scores)) rank_scores <- rep(1, m)
  if (m == 1L) {
    slug <- gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", gsub("\\*$", "", terms[1L])))
    return(paste0("topic_", if (nzchar(slug)) substr(slug, 1L, 18L) else "1"))
  }

  n_noise <- length(noise_cov_list[[1L]])
  if (n_noise == 0L) {
    return(rep("topic_1", m))
  }

  is_code <- grepl("^[a-z]+[0-9]+[a-z0-9]*$", gsub("\\*$", "", tolower(terms)))

  # Build optional MeSH + Author Keyword cosine similarity across noise articles
  meta_sim <- matrix(0, nrow = m, ncol = m)
  has_meta <- (!is.null(noise_mesh) && length(noise_mesh) == n_noise) ||
    (!is.null(noise_keywords) && length(noise_keywords) == n_noise)
  if (has_meta) {
    doc_tags <- lapply(seq_len(n_noise), function(d) {
      mh <- if (!is.null(noise_mesh)) noise_mesh[[d]] %||% character(0) else character(0)
      kw <- if (!is.null(noise_keywords)) noise_keywords[[d]] %||% character(0) else character(0)
      unique(tolower(c(mh, kw)))
    })
    tag_tab <- table(unlist(doc_tags))
    keep_tags <- names(tag_tab[tag_tab >= 2L])
    if (length(keep_tags) > 0L) {
      if (length(keep_tags) > 200L) keep_tags <- names(sort(tag_tab[keep_tags], decreasing = TRUE))[1:200]
      M_mat <- matrix(0, nrow = m, ncol = length(keep_tags))
      colnames(M_mat) <- keep_tags
      for (i in seq_len(m)) {
        hits <- which(noise_cov_list[[i]])
        if (length(hits) > 0L) {
          tb <- table(unlist(doc_tags[hits]))
          cv <- intersect(names(tb), keep_tags)
          if (length(cv) > 0L) M_mat[i, cv] <- as.numeric(tb[cv])
        }
      }
      norms <- sqrt(rowSums(M_mat^2))
      norms[norms < 1e-12] <- 1
      M_norm <- M_mat / norms
      meta_sim <- tcrossprod(M_norm)
    }
  }

  edge_from <- integer(0)
  edge_to <- integer(0)
  edge_w <- numeric(0)

  for (i in seq_len(m - 1L)) {
    ci <- noise_cov_list[[i]]
    pi <- sum(ci) / n_noise
    for (j in (i + 1L):m) {
      cj <- noise_cov_list[[j]]
      pj <- sum(cj) / n_noise
      n_both <- sum(ci & cj)
      n_union <- sum(ci | cj)
      jacc <- if (n_union > 0L) n_both / n_union else 0

      npmi <- 0
      if (n_both > 0L && pi > 0 && pj > 0) {
        pij <- n_both / n_noise
        denom <- -log(max(pij, 1e-9))
        if (denom > 1e-9) {
          npmi <- max(0, log(pij / (pi * pj)) / denom)
        }
      }

      stem_eq <- shares_stem(terms[i], terms[j]) || .is_code_family(terms[i], terms[j])
      code_aff <- is_code[i] && is_code[j]

      w <- 0.45 * npmi + 0.35 * jacc + 0.20 * meta_sim[i, j] +
        (if (stem_eq) 0.35 else 0) + (if (code_aff) 0.28 else 0)

      if (w >= 0.18) {
        edge_from <- c(edge_from, i)
        edge_to <- c(edge_to, j)
        edge_w <- c(edge_w, w)
      }
    }
  }

  g <- igraph::make_empty_graph(n = m, directed = FALSE)
  node_str <- rep(0, m)
  if (length(edge_w) > 0L) {
    g <- igraph::add_edges(g, as.vector(rbind(edge_from, edge_to)), attr = list(weight = edge_w))
    comm <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
    mem <- as.integer(igraph::membership(comm))
    raw_str <- as.numeric(igraph::strength(g, weights = igraph::E(g)$weight))
    if (max(raw_str) > 0) node_str <- raw_str / max(raw_str)
  } else {
    mem <- seq_len(m)
  }

  # Consolidate tiny isolated singleton clusters if > 5 clusters exist so UI stays concise
  c_sizes <- table(mem)
  c_max_score <- tapply(rank_scores, mem, max)
  ord_c <- as.integer(names(sort(c_max_score, decreasing = TRUE)))
  top_c <- utils::head(ord_c[c_sizes[as.character(ord_c)] >= 2L | seq_along(ord_c) <= 4L], 5L)

  labels <- character(m)
  used_lbls <- character(0)
  for (cid in top_c) {
    grp <- which(mem == cid)
    grp_ord <- grp[order(rank_scores[grp], decreasing = TRUE)]
    # Pick up to 2 stem-distinct representative words to build a descriptive topic slug
    rep_words <- character(0)
    rep_stems <- character(0)
    for (gi in grp_ord) {
      w_clean <- gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", gsub("\\*$", "", tolower(terms[gi]))))
      w_tok <- strsplit(w_clean, "_", fixed = TRUE)[[1L]][1L]
      st <- term_stems(w_tok)
      st_key <- if (length(st) > 0L) st[1L] else w_tok
      if (nzchar(w_tok) && !(st_key %in% rep_stems)) {
        rep_words <- c(rep_words, substr(w_tok, 1L, 12L))
        rep_stems <- c(rep_stems, st_key)
      }
      if (length(rep_words) >= 2L) break
    }
    base_lbl <- paste0("topic_", if (length(rep_words) > 0L) paste(rep_words, collapse = "_") else cid)
    lbl <- base_lbl
    sfx <- 2L
    while (lbl %in% used_lbls) {
      lbl <- paste0(base_lbl, "_", sfx)
      sfx <- sfx + 1L
    }
    used_lbls <- c(used_lbls, lbl)
    labels[grp] <- lbl
  }

  leftover <- which(!nzchar(labels))
  if (length(leftover) > 0L) {
    labels[leftover] <- "topic_other"
  }
  attr(labels, "node_strength") <- node_str
  labels
}

# Unsupervised prefix wildcard synthesis across sampled non-seed articles:
# finds 5..8-char prefixes matching >= 2 distinct noise words with 0 seed loss in tiab or title
.synthesize_exclusion_wildcards <- function(doc_grams, noise_texts, noise_titles,
                                            seed_texts, seed_titles, min_noise_cnt, stops) {
  all_u <- unique(unlist(doc_grams))
  is_gen <- if (exists(".is_generic_academic_prefix", mode = "function")) {
    .is_generic_academic_prefix(all_u)
  } else {
    rep(FALSE, length(all_u))
  }
  u <- all_u[grepl("^[a-z]{4,}$", all_u) & !(all_u %in% stops) & !is_gen]
  if (length(u) < 2L) {
    return(character(0))
  }

  u_counts <- table(unlist(lapply(doc_grams, function(g) intersect(g, u))))
  cand_prefixes <- unique(unlist(lapply(5:12, function(len) {
    long_u <- u[nchar(u) >= len]
    if (length(long_u) == 0L) character(0) else substr(long_u, 1L, len)
  })))
  is_gen_p <- if (exists(".is_generic_academic_prefix", mode = "function")) {
    .is_generic_academic_prefix(cand_prefixes)
  } else {
    rep(FALSE, length(cand_prefixes))
  }
  cand_prefixes <- cand_prefixes[!(cand_prefixes %in% stops) & !is_gen_p]
  if (length(cand_prefixes) == 0L) {
    return(character(0))
  }

  matched_sets <- lapply(cand_prefixes, function(p) sort(u[startsWith(u, p)]))
  keep_p <- vapply(seq_along(cand_prefixes), function(i) {
    p <- cand_prefixes[i]
    m_u <- matched_sets[[i]]
    if (length(m_u) < 2L) {
      return(FALSE)
    }
    sfx <- substring(m_u, nchar(p) + 1L)
    if (min(nchar(sfx)) > 4L) {
      return(FALSE)
    }
    if (nchar(p) == 5L) {
      ends_classical <- substr(p, 5L, 5L) %in% c("o", "i") && min(nchar(sfx)) <= 3L
      has_std_sfx <- any(sfx %in% c("", "e", "s", "ed", "al", "ar", "ic", "or", "er", "y", "ia", "is", "us", "um", "on"))
      if (!ends_classical && !has_std_sfx) {
        return(FALSE)
      }
    }
    TRUE
  }, logical(1))
  cand_prefixes <- cand_prefixes[keep_p]
  matched_sets <- matched_sets[keep_p]
  if (length(cand_prefixes) == 0L) {
    return(character(0))
  }

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

    loss_tiab <- sum(term_matches_text(w_term, seed_texts))
    loss_ti <- sum(term_matches_text(w_term, seed_titles))
    if (loss_tiab > 0L && loss_ti > 0L) next

    hits_tiab <- sum(term_matches_text(w_term, noise_texts))
    hits_ti <- sum(term_matches_text(w_term, noise_titles))
    eff_hits <- if (loss_tiab == 0L) hits_tiab else hits_ti

    if ((loss_tiab == 0L && hits_tiab >= min_noise_cnt && (hits_tiab > max_single || length(m_u) >= 3L)) ||
      (loss_tiab > 0L && loss_ti == 0L && hits_ti >= max(3L, floor(min_noise_cnt * 0.75)))) {
      if (exists(".register_wildcard_expansion", mode = "function")) {
        m_cnts <- as.integer(u_counts[m_u])
        best_w <- m_u[order(-m_cnts, nchar(m_u))][1L]
        .register_wildcard_expansion(p, best_w)
      }
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

#' Discover unsupervised, concept-clustered exclusion (`NOT` / `-`) terms from non-seed PubMed results
#'
#' Samples articles returned by `base_query` that are NOT in `seeds_df` across the full
#' publication timeline, mines frequent 1- and 2-grams plus morphological prefix wildcards
#' (`*`), checks seed loss locally against PubMed-resolved seeds (including PubMed Author
#' Keywords), partitions candidates into unsupervised noise topic clusters via
#' [cluster_exclusion_terms()], builds greedy set-cover exclusion blocks, and returns up
#' to `max_hints` ranked exclusion hints.
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
#' @export
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
  if (is.null(base_query) || length(base_query$blocks) == 0L) {
    return(.empty_exclusion_hints())
  }

  pm_base <- to_pubmed_query(base_query)
  if (!nzchar(pm_base)) {
    return(.empty_exclusion_hints())
  }

  seed_pmids <- unique(as.character(seeds_df$pmid[nzchar(seeds_df$pmid)]))
  retmax <- min(max(as.integer(sample_n) * 4L + length(seed_pmids), 600L), 2500L)

  if (is.function(progress_cb)) progress_cb("Sampling non-seed PubMed results for exclusion hints...")
  res_ids <- pm_esearch_fn(
    term = pm_base, retmax = retmax,
    email = email, api_key = ncbi_key, cache_dir = cache_dir, use_cache = use_cache
  )
  if (!isTRUE(res_ids$ok) || length(res_ids$idlist) == 0L) {
    return(.empty_exclusion_hints())
  }

  if (is.na(base_pm_hits)) base_pm_hits <- as.integer(res_ids$count %||% length(res_ids$idlist))
  noise_all <- setdiff(res_ids$idlist, seed_pmids)
  if (length(noise_all) < 5L) {
    return(.empty_exclusion_hints())
  }

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
  if (is.null(noise_meta) || nrow(noise_meta) < 5L) {
    return(.empty_exclusion_hints())
  }

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
    "subsequent", "subsequently", "exhibit", "exhibited", "exhibiting", "activate", "activated",
    "derive", "derived", "attribute", "attributed", "attributable", "determine", "determined",
    "represent", "represented", "confirm", "confirmed", "observe", "observed",
    "broad", "broadly", "diverse", "diversity", "robust", "robustly", "single", "multiple",
    "notable", "notably", "integrate", "integrated", "novel", "current", "common", "distinct",
    "major", "minor", "rapid", "effective", "critical", "essential", "initial", "final",
    "individual", "typical", "widely", "largely", "mostly", "strongly", "clearly", "directly",
    "design", "designed", "scale", "scales", "knowledge", "context", "content", "combine",
    "combined", "combining", "responsible", "responsibility", "function", "functional", "commercial"
  ))
  pos_terms <- unique(tolower(gsub("\\*", "", unlist(lapply(base_query$blocks, `[[`, "terms")))))
  pos_stems <- unique(unlist(lapply(pos_terms, term_stems)))

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

  # Include parenthesized compound alphanumeric codes (e.g. "A(H1N1)pdm09" -> "h1n1pdm09")
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

  cands <- Filter(function(t_str) {
    t_bare <- gsub("\\*", "", t_str)
    if (t_bare %in% pos_terms) {
      return(FALSE)
    }
    if (any(vapply(pos_terms, .is_code_family, logical(1), t2 = t_bare))) {
      return(FALSE)
    }
    if (any(vapply(pos_terms, function(pt) grepl(paste0("(?<![a-z0-9])", escape_regex(pt), "(?![a-z0-9])"), t_bare, perl = TRUE), logical(1)))) {
      return(FALSE)
    }
    if (length(intersect(term_stems(t_bare), pos_stems)) > 0L) {
      return(FALSE)
    }
    if (nchar(t_bare) >= 5L && any(vapply(pos_terms[nchar(pos_terms) >= 5L], function(pt) {
      identical(substr(pt, 1L, 5L), substr(t_bare, 1L, 5L)) || .char_ngram_sim(pt, t_bare, 4L) >= 0.40
    }, logical(1)))) {
      return(FALSE)
    }
    TRUE
  }, cands)
  if (length(cands) == 0L) {
    return(.empty_exclusion_hints())
  }

  in_gram <- intersect(cands, names(gram_df))
  cands_ord <- in_gram[order(gram_df[in_gram], decreasing = TRUE)]
  cands <- unique(c(intersect(wildcard_cands, cands), utils::head(cands_ord, 240L), intersect(cands, code_tokens)))

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
    scope_bonus <- if (use_scope$field == "tiab" && n_lost == 0L) 1.35 else 1.0
    code_bonus <- if (is_code && use_scope$field == "tiab" && n_lost == 0L) 1.50 else 1.0
    wild_bonus <- if (is_wild && n_lost == 0L) 1.35 else 1.0
    unigram_bonus <- if (!grepl(" ", t_str, fixed = TRUE)) 1.15 else 1.0

    scored[[length(scored) + 1L]] <- list(
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
      rank_score = (est_rem * scope_bonus * code_bonus * wild_bonus * unigram_bonus) / (1 + 6 * n_lost)
    )
  }

  if (length(scored) == 0L) {
    return(.empty_exclusion_hints())
  }

  # Deduplicate literal words already subsumed by a safe wildcard in the same scope
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
  if (length(scored_ord) > 65L) scored_ord <- scored_ord[1:65]

  # Unsupervised noise co-occurrence graph clustering (`cluster_exclusion_terms`)
  cpt_labels <- cluster_exclusion_terms(
    terms = vapply(scored_ord, `[[`, character(1), "term"),
    noise_cov_list = lapply(scored_ord, `[[`, "noise_cov"),
    rank_scores = vapply(scored_ord, `[[`, numeric(1), "rank_score"),
    noise_mesh = noise_meta$pm_mesh,
    noise_keywords = noise_meta$pm_keywords
  )
  node_str_vec <- attr(cpt_labels, "node_strength") %||% rep(0.5, length(scored_ord))
  cpt_labels_plain <- as.character(cpt_labels)
  for (i in seq_along(scored_ord)) {
    scored_ord[[i]]$concept <- cpt_labels_plain[i]
    scored_ord[[i]]$node_strength <- node_str_vec[i]
  }
  cpt_labels <- cpt_labels_plain

  # Select diverse candidates across all discovered noise topic clusters (up to max_hints)
  code_idx <- which(vapply(scored_ord, function(x) isTRUE(x$is_code) && x$seeds_lost_n == 0L, logical(1)))
  wild_idx <- which(vapply(scored_ord, function(x) isTRUE(x$is_wild) && x$seeds_lost_n == 0L, logical(1)))
  all_cpts <- unique(cpt_labels)
  concept_picks <- integer(0)
  for (cpt in all_cpts) {
    c_rows <- which(cpt_labels == cpt)
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

  spec_map <- stats::setNames(lapply(sel_cands, `[[`, "noise_cov"), vapply(sel_cands, `[[`, character(1), "exclude_spec"))
  coh_map <- stats::setNames(
    vapply(sel_cands, function(x) as.numeric(x$node_strength %||% 0.5), numeric(1)),
    vapply(sel_cands, `[[`, character(1), "exclude_spec")
  )
  safe_rows <- which(out_df$seeds_lost_n == 0L & out_df$exact_hits_removed > 0L)

  # 1. Build Concept-Grouped Exclusion Blocks (greedy set-cover of safe exclusions within each discovered noise topic)
  exclusion_blocks <- list()
  out_cpts <- unique(out_df$concept[safe_rows])
  for (cpt in out_cpts) {
    c_idx <- safe_rows[out_df$concept[safe_rows] == cpt]
    if (length(c_idx) == 0L) next
    cov_cur <- rep(FALSE, n_noise)
    picked_specs <- character(0)
    while (length(picked_specs) < 5L) {
      rem <- setdiff(out_df$exclude_spec[c_idx], picked_specs)
      if (length(rem) == 0L) break
      gains <- vapply(rem, function(sp) {
        cv <- spec_map[[sp]]
        if (is.null(cv)) {
          return(0L)
        }
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
  if (length(exclusion_blocks) > 1L) {
    eb_ord <- order(vapply(exclusion_blocks, `[[`, integer(1), "noise_covered_n"), decreasing = TRUE)
    exclusion_blocks <- exclusion_blocks[eb_ord]
  }

  # 2. Build global greedy set-cover exclusion sequence across all safe exclusion candidates
  global_setcover <- character(0)
  if (length(safe_rows) > 0L) {
    cov_g <- rep(FALSE, n_noise)
    avail_sp <- out_df$exclude_spec[safe_rows]
    while (length(global_setcover) < min(18L, length(avail_sp))) {
      rem <- setdiff(avail_sp, global_setcover)
      if (length(rem) == 0L) break
      gains <- vapply(rem, function(sp) {
        cv <- spec_map[[sp]]
        if (is.null(cv)) {
          return(-1)
        }
        new_n <- sum(cv & !cov_g)
        if (new_n == 0L) {
          return(-1)
        }
        row_i <- match(sp, out_df$exclude_spec)
        gs_len <- nchar(out_df$gs_clause[row_i])
        t_sp <- out_df$term[row_i]
        rem_hits <- pmax(as.numeric(out_df$exact_hits_removed[row_i]), 1)
        coh_mult <- 0.45 + 0.55 * (coh_map[[sp]] %||% 0.5)
        type_mult <- if (grepl("^[a-z]+[0-9]+[a-z0-9]*$", t_sp)) 1.25 else if (endsWith(t_sp, "*")) 1.15 else 1.0
        ((new_n^0.75) * (rem_hits^0.30) * coh_mult * type_mult) / (gs_len^0.35)
      }, numeric(1))
      if (max(gains) <= 0) break
      best_sp <- rem[which.max(gains)]
      global_setcover <- c(global_setcover, best_sp)
      cov_g <- cov_g | spec_map[[best_sp]]
    }
  }

  # 3. Build sibling alphanumeric code exclusion bundle (top 2 verified zero-loss tiab codes + prefix extensions)
  verified_safe_codes <- out_df$term[
    out_df$seeds_lost_n == 0L &
      out_df$field == "tiab" &
      grepl("^[a-z]+[0-9]+[a-z0-9]*$", out_df$term) &
      out_df$exact_hits_removed > 0L
  ]
  if (length(verified_safe_codes) > 0L) {
    top2_codes <- utils::head(verified_safe_codes, 2L)
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

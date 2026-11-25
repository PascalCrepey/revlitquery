# R/04_term_extraction.R
# N-gram extraction, academic stopword filtering, PubMed background lift scoring,
# and stem-disjoint automatic concept block discovery.

#' Stopword list (standard English + academic/methods filler)
get_stopwords <- function() {
  unique(c(
    # Standard English
    "a", "an", "and", "are", "as", "at", "be", "been", "being", "but", "by",
    "can", "could", "did", "do", "does", "doing", "done", "each", "either",
    "else", "for", "from", "further", "had", "has", "have", "having", "he",
    "her", "here", "hers", "herself", "him", "himself", "his", "how", "i",
    "if", "in", "into", "is", "it", "its", "itself", "just", "me", "more",
    "most", "my", "myself", "no", "nor", "not", "now", "of", "off", "on",
    "once", "only", "or", "other", "our", "ours", "ourselves", "out", "over",
    "own", "same", "she", "should", "so", "some", "such", "than", "that",
    "the", "their", "theirs", "them", "themselves", "then", "there", "these",
    "they", "this", "those", "through", "to", "too", "under", "until", "up",
    "very", "was", "we", "were", "what", "when", "where", "which", "while",
    "who", "whom", "why", "will", "with", "would", "you", "your", "yours",
    # Academic / boilerplate stopwords
    "study", "studies", "results", "result", "methods", "method", "background",
    "conclusions", "conclusion", "introduction", "objective", "objectives",
    "purpose", "data", "analysis", "analyses", "using", "used", "use",
    "based", "however", "although", "among", "between", "during", "within",
    "across", "including", "include", "included", "associated", "reported",
    "report", "reports", "found", "showed", "shown", "total", "number",
    "numbers", "high", "higher", "low", "lower", "may", "also", "one",
    "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten",
    "per", "et", "al", "ci", "confidence", "interval", "significantly",
    "significant", "respectively", "copyright", "rights", "reserved",
    "elsevier", "springer", "wiley", "author", "authors", "article", "paper",
    "review", "table", "figure", "group", "groups", "year", "years", "month",
    "months", "day", "days", "time", "median", "mean", "range", "versus", "vs",
    "compared", "compare", "observed", "estimated", "estimate", "estimates",
    "available", "well", "since", "first", "overall", "all", "both", "many",
    "several", "large", "small", "new", "potential", "important", "furthermore",
    "moreover", "therefore", "thus", "hence", "whereas", "whether", "less",
    "least", "much", "about", "after", "before", "above", "below", "against"
  ))
}

#' Normalize text for n-gram extraction (preserves alphanumeric tokens like h5n1)
normalize_search_text <- function(x) {
  x <- tolower(as.character(x %||% ""))
  x[is.na(x)] <- ""
  # Replace hyphens and slashes with spaces so "case-fatality" -> "case fatality"
  x <- stringr::str_replace_all(x, "[-/]+", " ")
  # Strip non-alphanumeric punctuation
  x <- stringr::str_replace_all(x, "[^a-z0-9\\s]", " ")
  stringr::str_squish(x)
}

#' Escape regex metacharacters in literal term
escape_regex <- function(x) {
  stringr::str_replace_all(x, "([.\\^$|()*+?{}\\[\\]\\\\])", "\\\\\\1")
}

#' Check whether a candidate term matches each text in `texts` (case-insensitive, hyphen-flexible)
#'
#' @param term Single candidate term.
#' @param texts Character vector of texts.
#' @return Logical vector of length `length(texts)`.
term_matches_text <- function(term, texts) {
  if (!nzchar(term) || length(texts) == 0L) return(rep(FALSE, length(texts)))
  t_norm <- stringr::str_squish(tolower(stringr::str_replace_all(term, "[-/]+", " ")))
  if (!nzchar(t_norm)) return(rep(FALSE, length(texts)))
  pat <- escape_regex(t_norm)
  pat <- stringr::str_replace_all(pat, "\\\\\\*", "[a-z0-9]*")
  pat <- stringr::str_replace_all(pat, " ", "[\\\\s-]+")
  full_pat <- paste0("(?i)(?<![a-z0-9])", pat, "(?![a-z0-9])")
  texts_norm <- stringr::str_squish(tolower(as.character(texts %||% "")))
  grepl(full_pat, texts_norm, perl = TRUE)
}

#' Simple token stemmer for lexical overlap detection across concept blocks
term_stems <- function(term) {
  stops <- get_stopwords()
  w <- unlist(stringr::str_split(normalize_search_text(term), " "))
  w <- w[nzchar(w) & !(w %in% stops) & nchar(w) >= 2L]
  if (length(w) == 0L) return(character(0))
  # Light suffix stripping (careful with short acronyms containing digits)
  stems <- vapply(w, function(tok) {
    if (grepl("[0-9]", tok) || nchar(tok) <= 4L) return(tok)
    tok <- sub("alities$", "al", tok)
    tok <- sub("ality$", "al", tok)
    tok <- sub("ally$", "al", tok)
    tok <- sub("ities$", "y", tok)
    tok <- sub("ity$", "", tok)
    tok <- sub("ies$", "y", tok)
    tok <- sub("ing$", "", tok)
    tok <- sub("ed$", "", tok)
    tok <- sub("es$", "", tok)
    tok <- sub("s$", "", tok)
    if (nchar(tok) >= 5L) substr(tok, 1L, min(nchar(tok), 6L)) else tok
  }, character(1), USE.NAMES = FALSE)
  unique(stems)
}

#' Check whether two terms share an alphanumeric subtype/code prefix (e.g. h5n1 and h5)
.is_code_family <- function(t1, t2) {
  t1 <- tolower(stringr::str_squish(as.character(t1 %||% "")))
  t2 <- tolower(stringr::str_squish(as.character(t2 %||% "")))
  if (!grepl("^[a-z]+[0-9]+", t1) || !grepl("^[a-z]+[0-9]+", t2)) return(FALSE)
  p1 <- regmatches(t1, regexpr("^[a-z]+[0-9]+", t1))
  p2 <- regmatches(t2, regexpr("^[a-z]+[0-9]+", t2))
  identical(p1, p2) && (startsWith(t1, t2) || startsWith(t2, t1))
}

#' Check whether term `t_broad` lexically subsumes `t_narrow` inside an OR block
.term_subsumes <- function(t_broad, t_narrow) {
  t_broad <- tolower(stringr::str_squish(as.character(t_broad %||% "")))
  t_narrow <- tolower(stringr::str_squish(as.character(t_narrow %||% "")))
  if (!nzchar(t_broad) || !nzchar(t_narrow) || identical(t_broad, t_narrow)) return(FALSE)
  pat <- escape_regex(stringr::str_squish(stringr::str_replace_all(t_broad, "[-/]+", " ")))
  pat <- stringr::str_replace_all(pat, "\\\\\\*", "[a-z0-9]*")
  pat <- stringr::str_replace_all(pat, " ", "[\\\\s-]+")
  t_narrow_norm <- stringr::str_squish(stringr::str_replace_all(t_narrow, "[-/]+", " "))
  grepl(paste0("(?<![a-z0-9])", pat, "(?![a-z0-9])"), t_narrow_norm, perl = TRUE)
}

#' Check whether two terms share any significant token stem
shares_stem <- function(t1, t2) {
  s1 <- term_stems(t1)
  s2 <- term_stems(t2)
  length(intersect(s1, s2)) > 0L
}

# Extract 1..3-grams from a single normalized text
.extract_doc_ngrams <- function(text, max_n = 3L, stops = get_stopwords()) {
  if (!nzchar(text)) return(character(0))
  toks <- unlist(stringr::str_split(text, " "))
  toks <- toks[nzchar(toks)]
  nt <- length(toks)
  if (nt == 0L) return(character(0))

  grams <- character(0)
  for (n_g in seq_len(min(max_n, nt))) {
    if (n_g == 1L) {
      keep <- !(toks %in% stops) & nchar(toks) >= 2L & !grepl("^[0-9]+$", toks)
      grams <- c(grams, toks[keep])
    } else {
      idx <- seq_len(nt - n_g + 1L)
      first_w <- toks[idx]
      last_w <- toks[idx + n_g - 1L]
      # Do not start or end with a stopword or pure number
      valid <- !(first_w %in% stops) & !(last_w %in% stops) &
        !grepl("^[0-9]+$", first_w) & !grepl("^[0-9]+$", last_w)
      if (n_g == 3L) {
        mid_w <- toks[idx + 1L]
        mid_bad <- c("with", "of", "in", "on", "for", "and", "or", "the", "a", "an", "to", "by", "at", "as", "from", "into")
        valid <- valid & !(mid_w %in% mid_bad)
      }
      if (any(valid)) {
        cand <- vapply(idx[valid], function(k) {
          paste(toks[k:(k + n_g - 1L)], collapse = " ")
        }, character(1))
        grams <- c(grams, cand)
      }
    }
  }
  unique(grams)
}

# Detect uppercase acronyms in original un-lowercased text
.detect_acronyms <- function(raw_texts) {
  matches <- stringr::str_extract_all(paste(raw_texts, collapse = " "), "\\b[A-Z][A-Z0-9]{1,9}\\b")
  unique(tolower(unlist(matches)))
}

#' Classify a candidate search term into a coherent semantic PICO family
#'
#' @param term Single candidate term.
#' @param block1_coupled Optional character vector of terms coupled to Block 1.
#' @return One of `"primary_topic"`, `"outcome_severity"`, `"population_host"`,
#'   `"epi_methods"`, or `"other"`.
assign_term_family <- function(term, block1_coupled = character(0)) {
  t_raw <- tolower(stringr::str_squish(as.character(term %||% "")))
  if (!nzchar(t_raw)) return("other")
  t_low <- gsub("\\*", "", t_raw)
  if (length(block1_coupled) > 0L) {
    if (t_raw %in% block1_coupled || t_low %in% block1_coupled) return("primary_topic")
    toks <- unlist(strsplit(t_low, "\\s+"))
    if (length(toks) > 1L && any(toks %in% block1_coupled)) return("primary_topic")
  }
  if (grepl("\\b(case fatality|fatality|fatal[a-z]*|mortalit[a-z]*|morbidit[a-z]*|death|deaths|died|dying|deadly|lethal[a-z]*|cfr|ifr|hfr|surviv[a-z]*|severit[a-z]*|severe disease)\\b", t_low, perl = TRUE)) {
    return("outcome_severity")
  }
  if (grepl("\\b(human|humans|patient|patients|children|pediatric|paediatric|adult|adults|workers|individuals|people|persons|confirmed cases|human cases|cases)\\b", t_low, perl = TRUE)) {
    return("population_host")
  }
  if (grepl("\\b(seroprevalence|seroprev[a-z]*|serologic[a-z]*|serol[a-z]*|serosurvey[a-z]*|seroepidemiolog[a-z]*|epidemiolog[a-z]*|epidemiol[a-z]*|surveillance|outbreak|outbreaks|transmission|transmi[a-z]*|incidence|prevalence|cluster|clusters|virulence|virulen[a-z]*|zoonot[a-z]*|zoonos[a-z]*|model|modeling|modelling)\\b", t_low, perl = TRUE)) {
    return("epi_methods")
  }
  "other"
}

#' Vectorized wrapper for `assign_term_family()`
assign_term_families <- function(terms, block1_coupled = character(0)) {
  if (length(terms) == 0L) return(stats::setNames(character(0), character(0)))
  fams <- vapply(terms, assign_term_family, character(1), block1_coupled = block1_coupled, USE.NAMES = FALSE)
  stats::setNames(fams, terms)
}

# Synthesize morphological prefix wildcards (e.g. fatal*, mortalit*, epidemiol*, human*)
# for PICO terms when >= 2 corpus unigrams share a 5..8-char prefix and the wildcard rescues >= 2 seeds
.synthesize_prefix_wildcards <- function(doc_grams, raw_texts, min_df = 2L, stops = get_stopwords()) {
  k <- length(raw_texts)
  if (k < 6L) return(character(0))
  all_u <- unique(unlist(doc_grams))
  u <- all_u[grepl("^[a-z]{4,}$", all_u) & !(all_u %in% stops)]
  if (length(u) < 2L) return(character(0))
  # Restrict to PICO vocabulary so generic verbs/nouns (e.g. result/resulting) are not wildcarded
  u_fams <- assign_term_families(u)
  u <- u[u_fams %in% c("outcome_severity", "population_host", "epi_methods")]
  if (length(u) < 2L) return(character(0))

  u_counts <- table(unlist(lapply(doc_grams, function(g) intersect(g, u))))
  cand_prefixes <- unique(unlist(lapply(5:8, function(len) {
    long_u <- u[nchar(u) >= len]
    if (length(long_u) == 0L) character(0) else substr(long_u, 1L, len)
  })))
  cand_prefixes <- cand_prefixes[!(cand_prefixes %in% stops)]
  if (length(cand_prefixes) == 0L) return(character(0))

  # Group prefixes by the exact set of corpus unigrams they match; keep longest prefix per set
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

  kept_wildcards <- character(0)
  kept_roots <- character(0)
  kept_dfs <- integer(0)

  for (i in seq_along(cand_prefixes)) {
    p <- cand_prefixes[i]
    m_u <- matched_sets[[i]]
    max_single_df <- max(as.integer(u_counts[m_u]), na.rm = TRUE)
    w_term <- paste0(p, "*")
    df_star <- sum(term_matches_text(w_term, raw_texts))
    if (df_star >= max(3L, min_df) && (df_star - max_single_df) >= 2L) {
      root5 <- substr(p, 1L, 5L)
      existing_idx <- match(root5, kept_roots)
      if (is.na(existing_idx)) {
        kept_wildcards <- c(kept_wildcards, w_term)
        kept_roots <- c(kept_roots, root5)
        kept_dfs <- c(kept_dfs, df_star)
      } else if (df_star > kept_dfs[existing_idx] ||
                 (df_star == kept_dfs[existing_idx] && nchar(w_term) > nchar(kept_wildcards[existing_idx]))) {
        kept_wildcards[existing_idx] <- w_term
        kept_dfs[existing_idx] <- df_star
      }
    }
  }
  kept_wildcards
}

#' Extract and score candidate search terms from enriched seed papers
#'
#' @param seeds_df Enriched tibble from `enrich_seeds()`.
#' @param config Optional config list (`banned_terms`, `required_terms`, `custom_blocks`).
#' @param max_candidates Max terms to probe against PubMed background frequency.
#' @param probe_pubmed Logical; whether to probe PubMed for background `bg_hits`.
#' @param eval_bg_fn Optional function `(term) -> integer` for testing without network.
#' @return Tibble of scored candidate terms with `cov_tiab` and `cov_title` list-columns.
extract_candidate_terms <- function(seeds_df,
                                    config = list(),
                                    max_candidates = 45L,
                                    probe_pubmed = TRUE,
                                    eval_bg_fn = NULL,
                                    email = NULL,
                                    ncbi_key = NULL,
                                    cache_dir = ".cache",
                                    use_cache = TRUE,
                                    progress_cb = NULL) {
  k <- nrow(seeds_df)
  stops <- get_stopwords()
  banned <- tolower(stringr::str_squish(as.character(config$banned_terms %||% character(0))))
  required <- tolower(stringr::str_squish(as.character(config$required_terms %||% character(0))))
  custom_terms <- tolower(stringr::str_squish(unlist(config$custom_blocks %||% list())))
  pinned_terms <- unique(c(required, custom_terms))
  pinned_terms <- pinned_terms[nzchar(pinned_terms)]

  raw_titles <- as.character(seeds_df$full_title %||% seeds_df$title)
  raw_texts <- as.character(seeds_df$full_text %||% paste(raw_titles, seeds_df$abstract))
  norm_texts <- normalize_search_text(raw_texts)
  acronyms <- .detect_acronyms(raw_texts)

  # Collect per-document n-grams
  doc_grams <- lapply(norm_texts, .extract_doc_ngrams, max_n = 3L, stops = stops)

  # Also add cleaned author keywords and MeSH terms that appear in seed full_text
  extra_terms <- unique(c(
    unlist(seeds_df$all_keywords %||% list()),
    unlist(seeds_df$mesh_terms %||% list())
  ))
  extra_terms <- normalize_search_text(extra_terms)
  extra_terms <- extra_terms[nzchar(extra_terms) & !(extra_terms %in% stops) & nchar(extra_terms) >= 2L]

  # Count raw document frequency from n-grams
  all_grams <- unique(unlist(doc_grams))
  all_terms_pool <- unique(c(all_grams, extra_terms, pinned_terms))
  all_terms_pool <- all_terms_pool[nzchar(all_terms_pool) & !(all_terms_pool %in% banned)]
  # Limit extreme unigrams/ngrams before exact regex matching if huge corpus
  gram_counts <- table(unlist(doc_grams))
  min_df <- if (k <= 4L) 1L else 2L

  frequent_grams <- names(gram_counts[gram_counts >= min_df])
  prefix_wildcards <- .synthesize_prefix_wildcards(doc_grams, raw_texts, min_df = min_df, stops = stops)
  # Guarantee each seed contributes at least its top title/text n-grams and any PICO outcome terms (even if df==1)
  per_seed_rescue <- unique(unlist(lapply(seq_len(k), function(i) {
    t_grams <- .extract_doc_ngrams(normalize_search_text(raw_titles[i]), max_n = 3L, stops = stops)
    d_grams <- doc_grams[[i]]
    out_grams <- d_grams[assign_term_families(d_grams) == "outcome_severity"]
    c(utils::head(t_grams, 6L), utils::head(out_grams, 4L))
  })))

  keep_pool <- unique(c(frequent_grams, prefix_wildcards, per_seed_rescue, extra_terms, pinned_terms))
  keep_pool <- keep_pool[nzchar(keep_pool) & !(keep_pool %in% banned) & !(keep_pool %in% stops)]
  # Drop single-character or pure-digit terms
  keep_pool <- keep_pool[nchar(keep_pool) >= 2L & !grepl("^[0-9\\s]+$", keep_pool)]

  # Compute exact coverage across titles and tiab
  cov_tiab_list <- lapply(keep_pool, term_matches_text, texts = raw_texts)
  cov_title_list <- lapply(keep_pool, term_matches_text, texts = raw_titles)

  df_tiab <- vapply(cov_tiab_list, sum, integer(1))
  df_title <- vapply(cov_title_list, sum, integer(1))

  # Filter terms with 0 matches unless pinned
  keep_idx <- (df_tiab >= 1L) | (keep_pool %in% pinned_terms)
  keep_pool <- keep_pool[keep_idx]
  cov_tiab_list <- cov_tiab_list[keep_idx]
  cov_title_list <- cov_title_list[keep_idx]
  df_tiab <- df_tiab[keep_idx]
  df_title <- df_title[keep_idx]

  n_words <- vapply(strsplit(keep_pool, " "), length, integer(1))
  is_acro <- grepl("[0-9]", keep_pool) | (keep_pool %in% acronyms)
  is_wildcard <- endsWith(keep_pool, "*")

  # Local pre-score to shortlist top `max_candidates` before background probing
  local_score <- (df_tiab / k) *
    (1 + 0.5 * (df_title / pmax(df_tiab, 1L))) *
    (1 + 0.25 * (n_words > 1L)) *
    (1 + 0.25 * is_acro) *
    (1 + 0.20 * is_wildcard)

  # Also retain top 3 terms per seed paper so every paper has rescue candidates
  rescue_idx <- unique(unlist(lapply(seq_len(k), function(i) {
    covers_i <- which(vapply(cov_tiab_list, `[`, logical(1), i))
    if (length(covers_i) == 0L) return(integer(0))
    ord <- covers_i[order(local_score[covers_i], decreasing = TRUE)]
    utils::head(ord, 3L)
  })))

  # Enforce stem-diversity cap AND preserve stem-diverse representatives from each PICO semantic family
  ord_local <- order(local_score, decreasing = TRUE)
  fams_local <- assign_term_families(keep_pool)
  pico_idx <- integer(0)
  for (fam in c("outcome_severity", "population_host", "epi_methods")) {
    fam_rows <- which(fams_local == fam & df_tiab >= min_df)
    if (length(fam_rows) > 0L) {
      fam_ord <- fam_rows[order(local_score[fam_rows], decreasing = TRUE)]
      fam_sel <- integer(0)
      fam_st_cnt <- integer(0)
      for (fi in fam_ord) {
        if (length(fam_sel) >= 8L) break
        st_f <- term_stems(keep_pool[fi])
        pst <- if (length(st_f) > 0L) st_f[1L] else keep_pool[fi]
        cc <- if (pst %in% names(fam_st_cnt)) fam_st_cnt[[pst]] else 0L
        if (cc < 2L) {
          fam_sel <- c(fam_sel, fi)
          fam_st_cnt[[pst]] <- cc + 1L
        }
      }
      pico_idx <- c(pico_idx, fam_sel)
    }
  }

  # Specifically rescue acronyms/subtypes for seeds missed by the top primary term,
  # and PICO terms for seeds missed by the top outcome_severity term
  uncov_rescue_idx <- integer(0)
  if (length(ord_local) > 0L) {
    top1_cov <- cov_tiab_list[[ord_local[1L]]]
    uncov1 <- which(!top1_cov)
    if (length(uncov1) > 0L) {
      c_u1 <- which(vapply(cov_tiab_list, function(cv) any(cv[uncov1]), logical(1)) & df_tiab >= min_df)
      if (length(c_u1) > 0L) {
        u1_hits <- vapply(cov_tiab_list[c_u1], function(cv) sum(cv[uncov1]), integer(1))
        u1_score <- u1_hits * (1 + 3 * is_acro[c_u1]) * (1 + (df_title[c_u1] > 0L)) / log10(df_tiab[c_u1] + 2)
        uncov_rescue_idx <- c(uncov_rescue_idx, utils::head(c_u1[order(u1_score, decreasing = TRUE)], 4L))
      }
    }
    out_rows <- which(fams_local == "outcome_severity")
    if (length(out_rows) > 0L) {
      top_out <- out_rows[which.max(local_score[out_rows])]
      uncov_out <- which(!cov_tiab_list[[top_out]])
      if (length(uncov_out) > 0L) {
        c_uo <- which(
          (fams_local %in% c("outcome_severity", "epi_methods")) &
            vapply(cov_tiab_list, function(cv) any(cv[uncov_out]), logical(1)) &
            (df_tiab >= min_df | fams_local == "outcome_severity")
        )
        if (length(c_uo) > 0L) {
          uo_hits <- vapply(cov_tiab_list[c_uo], function(cv) sum(cv[uncov_out]), integer(1))
          uo_score <- uo_hits * (1 + 0.8 * (fams_local[c_uo] == "outcome_severity")) *
            (1 + 0.5 * (n_words[c_uo] > 1L)) / log10(df_tiab[c_uo] + 2)
          uncov_rescue_idx <- c(uncov_rescue_idx, utils::head(c_uo[order(uo_score, decreasing = TRUE)], 8L))
        }
      }
    }
  }

  diverse_idx <- integer(0)
  stem_counts <- integer(0)
  for (idx in ord_local) {
    if (length(diverse_idx) >= max_candidates) break
    st <- term_stems(keep_pool[idx])
    if (length(st) == 0L) next
    primary_st <- st[1L]
    cur_cnt <- if (primary_st %in% names(stem_counts)) stem_counts[[primary_st]] else 0L
    if (cur_cnt < 4L) {
      diverse_idx <- c(diverse_idx, idx)
      stem_counts[[primary_st]] <- cur_cnt + 1L
    }
  }
  top_idx <- unique(c(uncov_rescue_idx, pico_idx, diverse_idx, utils::head(ord_local, max_candidates)))
  top_idx <- utils::head(top_idx, max_candidates + 5L)
  pin_idx <- which(keep_pool %in% pinned_terms)
  sel_idx <- unique(c(pin_idx, top_idx, rescue_idx))
  sel_idx <- sel_idx[seq_len(min(length(sel_idx), max_candidates + 15L))]

  terms_sel <- keep_pool[sel_idx]
  cov_tiab_sel <- cov_tiab_list[sel_idx]
  cov_title_sel <- cov_title_list[sel_idx]
  df_tiab_sel <- df_tiab[sel_idx]
  df_title_sel <- df_title[sel_idx]
  n_words_sel <- n_words[sel_idx]
  is_acro_sel <- is_acro[sel_idx]

  # Probe background PubMed frequency
  probed_live <- is.null(eval_bg_fn) && isTRUE(probe_pubmed)
  if (is.null(eval_bg_fn)) {
    eval_bg_fn <- function(t_str) {
      res <- pubmed_esearch(
        term = paste0('"', t_str, '"[tiab]'), retmax = 0L,
        email = email, api_key = ncbi_key, cache_dir = cache_dir, use_cache = use_cache
      )
      if (isTRUE(res$ok)) res$count else NA_integer_
    }
  }

  bg_hits <- rep(NA_integer_, length(terms_sel))
  if (isTRUE(probe_pubmed) && length(terms_sel) > 0L) {
    for (idx in seq_along(terms_sel)) {
      t_str <- terms_sel[idx]
      if (is.function(progress_cb)) {
        progress_cb(sprintf("Probing term %d/%d: %s", idx, length(terms_sel), t_str))
      }
      val <- suppressWarnings(as.integer(eval_bg_fn(t_str)))
      bg_hits[idx] <- if (length(val) != 1L || is.na(val)) NA_integer_ else val
    }
  }
  # Fallback estimate if probing failed or disabled
  bg_missing <- is.na(bg_hits)
  if (any(bg_missing)) {
    bg_hits[bg_missing] <- ifelse(n_words_sel[bg_missing] >= 3L, 2000L,
                                  ifelse(n_words_sel[bg_missing] == 2L, 25000L, 200000L))
  }

  n_pubmed_corpus <- 38000000
  seed_prop <- df_tiab_sel / max(k, 1L)
  bg_prop <- (bg_hits + 1) / n_pubmed_corpus
  lift <- seed_prop / bg_prop

  term_score <- seed_prop *
    log10(pmax(lift, 1.1)) *
    (1 + 0.35 * (df_title_sel / pmax(df_tiab_sel, 1L))) *
    (1 + 0.15 * (n_words_sel > 1L)) *
    (1 + 0.15 * is_acro_sel) *
    (1 + 0.20 * endsWith(terms_sel, "*"))

  tbl <- tibble::tibble(
    term = terms_sel,
    n_words = n_words_sel,
    is_acronym = is_acro_sel,
    is_pinned = terms_sel %in% pinned_terms,
    df_tiab = df_tiab_sel,
    df_title = df_title_sel,
    bg_hits = bg_hits,
    seed_prop = seed_prop,
    lift = lift,
    term_score = term_score,
    cov_tiab = cov_tiab_sel,
    cov_title = cov_title_sel
  ) |>
    dplyr::arrange(dplyr::desc(term_score), bg_hits)

  attr(tbl, "probed_live") <- probed_live
  attr(tbl, "probe_meta") <- list(email = email, ncbi_key = ncbi_key, cache_dir = cache_dir, use_cache = use_cache)
  tbl
}

#' Build stem-disjoint, conditionally selective concept blocks from scored candidate terms
#'
#' @param terms_df Scored candidate tibble from `extract_candidate_terms()`.
#' @param seeds_df Enriched seeds tibble.
#' @param config Optional config list (`custom_blocks`, `max_blocks`, `max_terms_per_block`).
#' @param max_blocks Maximum number of concept blocks.
#' @param max_terms_per_block Maximum terms per block.
#' @param probe_selectivity Logical or NULL (defaults to TRUE when `terms_df` was probed live).
#' @param eval_joint_fn Optional function `(t1, t2) -> integer` returning PubMed joint hit count.
#' @return List of block lists `list(name = ..., terms = ...)`.
build_concept_blocks <- function(terms_df,
                                 seeds_df,
                                 config = list(),
                                 max_blocks = NULL,
                                 max_terms_per_block = NULL,
                                 probe_selectivity = NULL,
                                 eval_joint_fn = NULL) {
  max_blocks <- as.integer(max_blocks %||% config$max_blocks %||% 3L)
  max_terms_per_block <- as.integer(max_terms_per_block %||% config$max_terms_per_block %||% 5L)

  # Manual override via custom_blocks
  if (!is.null(config$custom_blocks) && length(config$custom_blocks) > 0L) {
    cb <- config$custom_blocks
    b_names <- names(cb) %||% paste0("block_", seq_along(cb))
    blocks <- lapply(seq_along(cb), function(i) {
      t_vec <- unique(tolower(stringr::str_squish(as.character(cb[[i]]))))
      t_vec <- t_vec[nzchar(t_vec)]
      list(name = b_names[i], terms = utils::head(t_vec, max_terms_per_block))
    })
    return(Filter(function(b) length(b$terms) > 0L, blocks))
  }

  k <- nrow(seeds_df)
  if (nrow(terms_df) == 0L) return(list())

  if (is.null(probe_selectivity)) {
    probe_selectivity <- isTRUE(attr(terms_df, "probed_live")) || !is.null(eval_joint_fn)
  }
  if (isTRUE(probe_selectivity) && is.null(eval_joint_fn)) {
    pm <- attr(terms_df, "probe_meta") %||% list()
    eval_joint_fn <- function(t1, t2) {
      res <- pubmed_esearch(
        term = sprintf('"%s"[tiab] AND "%s"[tiab]', t1, t2), retmax = 0L,
        email = pm$email, api_key = pm$ncbi_key,
        cache_dir = pm$cache_dir %||% ".cache", use_cache = pm$use_cache %||% TRUE
      )
      if (isTRUE(res$ok)) res$count else NA_integer_
    }
  }

  min_anchor_df <- if (k <= 3L) 1L else max(2L, ceiling(0.25 * k))
  blocks <- list()

  # --- Block 1: primary_topic ---
  b1_elig <- which(terms_df$df_tiab >= min_anchor_df & terms_df$lift >= 5)
  if (length(b1_elig) == 0L) b1_elig <- which(terms_df$df_tiab >= 1L)
  if (length(b1_elig) == 0L) return(list())

  a1_idx <- b1_elig[which.max(terms_df$term_score[b1_elig])]
  a1 <- terms_df$term[a1_idx]
  a1_bg <- as.numeric(terms_df$bg_hits[a1_idx])
  a1_cov <- terms_df$cov_tiab[[a1_idx]]

  # Probe joint hits P(a1, t) for top 20 stem-disjoint candidates to detect domain-coupled synonyms
  j_hits_vec <- stats::setNames(numeric(0), character(0))
  probe_cands <- terms_df$term[
    terms_df$term != a1 &
      !vapply(terms_df$term, shares_stem, logical(1), t2 = a1)
  ]
  probe_cands <- utils::head(probe_cands, 20L)

  if (isTRUE(probe_selectivity) && is.function(eval_joint_fn) && length(probe_cands) > 0L) {
    jh_vals <- vapply(probe_cands, function(tc) {
      jh <- suppressWarnings(as.numeric(eval_joint_fn(a1, tc)))
      if (length(jh) == 1L && !is.na(jh)) jh else NA_real_
    }, numeric(1))
    j_hits_vec <- jh_vals[!is.na(jh_vals)]
  }

  raw_families <- assign_term_families(terms_df$term)
  coupled_mask <- vapply(seq_len(nrow(terms_df)), function(idx) {
    tc <- terms_df$term[idx]
    if (tc == a1 || shares_stem(tc, a1) || .is_code_family(a1, tc)) return(TRUE)
    if (unname(raw_families[idx]) %in% c("outcome_severity", "population_host", "epi_methods")) {
      return(FALSE)
    }
    c_bg <- as.numeric(terms_df$bg_hits[idx])
    if (tc %in% names(j_hits_vec)) {
      jh <- j_hits_vec[[tc]]
      if ((jh / max(a1_bg, 1)) > 0.55) return(TRUE)
      if ((jh / max(c_bg, 1)) > 0.15) return(TRUE)
      return(FALSE)
    }
    c_cov <- terms_df$cov_tiab[[idx]]
    co_prop <- sum(a1_cov & c_cov) / max(sum(a1_cov), 1L)
    c_bg > 50000 && co_prop > 0.85 && terms_df$lift[idx] < 500
  }, logical(1))
  block1_coupled <- unique(terms_df$term[coupled_mask])

  b1_terms <- a1
  b1_cov <- a1_cov
  b1_max <- min(3L, max_terms_per_block)
  if (mean(b1_cov) < 0.98) {
    b1_syns <- which(
      terms_df$term != a1 &
        (terms_df$term %in% block1_coupled) &
        (terms_df$lift >= 10 | vapply(terms_df$term, .is_code_family, logical(1), t1 = a1)) &
        terms_df$bg_hits <= 1500000
    )
    while (length(b1_terms) < b1_max && mean(b1_cov) < 0.98 && length(b1_syns) > 0L) {
      rem_syns <- setdiff(b1_syns, match(b1_terms, terms_df$term))
      if (length(rem_syns) == 0L) break
      syn_scores <- vapply(rem_syns, function(si) {
        c_cov <- terms_df$cov_tiab[[si]]
        new_cov <- sum(c_cov & !b1_cov)
        if (new_cov == 0L) return(-1)
        tc <- terms_df$term[si]
        is_code_rel <- .is_code_family(a1, tc) || shares_stem(a1, tc)
        new_cov * log10(pmax(terms_df$lift[si], 1.1)) / log10(pmax(terms_df$bg_hits[si], 10)) *
          (if (is_code_rel) 5.0 else 1.0)
      }, numeric(1))
      if (max(syn_scores) <= 0) break
      best_si <- rem_syns[which.max(syn_scores)]
      b1_terms <- c(b1_terms, terms_df$term[best_si])
      b1_cov <- b1_cov | terms_df$cov_tiab[[best_si]]
    }
  }
  blocks[[1L]] <- list(name = "primary_topic", terms = b1_terms)

  # Identify universal domain context qualifier (e.g. "influenza") coupled to Block 1 with >=85% seed coverage
  dq_rows <- which(
    (terms_df$term %in% block1_coupled) &
      !(terms_df$term %in% b1_terms) &
      !vapply(terms_df$term, shares_stem, logical(1), t2 = a1) &
      !vapply(terms_df$term, .is_code_family, logical(1), t1 = a1) &
      terms_df$seed_prop >= 0.85 &
      terms_df$lift >= 40 &
      terms_df$n_words == 1L
  )
  domain_qualifier <- if (length(dq_rows) > 0L) {
    terms_df$term[dq_rows[which.max(terms_df$term_score[dq_rows])]]
  } else {
    character(0)
  }

  if (max_blocks <= 1L) {
    attr(blocks, "block1_coupled") <- block1_coupled
    attr(blocks, "joint_hits") <- j_hits_vec
    attr(blocks, "term_families") <- assign_term_families(terms_df$term, block1_coupled)
    attr(blocks, "domain_qualifier") <- domain_qualifier
    return(blocks)
  }

  # --- Blocks 2..max_blocks: semantically pure PICO families + igraph fallback for "other" ---
  term_families <- assign_term_families(terms_df$term, block1_coupled)

  # Helper: check conditional selectivity for candidate idx relative to a1
  is_selective_vs_a1 <- function(idx, max_bg = 1500000) {
    tc <- terms_df$term[idx]
    c_bg <- as.numeric(terms_df$bg_hits[idx])
    is_pico <- unname(raw_families[idx]) %in% c("outcome_severity", "population_host", "epi_methods")
    if (tc %in% names(j_hits_vec)) {
      jh <- j_hits_vec[[tc]]
      if (is_pico) return((jh / max(a1_bg, 1)) <= 0.65 && c_bg <= max_bg)
      return((jh / max(a1_bg, 1)) <= 0.55 && (jh / max(c_bg, 1)) <= 0.15)
    }
    if (!is_pico) {
      c_cov <- terms_df$cov_tiab[[idx]]
      co_prop <- sum(a1_cov & c_cov) / max(sum(a1_cov), 1L)
      if (c_bg > 50000 && co_prop > 0.85 && terms_df$lift[idx] < 500) return(FALSE)
    }
    c_bg <= max_bg
  }

  cand_blocks <- list()
  used_terms <- b1_terms

  # 1. Build pure blocks for universal PICO families
  for (fam in c("outcome_severity", "population_host", "epi_methods")) {
    max_fam_bg <- if (fam == "population_host") 4000000 else 1500000
    fam_idx <- which(
      term_families == fam &
        !(terms_df$term %in% used_terms) &
        !(terms_df$term %in% block1_coupled) &
        terms_df$df_tiab >= min_anchor_df &
        terms_df$lift >= 5
    )
    fam_idx <- Filter(function(idx) {
      !any(vapply(b1_terms, shares_stem, logical(1), t2 = terms_df$term[idx])) &&
        is_selective_vs_a1(idx, max_bg = max_fam_bg)
    }, fam_idx)
    if (length(fam_idx) == 0L) next

    # Score family anchor by seed coverage * conditional selectivity vs a1
    anchor_scores <- vapply(fam_idx, function(idx) {
      tc <- terms_df$term[idx]
      jh <- if (tc %in% names(j_hits_vec)) j_hits_vec[[tc]] else NA_real_
      red <- if (!is.na(jh)) log10(pmax(a1_bg / pmax(jh, 1), 1.1)) else log10(pmax(38000000 / pmax(terms_df$bg_hits[idx], 100), 1.1))
      (terms_df$df_tiab[idx] / k) * red *
        (1 + 0.15 * (terms_df$n_words[idx] > 1L)) *
        (1 + 0.25 * endsWith(tc, "*"))
    }, numeric(1))

    a_idx <- fam_idx[which.max(anchor_scores)]
    a_term <- terms_df$term[a_idx]
    a_bg <- terms_df$bg_hits[a_idx]
    a_df <- terms_df$df_tiab[a_idx]
    fam_terms <- a_term
    fam_cov <- terms_df$cov_tiab[[a_idx]]

    # Grow purely within the same semantic family by greedy marginal uncovered seed coverage
    pool_idx <- which(
      term_families == fam &
        terms_df$term != a_term &
        !(terms_df$term %in% used_terms) &
        !(terms_df$term %in% block1_coupled) &
        terms_df$lift >= 5
    )
    pool_idx <- Filter(function(idx) {
      tc <- terms_df$term[idx]
      c_bg <- terms_df$bg_hits[idx]
      if (any(vapply(b1_terms, shares_stem, logical(1), t2 = tc))) return(FALSE)
      if (fam != "population_host" && c_bg > 1500000) return(FALSE)
      if (tc %in% names(j_hits_vec) && fam != "population_host" && j_hits_vec[[tc]] > 2500) return(FALSE)
      is_selective_vs_a1(idx, max_bg = max_fam_bg)
    }, pool_idx)

    while (length(fam_terms) < max_terms_per_block && length(pool_idx) > 0L) {
      rem_pool <- Filter(function(pi) {
        tc <- terms_df$term[pi]
        !(tc %in% fam_terms) && !any(vapply(fam_terms, .term_subsumes, logical(1), t_narrow = tc))
      }, pool_idx)
      if (length(rem_pool) == 0L) break
      gains <- vapply(rem_pool, function(pi) {
        c_cov <- terms_df$cov_tiab[[pi]]
        new_cov <- sum(c_cov & !fam_cov)
        if (new_cov == 0L) return(-1)
        new_cov * log10(pmax(terms_df$lift[pi], 1.1)) / log10(pmax(terms_df$bg_hits[pi], 100))
      }, numeric(1))
      if (max(gains) <= 0) break
      best_pi <- rem_pool[which.max(gains)]
      fam_terms <- c(fam_terms, terms_df$term[best_pi])
      fam_cov <- fam_cov | terms_df$cov_tiab[[best_pi]]
    }

    # If no complementary terms had new_cov > 0 (e.g. anchor already covered 100%), add at most 1 narrower high-precision term
    if (length(fam_terms) == 1L && length(pool_idx) > 0L) {
      narrow_cands <- pool_idx[terms_df$bg_hits[pool_idx] < a_bg & terms_df$df_tiab[pool_idx] >= ceiling(0.4 * a_df)]
      if (length(narrow_cands) > 0L) {
        best_n <- narrow_cands[which.max(terms_df$term_score[narrow_cands])]
        fam_terms <- c(fam_terms, terms_df$term[best_n])
      }
    }

    # Also record narrowest high-coverage term in this family for ultra-specific 1-term tiers
    narrow_pool <- fam_idx[terms_df$bg_hits[fam_idx] < a_bg & terms_df$df_tiab[fam_idx] >= ceiling(0.45 * a_df)]
    narrow_term <- if (length(narrow_pool) > 0L) {
      terms_df$term[narrow_pool[which.min(terms_df$bg_hits[narrow_pool])]]
    } else {
      character(0)
    }

    used_terms <- c(used_terms, fam_terms)

    best_jh <- min(c(a1_bg, j_hits_vec[intersect(fam_terms, names(j_hits_vec))]), na.rm = TRUE)
    blk_score <- mean(fam_cov) * log10(pmax(a1_bg / pmax(best_jh, 10), 1.5))

    cand_blocks[[length(cand_blocks) + 1L]] <- list(
      name = fam,
      terms = fam_terms,
      narrow_term = narrow_term,
      score = blk_score
    )
  }

  # 2. If fewer than (max_blocks - 1) PICO blocks exist, cluster remaining "other" terms via igraph
  if (length(cand_blocks) < (max_blocks - 1L)) {
    other_idx <- which(
      term_families == "other" &
        !(terms_df$term %in% used_terms) &
        !(terms_df$term %in% block1_coupled) &
        terms_df$df_tiab >= min_anchor_df &
        terms_df$lift >= 10 &
        terms_df$bg_hits <= 1500000
    )
    other_idx <- Filter(function(idx) {
      !any(vapply(b1_terms, shares_stem, logical(1), t2 = terms_df$term[idx])) &&
        is_selective_vs_a1(idx, max_bg = 1500000)
    }, other_idx)

    if (length(other_idx) >= 1L) {
      if (length(other_idx) == 1L) {
        oi <- other_idx[1L]
        cand_blocks[[length(cand_blocks) + 1L]] <- list(
          name = "concept_other_1",
          terms = terms_df$term[oi],
          score = terms_df$seed_prop[oi] * 0.5
        )
      } else {
        m_o <- length(other_idx)
        edges <- integer(0)
        for (ii in seq_len(m_o - 1L)) {
          for (jj in (ii + 1L):m_o) {
            ti <- terms_df$term[other_idx[ii]]
            tj <- terms_df$term[other_idx[jj]]
            ci <- terms_df$cov_tiab[[other_idx[ii]]]
            cj <- terms_df$cov_tiab[[other_idx[jj]]]
            jacc <- sum(ci & cj) / max(sum(ci | cj), 1L)
            if (shares_stem(ti, tj) || jacc >= 0.45) {
              edges <- c(edges, ii, jj)
            }
          }
        }
        g <- igraph::make_empty_graph(n = m_o, directed = FALSE)
        if (length(edges) > 0L) g <- igraph::add_edges(g, edges)
        mem <- if (igraph::ecount(g) > 0L) igraph::membership(igraph::cluster_walktrap(g)) else seq_len(m_o)
        for (cid in unique(mem)) {
          grp_rows <- other_idx[which(mem == cid)]
          grp_rows <- grp_rows[order(terms_df$term_score[grp_rows], decreasing = TRUE)]
          grp_rows <- utils::head(grp_rows, max_terms_per_block)
          grp_cov <- Reduce(`|`, terms_df$cov_tiab[grp_rows])
          cand_blocks[[length(cand_blocks) + 1L]] <- list(
            name = paste0("concept_other_", cid),
            terms = terms_df$term[grp_rows],
            score = mean(grp_cov) * 0.6
          )
        }
      }
    }
  }

  # Keep top (max_blocks - 1) orthogonal blocks by score
  if (length(cand_blocks) > 0L) {
    ord_cb <- order(vapply(cand_blocks, `[[`, numeric(1), "score"), decreasing = TRUE)
    cand_blocks <- utils::head(cand_blocks[ord_cb], max_blocks - 1L)
    for (cb in cand_blocks) {
      blk_item <- list(name = cb$name, terms = cb$terms)
      if (!is.null(cb$narrow_term) && length(cb$narrow_term) > 0L) {
        attr(blk_item, "narrow_term") <- cb$narrow_term
      }
      blocks[[length(blocks) + 1L]] <- blk_item
    }
  }

  attr(blocks, "block1_coupled") <- block1_coupled
  attr(blocks, "joint_hits") <- j_hits_vec
  attr(blocks, "term_families") <- term_families
  attr(blocks, "domain_qualifier") <- domain_qualifier
  blocks
}

# R/04_term_extraction.R
# Unsupervised n-gram & prefix-wildcard extraction, PubMed background lift scoring,
# distributional SVD + MeSH/keyword + complementarity graph clustering, and
# orthogonal set-cover concept block discovery.

#' Stopword list (standard English + generic academic reporting boilerplate)
#'
#' @return Character vector of stopwords.
#' @export
get_stopwords <- function() {
  unique(c(
    # Standard English function words
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
    # Generic cross-disciplinary academic reporting words
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
    "least", "much", "about", "after", "before", "above", "below", "against",
    "along", "alongside", "without", "throughout", "toward", "towards", "beyond", "despite"
  ))
}

#' Normalize text for n-gram extraction (preserves alphanumeric tokens)
#'
#' @param x Character vector of texts.
#' @return Normalized lowercase alphanumeric text vector.
#' @export
normalize_search_text <- function(x) {
  x <- tolower(as.character(x %||% ""))
  x[is.na(x)] <- ""
  x <- stringr::str_replace_all(x, "[-/]+", " ")
  x <- stringr::str_replace_all(x, "[^a-z0-9\\s]", " ")
  stringr::str_squish(x)
}

#' Escape regex metacharacters in literal term
#'
#' @param x Character string to escape.
#' @return Regex-escaped character string.
#' @export
escape_regex <- function(x) {
  stringr::str_replace_all(x, "([.\\^$|()*+?{}\\[\\]\\\\])", "\\\\\\1")
}

#' Check whether a candidate term matches each text in `texts` (case-insensitive, hyphen-flexible)
#'
#' @param term Single candidate term (supports trailing `*` prefix wildcard).
#' @param texts Character vector of texts.
#' @return Logical vector of length `length(texts)`.
#' @export
term_matches_text <- function(term, texts) {
  if (!nzchar(term) || length(texts) == 0L) {
    return(rep(FALSE, length(texts)))
  }
  t_norm <- stringr::str_squish(tolower(stringr::str_replace_all(term, "[-/]+", " ")))
  if (!nzchar(t_norm)) {
    return(rep(FALSE, length(texts)))
  }
  pat <- escape_regex(t_norm)
  pat <- stringr::str_replace_all(pat, "\\\\\\*", "[a-z0-9]*")
  pat <- stringr::str_replace_all(pat, " ", "[\\\\s-]+")
  full_pat <- paste0("(?i)(?<![a-z0-9])", pat, "(?![a-z0-9])")
  texts_norm <- stringr::str_squish(tolower(as.character(texts %||% "")))
  grepl(full_pat, texts_norm, perl = TRUE)
}

#' Simple token stemmer for lexical overlap detection across concept blocks
#'
#' @param term Candidate term string.
#' @return Character vector of stemmed tokens.
#' @export
term_stems <- function(term) {
  stops <- get_stopwords()
  w <- unlist(stringr::str_split(normalize_search_text(gsub("\\*", "", term)), " "))
  w <- w[nzchar(w) & !(w %in% stops) & nchar(w) >= 2L]
  if (length(w) == 0L) {
    return(character(0))
  }
  stems <- vapply(w, function(tok) {
    if (grepl("[0-9]", tok) || nchar(tok) <= 4L) {
      return(tok)
    }
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

#' Check whether two terms share an alphanumeric code prefix (e.g. brca1/brca, h5n1/h5)
#'
#' @param t1 First candidate term.
#' @param t2 Second candidate term.
#' @return Logical scalar.
#' @keywords internal
.is_code_family <- function(t1, t2) {
  t1 <- tolower(stringr::str_squish(as.character(t1 %||% "")))
  t2 <- tolower(stringr::str_squish(as.character(t2 %||% "")))
  if (!grepl("^[a-z]+[0-9]+", t1) || !grepl("^[a-z]+[0-9]+", t2)) {
    return(FALSE)
  }
  p1 <- regmatches(t1, regexpr("^[a-z]+[0-9]+", t1))
  p2 <- regmatches(t2, regexpr("^[a-z]+[0-9]+", t2))
  identical(p1, p2) && (startsWith(t1, t2) || startsWith(t2, t1))
}

#' Check whether term `t_broad` lexically subsumes `t_narrow` inside an OR block
#'
#' @param t_broad Broader candidate term.
#' @param t_narrow Narrower candidate term.
#' @return Logical scalar.
#' @keywords internal
.term_subsumes <- function(t_broad, t_narrow) {
  t_broad <- tolower(stringr::str_squish(as.character(t_broad %||% "")))
  t_narrow <- tolower(stringr::str_squish(as.character(t_narrow %||% "")))
  if (!nzchar(t_broad) || !nzchar(t_narrow) || identical(t_broad, t_narrow)) {
    return(FALSE)
  }
  pat <- escape_regex(stringr::str_squish(stringr::str_replace_all(t_broad, "[-/]+", " ")))
  pat <- stringr::str_replace_all(pat, "\\\\\\*", "[a-z0-9]*")
  pat <- stringr::str_replace_all(pat, " ", "[\\\\s-]+")
  t_narrow_norm <- stringr::str_squish(stringr::str_replace_all(t_narrow, "[-/]+", " "))
  grepl(paste0("(?<![a-z0-9])", pat, "(?![a-z0-9])"), t_narrow_norm, perl = TRUE)
}

#' Check whether two terms share any significant token stem
#'
#' @param t1 First candidate term.
#' @param t2 Second candidate term.
#' @return Logical scalar.
#' @export
shares_stem <- function(t1, t2) {
  s1 <- term_stems(t1)
  s2 <- term_stems(t2)
  length(intersect(s1, s2)) > 0L
}

# Extract 1..3-grams from a single normalized text
.extract_doc_ngrams <- function(text, max_n = 3L, stops = get_stopwords()) {
  if (!nzchar(text)) {
    return(character(0))
  }
  toks <- unlist(stringr::str_split(text, " "))
  toks <- toks[nzchar(toks)]
  nt <- length(toks)
  if (nt == 0L) {
    return(character(0))
  }

  grams <- character(0)
  for (n_g in seq_len(min(max_n, nt))) {
    if (n_g == 1L) {
      keep <- !(toks %in% stops) & nchar(toks) >= 2L & !grepl("^[0-9]+$", toks)
      grams <- c(grams, toks[keep])
    } else {
      idx <- seq_len(nt - n_g + 1L)
      first_w <- toks[idx]
      last_w <- toks[idx + n_g - 1L]
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

# Generic cross-disciplinary reporting stems that should not form prefix wildcards
.is_generic_academic_prefix <- function(p) {
  grepl(
    "^(studi|result|method|analys|conclu|observ|compar|includ|report|associ|signif|differ|increa|decrea|presen|perfor|evalua|determ|invest|measur|estima|calcul|consid|sugges|indica|demons|develo|identi|provid|requir|obtain|potent|number|period|sample|effect|factor|change|level|group|value|relat|approa|descri|examin|assess|accord|follow|addition|genera|specif|import|possib|recent|previous|simila|variou|underst|attribut|deriv|repres|corres|approx|respec|subseq|suffic|charac|remain|reveal|confirm|explai|theor|broad|diver|robus|notab|integra|distinc|essenti|critic|multip|individ|comple|implem|facilit|intensi|contrib|establ|adminis|design|combin|responsi|commerci|functio|contex|seriou|detect|optimi|moderat)",
    p,
    perl = TRUE
  )
}

# Morphological guard preventing accidental Latin/Greek prefix collisions (e.g. inter*, trans*, speci*, count*)
.is_valid_morphological_prefix <- function(p, m_u) {
  if (length(m_u) < 2L) {
    return(FALSE)
  }
  sfx <- substring(m_u, nchar(p) + 1L)
  if (any(nchar(sfx) > 5L)) {
    return(FALSE)
  }
  if (nchar(p) <= 6L) {
    has_base <- any(sfx %in% c("", "e", "s", "ed", "y", "al", "ic", "or", "er"))
    if (!has_base) {
      return(FALSE)
    }
    valid_inflect <- grepl(
      "^(|e|es|s|ed|d|ing|er|ers|or|ors|al|ally|ity|ities|y|ies|ly|ely|ic|ics|ical|ist|ists|ism|ive|ous|ious|ion|ions|ure|ures|ation|ations|ated|ating|ment|ments|ness|able)$",
      sfx
    )
    if (!all(valid_inflect)) {
      return(FALSE)
    }
  }
  TRUE
}

# Unsupervised morphological prefix wildcard synthesis across any domain vocabulary:
# finds 5..12-char prefixes matching >= 2 distinct corpus unigrams that rescue >= 2 additional seed documents
.synthesize_prefix_wildcards <- function(doc_grams, raw_texts, min_df = 2L, stops = get_stopwords()) {
  k <- length(raw_texts)
  if (k < 6L) {
    return(character(0))
  }
  all_u <- unique(unlist(doc_grams))
  u <- all_u[grepl("^[a-z]{4,}$", all_u) & !(all_u %in% stops) & !.is_generic_academic_prefix(all_u)]
  if (length(u) < 2L) {
    return(character(0))
  }

  u_counts <- table(unlist(lapply(doc_grams, function(g) intersect(g, u))))
  cand_prefixes <- unique(unlist(lapply(5:12, function(len) {
    long_u <- u[nchar(u) >= len]
    if (length(long_u) == 0L) character(0) else substr(long_u, 1L, len)
  })))
  cand_prefixes <- cand_prefixes[!(cand_prefixes %in% stops) & !.is_generic_academic_prefix(cand_prefixes)]
  if (length(cand_prefixes) == 0L) {
    return(character(0))
  }

  matched_sets <- lapply(cand_prefixes, function(p) sort(u[startsWith(u, p)]))
  keep_p <- vapply(seq_along(cand_prefixes), function(i) {
    .is_valid_morphological_prefix(cand_prefixes[i], matched_sets[[i]])
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
      if (exists(".register_wildcard_expansion", mode = "function")) {
        m_cnts <- as.integer(u_counts[m_u])
        best_w <- m_u[order(-m_cnts, nchar(m_u))][1L]
        .register_wildcard_expansion(p, best_w)
      }
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

# Multi-layer stem-disjoint greedy set-cover shortlisting before PubMed probing:
# guarantees the shortlist spans multiple orthogonal conceptual facets AND includes
# low-DF complementary synonyms needed to rescue outlier seed papers.
.multilayer_setcover_shortlist <- function(keep_pool, cov_tiab_list, local_score, df_tiab,
                                           is_acro, df_title, n_layers = 6L, max_per_layer = 8L) {
  n_terms <- length(keep_pool)
  if (n_terms == 0L) {
    return(integer(0))
  }
  k <- length(cov_tiab_list[[1L]])
  stems_list <- lapply(keep_pool, term_stems)
  max_gap_df <- max(2L, ceiling(0.60 * k))

  blocked_stems <- character(0)
  picked_all <- integer(0)

  for (layer in seq_len(n_layers)) {
    elig <- which(vapply(seq_len(n_terms), function(i) {
      if (i %in% picked_all) {
        return(FALSE)
      }
      st <- stems_list[[i]]
      if (length(st) == 0L) {
        return(FALSE)
      }
      !any(st %in% blocked_stems)
    }, logical(1)))
    if (length(elig) == 0L) break

    # Start layer from highest local_score term in `elig`
    a_idx <- elig[which.max(local_score[elig])]
    a_term <- keep_pool[a_idx]
    layer_picks <- a_idx
    layer_cov <- cov_tiab_list[[a_idx]]

    # Greedily add selective gap-filling terms (preferring non-ubiquitous terms so pillars don't mask synonyms)
    while (length(layer_picks) < max_per_layer && !all(layer_cov)) {
      rem_all <- setdiff(elig, layer_picks)
      if (length(rem_all) == 0L) break
      uncov <- !layer_cov
      rem_pref <- rem_all[
        df_tiab[rem_all] <= max_gap_df |
          vapply(keep_pool[rem_all], .is_code_family, logical(1), t1 = a_term)
      ]
      rem <- if (length(rem_pref) > 0L) rem_pref else rem_all
      score_fn <- function(cand_vec) {
        vapply(cand_vec, function(idx) {
          new_hits <- sum(cov_tiab_list[[idx]][uncov])
          if (new_hits == 0L) {
            return(-1)
          }
          code_bonus <- if (.is_code_family(a_term, keep_pool[idx])) 2.5 else 1.0
          new_hits * (1 + 0.5 * local_score[idx]) * (1 + 0.8 * is_acro[idx]) *
            (1 + 0.3 * (df_title[idx] > 0L)) * code_bonus
        }, numeric(1))
      }
      gains <- score_fn(rem)
      if (max(gains) <= 0 && length(rem) < length(rem_all)) {
        rem <- rem_all
        gains <- score_fn(rem)
      }
      if (max(gains) <= 0) break
      best_r <- rem[which.max(gains)]
      layer_picks <- c(layer_picks, best_r)
      layer_cov <- layer_cov | cov_tiab_list[[best_r]]
    }

    picked_all <- c(picked_all, layer_picks)
    # Block the anchor stem(s) and high-frequency stems of this layer so the next layer discovers an orthogonal facet
    layer_high <- layer_picks[df_tiab[layer_picks] >= max(2L, ceiling(0.30 * k))]
    if (length(layer_high) == 0L) layer_high <- layer_picks[1L]
    blocked_stems <- unique(c(blocked_stems, unlist(stems_list[layer_high])))
  }

  unique(picked_all)
}

#' Extract and score candidate search terms from enriched seed papers
#'
#' @param seeds_df Enriched tibble from `enrich_seeds()`.
#' @param config Optional config list (`banned_terms`, `required_terms`, `custom_blocks`).
#' @param max_candidates Max terms to probe against PubMed background frequency.
#' @param probe_pubmed Logical; whether to probe PubMed for background `bg_hits`.
#' @param eval_bg_fn Optional function `(term) -> integer` for testing without network.
#' @param email Optional email for NCBI polite pool.
#' @param ncbi_key Optional NCBI API key.
#' @param cache_dir Cache directory path.
#' @param use_cache Logical; whether to read/write cache.
#' @param progress_cb Optional progress callback function.
#' @return Tibble of scored candidate terms with `cov_tiab` and `cov_title` list-columns.
#' @export
extract_candidate_terms <- function(seeds_df,
                                    config = list(),
                                    max_candidates = 50L,
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

  doc_grams <- lapply(norm_texts, .extract_doc_ngrams, max_n = 3L, stops = stops)

  extra_terms <- unique(c(
    unlist(seeds_df$all_keywords %||% list()),
    unlist(seeds_df$mesh_terms %||% list())
  ))
  extra_terms <- normalize_search_text(extra_terms)
  extra_terms <- extra_terms[nzchar(extra_terms) & !(extra_terms %in% stops) & nchar(extra_terms) >= 2L]

  gram_counts <- table(unlist(doc_grams))
  min_df <- if (k <= 4L) 1L else 2L

  frequent_grams <- names(gram_counts[gram_counts >= min_df])
  prefix_wildcards <- .synthesize_prefix_wildcards(doc_grams, raw_texts, min_df = min_df, stops = stops)
  per_seed_rescue <- unique(unlist(lapply(seq_len(k), function(i) {
    t_grams <- .extract_doc_ngrams(normalize_search_text(raw_titles[i]), max_n = 3L, stops = stops)
    kw_i <- normalize_search_text(c(seeds_df$all_keywords[[i]] %||% character(0), seeds_df$mesh_terms[[i]] %||% character(0)))
    kw_i <- kw_i[nzchar(kw_i) & !(kw_i %in% stops)]
    dg_i <- doc_grams[[i]]
    u_i <- dg_i[!grepl(" ", dg_i, fixed = TRUE) & nchar(dg_i) >= 4L & !.is_generic_academic_prefix(dg_i)]
    if (length(u_i) > 0L) {
      u_i <- u_i[order(as.integer(gram_counts[u_i]), -nchar(u_i))]
    }
    c(utils::head(t_grams, 8L), utils::head(kw_i, 4L), utils::head(u_i, 10L))
  })))

  keep_pool <- unique(c(frequent_grams, prefix_wildcards, per_seed_rescue, extra_terms, pinned_terms))
  keep_pool <- keep_pool[nzchar(keep_pool) & !(keep_pool %in% banned) & !(keep_pool %in% stops)]
  keep_pool <- keep_pool[nchar(keep_pool) >= 2L & !grepl("^[0-9\\s]+$", keep_pool)]

  cov_tiab_list <- lapply(keep_pool, term_matches_text, texts = raw_texts)
  cov_title_list <- lapply(keep_pool, term_matches_text, texts = raw_titles)

  df_tiab <- vapply(cov_tiab_list, sum, integer(1))
  df_title <- vapply(cov_title_list, sum, integer(1))

  keep_idx <- (df_tiab >= 1L) | (keep_pool %in% pinned_terms)
  keep_pool <- keep_pool[keep_idx]
  cov_tiab_list <- cov_tiab_list[keep_idx]
  cov_title_list <- cov_title_list[keep_idx]
  df_tiab <- df_tiab[keep_idx]
  df_title <- df_title[keep_idx]

  n_words <- vapply(strsplit(keep_pool, " "), length, integer(1))
  is_acro <- grepl("[0-9]", keep_pool) | (keep_pool %in% acronyms)
  is_wildcard <- endsWith(keep_pool, "*")

  local_score <- (df_tiab / max(k, 1L)) *
    (1 + 0.5 * (df_title / pmax(df_tiab, 1L))) *
    (1 + 0.20 * (n_words > 1L)) *
    (1 + 0.25 * is_acro) *
    (1 + 0.25 * is_wildcard)

  # Retain top 3 terms per seed paper
  rescue_idx <- unique(unlist(lapply(seq_len(k), function(i) {
    covers_i <- which(vapply(cov_tiab_list, `[`, logical(1), i))
    if (length(covers_i) == 0L) {
      return(integer(0))
    }
    ord <- covers_i[order(local_score[covers_i], decreasing = TRUE)]
    utils::head(ord, 3L)
  })))

  # Multi-layer orthogonal set-cover shortlist across disjoint stems
  layer_idx <- .multilayer_setcover_shortlist(
    keep_pool = keep_pool,
    cov_tiab_list = cov_tiab_list,
    local_score = local_score,
    df_tiab = df_tiab,
    is_acro = is_acro,
    df_title = df_title,
    n_layers = 6L,
    max_per_layer = 8L
  )

  # Explicitly shortlist selective rescue candidates for outlier seeds missed by >=75%-coverage pillars
  high_pillars <- which(df_tiab >= max(2L, ceiling(0.75 * k)))
  outlier_rescue_idx <- integer(0)
  if (length(high_pillars) > 0L && k >= 5L) {
    hp_ord <- high_pillars[order(local_score[high_pillars], decreasing = TRUE)]
    hp_keep <- integer(0)
    hp_stems <- character(0)
    for (hi in hp_ord) {
      st_h <- term_stems(keep_pool[hi])
      if (length(st_h) > 0L && !any(st_h %in% hp_stems)) {
        hp_keep <- c(hp_keep, hi)
        hp_stems <- c(hp_stems, st_h)
      }
      if (length(hp_keep) >= 4L) break
    }
    outlier_seeds <- unique(unlist(lapply(hp_keep, function(hi) which(!cov_tiab_list[[hi]]))))
    max_out_df <- max(2L, ceiling(0.35 * k))
    cand_has_hp_stem <- vapply(keep_pool, function(t) any(term_stems(t) %in% hp_stems), logical(1))
    outlier_rescue_idx <- unique(unlist(lapply(outlier_seeds, function(si) {
      cand_s <- which(
        vapply(cov_tiab_list, `[`, logical(1), si) &
          n_words <= 2L &
          df_tiab <= max_out_df &
          !cand_has_hp_stem &
          !.is_generic_academic_prefix(keep_pool)
      )
      if (length(cand_s) == 0L) {
        return(integer(0))
      }
      sc_s <- (1 + 1.5 * is_acro[cand_s]) * (1 + 0.6 * (n_words[cand_s] == 1L)) *
        (1 + 0.4 * (df_title[cand_s] > 0L)) * (1 + 0.15 * (nchar(keep_pool[cand_s]) >= 5L))
      utils::head(cand_s[order(sc_s, decreasing = TRUE)], 5L)
    })))
  }

  # Enforce stem-diversity cap across ALL stems in multi-word phrases so permutations of Block 1 don't crowd out other facets
  ord_local <- order(local_score, decreasing = TRUE)
  diverse_idx <- integer(0)
  stem_counts <- integer(0)
  for (idx in ord_local) {
    if (length(diverse_idx) >= max_candidates) break
    st <- term_stems(keep_pool[idx])
    if (length(st) == 0L) next
    cnts <- vapply(st, function(s) if (s %in% names(stem_counts)) stem_counts[[s]] else 0L, integer(1))
    # Accept if at least one stem in the term has been seen < 3 times
    if (min(cnts) < 3L) {
      diverse_idx <- c(diverse_idx, idx)
      for (s in st) {
        stem_counts[[s]] <- (if (s %in% names(stem_counts)) stem_counts[[s]] else 0L) + 1L
      }
    }
  }

  pin_idx <- which(keep_pool %in% pinned_terms)
  sel_idx <- unique(c(pin_idx, layer_idx, outlier_rescue_idx, diverse_idx, rescue_idx))
  sel_idx <- sel_idx[seq_len(min(length(sel_idx), max_candidates + 22L))]

  terms_sel <- keep_pool[sel_idx]
  cov_tiab_sel <- cov_tiab_list[sel_idx]
  cov_title_sel <- cov_title_list[sel_idx]
  df_tiab_sel <- df_tiab[sel_idx]
  df_title_sel <- df_title[sel_idx]
  n_words_sel <- n_words[sel_idx]
  is_acro_sel <- is_acro[sel_idx]

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
  bg_missing <- is.na(bg_hits)
  probed_live <- probed_live && any(!bg_missing)
  if (any(bg_missing)) {
    known_idx <- which(!bg_missing)
    for (mi in which(bg_missing)) {
      sub_known <- if (length(known_idx) > 0L) {
        known_idx[vapply(terms_sel[known_idx], function(kt) .term_subsumes(terms_sel[mi], kt), logical(1))]
      } else {
        integer(0)
      }
      if (length(sub_known) > 0L) {
        bg_hits[mi] <- as.integer(min(30000000, ceiling(max(bg_hits[sub_known]) * 1.15)))
      } else {
        bg_hits[mi] <- if (n_words_sel[mi] >= 3L) 2000L else if (n_words_sel[mi] == 2L) 25000L else 200000L
      }
    }
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

# Build second-order distributional context + MeSH/keyword SVD embeddings for candidate terms
.compute_term_context_similarity <- function(cand_terms, cov_list, seeds_df) {
  m <- length(cand_terms)
  if (m <= 1L) {
    return(matrix(1, nrow = m, ncol = m))
  }

  stops <- get_stopwords()
  raw_texts <- as.character(seeds_df$full_text %||% paste(seeds_df$title, seeds_df$abstract))
  norm_docs <- normalize_search_text(raw_texts)
  doc_tokens <- lapply(norm_docs, function(txt) {
    w <- unlist(stringr::str_split(txt, " "))
    w[nzchar(w) & !(w %in% stops) & nchar(w) >= 3L]
  })

  # Collect local +-5 word window context features + MeSH/keyword features for each candidate term
  ctx_bags <- vector("list", m)
  for (i in seq_len(m)) {
    t_raw <- gsub("\\*$", "", cand_terms[i])
    t_toks <- unlist(stringr::str_split(normalize_search_text(t_raw), " "))
    t_first <- t_toks[1L]
    is_pref <- endsWith(cand_terms[i], "*")
    feats <- character(0)

    hit_docs <- which(cov_list[[i]])
    for (d in hit_docs) {
      w <- doc_tokens[[d]]
      if (length(w) > 0L && nzchar(t_first)) {
        pos <- if (is_pref) which(startsWith(w, t_first)) else which(w == t_first)
        for (p in pos) {
          lo <- max(1L, p - 5L)
          hi <- min(length(w), p + length(t_toks) + 4L)
          win <- w[setdiff(lo:hi, p:min(length(w), p + length(t_toks) - 1L))]
          feats <- c(feats, paste0("ctx:", win))
        }
      }
      # Add bipartite MeSH descriptor & author keyword tags from matching documents
      mh <- seeds_df$mesh_terms[[d]] %||% character(0)
      kw <- seeds_df$all_keywords[[d]] %||% character(0)
      if (length(mh) > 0L) feats <- c(feats, paste0("mesh:", mh))
      if (length(kw) > 0L) feats <- c(feats, paste0("kw:", kw))
    }
    ctx_bags[[i]] <- feats
  }

  all_feats <- table(unlist(ctx_bags))
  # Keep informative features shared by at least 2 candidate terms (but < 95% of terms if m >= 6)
  vocab <- names(all_feats[all_feats >= 2L])
  if (length(vocab) == 0L) {
    return(diag(1, m))
  }
  if (length(vocab) > 300L) {
    vocab <- names(sort(all_feats[vocab], decreasing = TRUE))[1:300]
  }

  # Term x Context TF-IDF matrix
  W <- matrix(0, nrow = m, ncol = length(vocab))
  colnames(W) <- vocab
  for (i in seq_len(m)) {
    tab_i <- table(ctx_bags[[i]])
    common_v <- intersect(names(tab_i), vocab)
    if (length(common_v) > 0L) {
      W[i, common_v] <- log1p(as.numeric(tab_i[common_v]))
    }
  }
  df_v <- colSums(W > 0)
  idf_v <- log((m + 1) / (df_v + 1)) + 1
  W <- sweep(W, 2L, idf_v, `*`)

  # Center columns and apply truncated SVD (Latent Semantic Analysis) when dimensions permit
  if (m >= 4L && ncol(W) >= 4L) {
    W_centered <- scale(W, center = TRUE, scale = FALSE)
    rnk <- min(8L, nrow(W_centered) - 1L, ncol(W_centered) - 1L)
    if (rnk >= 2L) {
      sv <- tryCatch(svd(W_centered, nu = rnk, nv = 0L), error = function(e) NULL)
      if (!is.null(sv) && any(sv$d[seq_len(rnk)] > 1e-8)) {
        emb <- sweep(sv$u[, seq_len(rnk), drop = FALSE], 2L, sv$d[seq_len(rnk)], `*`)
        norms <- sqrt(rowSums(emb^2))
        norms[norms < 1e-12] <- 1
        emb_norm <- emb / norms
        sim <- tcrossprod(emb_norm)
        sim[sim < 0] <- 0
        diag(sim) <- 1
        return(sim)
      }
    }
  }

  norms <- sqrt(rowSums(W^2))
  norms[norms < 1e-12] <- 1
  W_norm <- W / norms
  sim <- tcrossprod(W_norm)
  sim[sim < 0] <- 0
  diag(sim) <- 1
  sim
}

# Character 4-gram Jaccard similarity for morphological / affix similarity
.char_ngram_sim <- function(t1, t2, n = 4L) {
  s1 <- gsub("[^a-z0-9]", "", tolower(t1))
  s2 <- gsub("[^a-z0-9]", "", tolower(t2))
  if (nchar(s1) < n || nchar(s2) < n) {
    return(0)
  }
  g1 <- unique(substring(s1, 1:(nchar(s1) - n + 1L), n:nchar(s1)))
  g2 <- unique(substring(s2, 1:(nchar(s2) - n + 1L), n:nchar(s2)))
  length(intersect(g1, g2)) / max(length(union(g1, g2)), 1L)
}

#' Unsupervised concept cluster discovery across candidate search terms
#'
#' Partitions candidate terms into coherent conceptual facets using a hybrid
#' affinity graph combining:
#' 1. Second-order distributional context & MeSH/keyword SVD similarity (`LSA`),
#' 2. Seed coverage complementarity (linking substitutable `OR` synonyms that cover
#'    each other's missed seed papers while keeping high-frequency orthogonal `AND`
#'    pillars separate), and
#' 3. Lexical stem, code-family, and character n-gram similarity,
#' followed by Louvain community detection (`igraph::cluster_louvain()`).
#'
#' @param terms_df Scored candidate terms tibble from [extract_candidate_terms()].
#' @param seeds_df Enriched seeds tibble.
#' @param exclude_terms Optional character vector of terms (e.g., Block 1 coupled terms)
#'   to assign directly to `"primary_topic"` before clustering remaining facets.
#' @return Named character vector mapping each term in `terms_df$term` to its
#'   automatically discovered cluster label (e.g., `"concept_fatal"`, `"concept_cases"`).
#' @export
discover_term_clusters <- function(terms_df, seeds_df, exclude_terms = character(0)) {
  if (nrow(terms_df) == 0L) {
    return(stats::setNames(character(0), character(0)))
  }

  all_terms <- terms_df$term
  out_clusters <- stats::setNames(rep("", length(all_terms)), all_terms)
  if (length(exclude_terms) > 0L) {
    is_b1 <- all_terms %in% exclude_terms
    out_clusters[is_b1] <- "primary_topic"
  }

  rem_idx <- which(!nzchar(out_clusters))
  m <- length(rem_idx)
  if (m == 0L) {
    return(out_clusters)
  }
  if (m == 1L) {
    slug <- gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", gsub("\\*$", "", all_terms[rem_idx])))
    out_clusters[rem_idx] <- paste0("concept_", if (nzchar(slug)) slug else "1")
    return(out_clusters)
  }

  rem_terms <- all_terms[rem_idx]
  rem_cov <- terms_df$cov_tiab[rem_idx]
  rem_prop <- terms_df$seed_prop[rem_idx]
  rem_scores <- terms_df$term_score[rem_idx]
  k <- nrow(seeds_df)

  ctx_sim <- .compute_term_context_similarity(rem_terms, rem_cov, seeds_df)

  edge_from <- integer(0)
  edge_to <- integer(0)
  edge_w <- numeric(0)

  for (i in seq_len(m - 1L)) {
    ci <- rem_cov[[i]]
    pi <- rem_prop[i]
    ti <- rem_terms[i]
    for (j in (i + 1L):m) {
      cj <- rem_cov[[j]]
      pj <- rem_prop[j]
      tj <- rem_terms[j]

      lex_stem <- shares_stem(ti, tj) || .is_code_family(ti, tj)
      cgram <- if (lex_stem) 1.0 else .char_ngram_sim(ti, tj, 4L)
      s_lex <- if (lex_stem) 1.0 else if (cgram >= 0.40) cgram else 0.0

      # Coverage complementarity: does the smaller term specifically fill gaps left by the larger term?
      if (sum(ci) >= sum(cj)) {
        c_big <- ci
        c_sml <- cj
      } else {
        c_big <- cj
        c_sml <- ci
      }
      n_miss <- sum(!c_big)
      n_sml <- sum(c_sml)
      s_comp <- if (n_miss > 0L && n_sml > 0L) {
        (sum((!c_big) & c_sml) / min(n_miss, n_sml)) * (min(n_sml, 3L) / 3.0)
      } else {
        0.0
      }

      s_ctx <- ctx_sim[i, j]

      # Prevent two distinct high-coverage orthogonal pillars (both >= 55% seed coverage, no lexical overlap)
      # from collapsing into one OR block unless their local window context is nearly identical
      both_pillars <- (pi >= 0.55 && pj >= 0.55 && !lex_stem)
      if (both_pillars && s_ctx < 0.72) {
        w <- 0.0
      } else {
        s_sem <- 0.50 * s_ctx + 0.40 * s_comp * (if (s_ctx >= 0.12) 1.0 else 0.40)
        w <- max(0.65 * s_lex, s_sem) + 0.20 * min(s_lex, s_ctx)
      }

      if (w >= 0.24) {
        edge_from <- c(edge_from, i)
        edge_to <- c(edge_to, j)
        edge_w <- c(edge_w, w)
      }
    }
  }

  g <- igraph::make_empty_graph(n = m, directed = FALSE)
  if (length(edge_w) > 0L) {
    g <- igraph::add_edges(g, as.vector(rbind(edge_from, edge_to)), attr = list(weight = edge_w))
    comm <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
    mem <- as.integer(igraph::membership(comm))
  } else {
    mem <- seq_len(m)
  }

  # Name each community after its highest-scoring representative term
  used_names <- "primary_topic"
  for (cid in unique(mem)) {
    grp <- which(mem == cid)
    best_local <- grp[which.max(rem_scores[grp])]
    rep_t <- gsub("\\*$", "", rem_terms[best_local])
    slug <- gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", rep_t))
    if (!nzchar(slug)) slug <- as.character(cid)
    base_nm <- paste0("concept_", substr(slug, 1L, 18L))
    nm <- base_nm
    suffix <- 2L
    while (nm %in% used_names) {
      nm <- paste0(base_nm, "_", suffix)
      suffix <- suffix + 1L
    }
    used_names <- c(used_names, nm)
    out_clusters[rem_idx[grp]] <- nm
  }

  out_clusters
}

#' Build orthogonal, unsupervised concept blocks from scored candidate terms
#'
#' Discovers Boolean `AND` concept blocks and their `OR` synonyms directly from
#' the seed papers and PubMed co-occurrence statistics without hardcoded word lists
#' or domain-specific ontologies:
#' 1. Identifies the primary high-lift anchor concept (`primary_topic`) and expands
#'    stem/code synonyms covering remaining seeds.
#' 2. Separates any universal domain-coupled context term (`domain_context`).
#' 3. Partitions remaining selective terms via [discover_term_clusters()] (SVD local
#'    context + MeSH/keyword co-assignment + coverage complementarity Louvain graph)
#'    and grows each concept block via greedy minimum-hit seed set-cover.
#'
#' @param terms_df Scored candidate tibble from [extract_candidate_terms()].
#' @param seeds_df Enriched seeds tibble.
#' @param config Optional config list (`custom_blocks`, `max_blocks`, `max_terms_per_block`).
#' @param max_blocks Maximum number of concept blocks.
#' @param max_terms_per_block Maximum terms per block.
#' @param probe_selectivity Logical or NULL (defaults to TRUE when `terms_df` was probed live).
#' @param eval_joint_fn Optional function `(t1, t2) -> integer` returning PubMed joint hit count.
#' @return List of block lists `list(name = ..., terms = ...)`.
#' @export
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
  if (nrow(terms_df) == 0L) {
    return(list())
  }

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
  if (length(b1_elig) == 0L) {
    return(list())
  }

  a1_idx <- b1_elig[which.max(terms_df$term_score[b1_elig])]
  a1 <- terms_df$term[a1_idx]
  a1_bg <- as.numeric(terms_df$bg_hits[a1_idx])
  a1_cov <- terms_df$cov_tiab[[a1_idx]]

  # Probe conditional PubMed co-occurrence P(a1 AND t) for top 25 stem-disjoint candidates
  j_hits_vec <- stats::setNames(numeric(0), character(0))
  probe_cands <- terms_df$term[
    terms_df$term != a1 &
      !vapply(terms_df$term, shares_stem, logical(1), t2 = a1)
  ]
  probe_cands <- utils::head(probe_cands, 25L)

  if (isTRUE(probe_selectivity) && is.function(eval_joint_fn) && length(probe_cands) > 0L) {
    jh_vals <- vapply(probe_cands, function(tc) {
      jh <- suppressWarnings(as.numeric(eval_joint_fn(a1, tc)))
      if (length(jh) == 1L && !is.na(jh)) jh else NA_real_
    }, numeric(1))
    j_hits_vec <- jh_vals[!is.na(jh_vals)]
  }

  get_jh <- function(tc) {
    if (tc %in% names(j_hits_vec)) {
      return(j_hits_vec[[tc]])
    }
    if (isTRUE(probe_selectivity) && is.function(eval_joint_fn)) {
      jh_val <- suppressWarnings(as.numeric(eval_joint_fn(a1, tc)))
      if (length(jh_val) == 1L && !is.na(jh_val)) {
        j_hits_vec[[tc]] <<- jh_val
        return(jh_val)
      }
    }
    NA_real_
  }

  # Identify terms taxonomically redundant with or coupled to Block 1 purely from lexical & PubMed conditional probability
  tax_parents <- names(j_hits_vec)[(j_hits_vec / max(a1_bg, 1)) > 0.72]
  coupled_mask <- vapply(seq_len(nrow(terms_df)), function(idx) {
    tc <- terms_df$term[idx]
    if (tc == a1 || shares_stem(tc, a1) || .is_code_family(a1, tc)) {
      return(TRUE)
    }
    if (length(tax_parents) > 0L && any(vapply(tax_parents, function(tp) .term_subsumes(tc, tp), logical(1)))) {
      return(TRUE)
    }
    c_bg <- as.numeric(terms_df$bg_hits[idx])
    if (tc %in% names(j_hits_vec)) {
      jh <- j_hits_vec[[tc]]
      # Taxonomic parent/context (e.g. >72% of all PubMed `a1` papers already contain `tc`)
      if ((jh / max(a1_bg, 1)) > 0.72) {
        return(TRUE)
      }
      # Narrow synonym/sub-phrase of `a1` (>25% of all PubMed `tc` papers also contain `a1`)
      if ((jh / max(c_bg, 1)) > 0.25) {
        return(TRUE)
      }
    }
    FALSE
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
        if (new_cov == 0L) {
          return(-1)
        }
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

  # Identify universal domain context qualifier coupled to Block 1 with >=85% seed coverage
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

  term_clusters <- discover_term_clusters(terms_df, seeds_df, exclude_terms = block1_coupled)

  if (max_blocks <= 1L) {
    attr(blocks, "block1_coupled") <- block1_coupled
    attr(blocks, "joint_hits") <- j_hits_vec
    attr(blocks, "term_clusters") <- term_clusters
    attr(blocks, "domain_qualifier") <- domain_qualifier
    return(blocks)
  }

  # --- Blocks 2..max_blocks: Unsupervised Graph Communities + Orthogonal Set-Cover ---
  is_selective_vs_a1 <- function(idx, max_bg = 4000000) {
    tc <- terms_df$term[idx]
    c_bg <- as.numeric(terms_df$bg_hits[idx])
    if (c_bg > max_bg) {
      return(FALSE)
    }
    if (tc %in% names(j_hits_vec)) {
      jh <- j_hits_vec[[tc]]
      return((jh / max(a1_bg, 1)) <= 0.72)
    }
    TRUE
  }

  cand_blocks <- list()
  used_terms <- b1_terms
  cluster_ids <- setdiff(unique(unname(term_clusters)), "primary_topic")

  for (cid in cluster_ids) {
    c_idx <- which(
      unname(term_clusters) == cid &
        !(terms_df$term %in% used_terms) &
        !(terms_df$term %in% block1_coupled) &
        terms_df$df_tiab >= min_anchor_df &
        terms_df$lift >= 5
    )
    c_idx <- Filter(function(idx) {
      !any(vapply(b1_terms, shares_stem, logical(1), t2 = terms_df$term[idx])) &&
        is_selective_vs_a1(idx)
    }, c_idx)
    if (length(c_idx) == 0L) next

    anchor_scores <- vapply(c_idx, function(idx) {
      tc <- terms_df$term[idx]
      jh <- if (tc %in% names(j_hits_vec)) j_hits_vec[[tc]] else NA_real_
      red <- if (!is.na(jh)) {
        log10(pmax(a1_bg / pmax(jh, 1), 1.1))
      } else {
        log10(pmax(38000000 / pmax(terms_df$bg_hits[idx], 100), 1.1))
      }
      (terms_df$df_tiab[idx] / max(k, 1L)) * red *
        (1 + 0.15 * (terms_df$n_words[idx] > 1L)) *
        (1 + 0.25 * endsWith(tc, "*"))
    }, numeric(1))

    a_idx <- c_idx[which.max(anchor_scores)]
    a_term <- terms_df$term[a_idx]
    a_bg <- terms_df$bg_hits[a_idx]
    a_df <- terms_df$df_tiab[a_idx]
    jh_a_known <- get_jh(a_term)
    max_syn_jh <- if (is.finite(jh_a_known)) max(0.08 * a1_bg, 1.15 * jh_a_known) else 0.55 * a1_bg
    blk_terms <- a_term
    blk_cov <- terms_df$cov_tiab[[a_idx]]

    # Grow block within the discovered community via greedy marginal uncovered seed coverage
    pool_idx <- which(
      unname(term_clusters) == cid &
        terms_df$term != a_term &
        !(terms_df$term %in% used_terms) &
        !(terms_df$term %in% block1_coupled) &
        terms_df$lift >= 5
    )
    pool_idx <- Filter(function(idx) {
      tc <- terms_df$term[idx]
      c_bg <- terms_df$bg_hits[idx]
      if (any(vapply(b1_terms, shares_stem, logical(1), t2 = tc))) {
        return(FALSE)
      }
      if (c_bg > 2000000 && a_bg < 1000000) {
        return(FALSE)
      }
      if (tc %in% names(j_hits_vec) && j_hits_vec[[tc]] > max_syn_jh) {
        return(FALSE)
      }
      is_selective_vs_a1(idx)
    }, pool_idx)

    a_cov_prop <- mean(blk_cov)
    dq_stems <- if (length(domain_qualifier) > 0L) term_stems(domain_qualifier) else character(0)
    b1_all_stems <- unique(c(unlist(lapply(b1_terms, term_stems)), dq_stems))
    is_pillar_anchor <- (a_df / max(k, 1L)) >= 0.65

    # Unified minimum-joint-hit set-cover growth (in-community synonyms + selective unigram complements for pillar anchors)
    cross_idx <- if (is_pillar_anchor) {
      which(
        unname(term_clusters) != cid &
          terms_df$term != a_term &
          !(terms_df$term %in% used_terms) &
          !(terms_df$term %in% block1_coupled) &
          terms_df$n_words == 1L &
          terms_df$bg_hits <= 250000 &
          terms_df$lift >= 8
      )
    } else {
      integer(0)
    }
    cross_idx <- Filter(function(idx) {
      tc <- terms_df$term[idx]
      !any(term_stems(tc) %in% b1_all_stems) && is_selective_vs_a1(idx, max_bg = 250000)
    }, cross_idx)

    unified_pool <- unique(c(pool_idx, cross_idx))
    while (length(blk_terms) < max_terms_per_block && length(unified_pool) > 0L && !all(blk_cov)) {
      rem_pool <- Filter(function(pi) {
        tc <- terms_df$term[pi]
        !(tc %in% blk_terms) &&
          !any(vapply(blk_terms, .term_subsumes, logical(1), t_narrow = tc)) &&
          any(terms_df$cov_tiab[[pi]] & !blk_cov)
      }, unified_pool)
      if (length(rem_pool) == 0L) break
      gains <- vapply(rem_pool, function(pi) {
        c_cov <- terms_df$cov_tiab[[pi]]
        new_cov <- sum(c_cov & !blk_cov)
        if (new_cov == 0L) {
          return(-1)
        }
        tc <- terms_df$term[pi]
        jh_tc <- get_jh(tc)
        if (is.finite(jh_tc) && jh_tc > max_syn_jh) {
          return(-1)
        }
        in_cid_bonus <- if (identical(unname(term_clusters[pi]), cid)) 1.25 else 1.0
        df_bonus <- 1 + 0.15 * log2(pmax(terms_df$df_tiab[pi], 1L))
        hit_cost <- if (is.finite(jh_tc)) {
          (pmax(jh_tc, 1) + 15)^0.55
        } else {
          (pmax(terms_df$bg_hits[pi] / 250, 5) + 15)^0.55
        }
        ((new_cov^1.35) * df_bonus * in_cid_bonus) / hit_cost
      }, numeric(1))
      if (max(gains) <= 0) break
      best_pi <- rem_pool[which.max(gains)]
      blk_terms <- c(blk_terms, terms_df$term[best_pi])
      blk_cov <- blk_cov | terms_df$cov_tiab[[best_pi]]
    }

    # If anchor covered all seeds in community alone, add 1 narrower high-precision alternative for specificity tiers
    if (length(blk_terms) == 1L && length(pool_idx) > 0L) {
      narrow_cands <- pool_idx[terms_df$bg_hits[pool_idx] < a_bg & terms_df$df_tiab[pool_idx] >= ceiling(0.4 * a_df)]
      if (length(narrow_cands) > 0L) {
        best_n <- narrow_cands[which.max(terms_df$term_score[narrow_cands])]
        blk_terms <- c(blk_terms, terms_df$term[best_n])
      }
    }

    narrow_pool <- c_idx[terms_df$bg_hits[c_idx] < a_bg & terms_df$df_tiab[c_idx] >= ceiling(0.45 * a_df)]
    narrow_term <- if (length(narrow_pool) > 0L) {
      terms_df$term[narrow_pool[which.min(terms_df$bg_hits[narrow_pool])]]
    } else {
      character(0)
    }

    # Score block by anchor/community coverage and the anchor's true selectivity against Block 1
    jh_anchor <- if (a_term %in% names(j_hits_vec)) j_hits_vec[[a_term]] else NA_real_
    if (!is.finite(jh_anchor) || jh_anchor >= a1_bg) {
      est_jh <- a1_bg * min((a_bg / 38000000) * 8, 0.85)
      jh_anchor <- max(est_jh, 25)
    }
    eff_cov <- 0.60 * a_cov_prop + 0.40 * mean(blk_cov)
    blk_score <- (eff_cov^1.35) * log10(pmax(a1_bg / pmax(jh_anchor, 10), 1.5))

    # Name block cleanly after its chosen anchor term
    slug_a <- gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", gsub("\\*$", "", a_term)))
    blk_nm <- paste0("concept_", substr(if (nzchar(slug_a)) slug_a else "facet", 1L, 18L))

    cand_blocks[[length(cand_blocks) + 1L]] <- list(
      name = blk_nm,
      cluster_id = cid,
      terms = blk_terms,
      cov = blk_cov,
      narrow_term = narrow_term,
      score = blk_score
    )
  }

  # Greedily select top (max_blocks - 1) orthogonal blocks (avoiding blocks with overlapping stems)
  if (length(cand_blocks) > 0L) {
    ord_cb <- order(vapply(cand_blocks, `[[`, numeric(1), "score"), decreasing = TRUE)
    cand_blocks <- cand_blocks[ord_cb]
    picked_cb <- list()
    picked_stems <- unlist(lapply(b1_terms, term_stems))
    used_blk_names <- "primary_topic"

    for (cb in cand_blocks) {
      if (length(picked_cb) >= (max_blocks - 1L)) break
      cb_stems <- unlist(lapply(cb$terms, term_stems))
      if (length(intersect(cb_stems, picked_stems)) > 0L) next
      nm <- cb$name
      sfx <- 2L
      while (nm %in% used_blk_names) {
        nm <- paste0(cb$name, "_", sfx)
        sfx <- sfx + 1L
      }
      used_blk_names <- c(used_blk_names, nm)
      picked_stems <- unique(c(picked_stems, cb_stems))
      blk_item <- list(name = nm, terms = cb$terms)
      if (!is.null(cb$narrow_term) && length(cb$narrow_term) > 0L) {
        attr(blk_item, "narrow_term") <- cb$narrow_term
      }
      attr(blk_item, "cluster_id") <- cb$cluster_id
      picked_cb[[length(picked_cb) + 1L]] <- blk_item
    }
    blocks <- c(blocks, picked_cb)
  }

  attr(blocks, "block1_coupled") <- block1_coupled
  attr(blocks, "joint_hits") <- j_hits_vec
  attr(blocks, "term_clusters") <- term_clusters
  attr(blocks, "domain_qualifier") <- domain_qualifier
  blocks
}

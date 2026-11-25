# R/03_enrich_seeds.R
# Resolve BibTeX seeds to PubMed PMIDs, fetch missing abstracts,
# MeSH descriptors, and Author Keywords, with token-Jaccard title verification.

#' Compute alphanumeric lowercase token Jaccard similarity between two titles
#'
#' @param t1 First title string.
#' @param t2 Second title string.
#' @return Numeric Jaccard similarity in `[0, 1]`.
title_jaccard <- function(t1, t2) {
  tok <- function(s) {
    s <- tolower(as.character(s %||% ""))
    s <- stringr::str_replace_all(s, "[^a-z0-9]+", " ")
    w <- unlist(stringr::str_split(stringr::str_squish(s), " "))
    w[nzchar(w) & nchar(w) >= 2L]
  }
  a <- unique(tok(t1))
  b <- unique(tok(t2))
  if (length(a) == 0L || length(b) == 0L) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

#' Clean title string for search queries (strip quotes/brackets/colons that break E-utilities)
#'
#' @param title Input article title string.
#' @return Cleaned title string.
clean_title_for_query <- function(title) {
  s <- as.character(title %||% "")
  s <- stringr::str_replace_all(s, '["\\[\\]()]', " ")
  stringr::str_squish(s)
}

#' First N significant words from a title for fallback matching
#'
#' @param title Input article title string.
#' @param max_words Maximum number of significant words to retain.
#' @return Space-separated significant words string.
title_sig_words <- function(title, max_words = 12L) {
  s <- tolower(clean_title_for_query(title))
  s <- stringr::str_replace_all(s, "[^a-z0-9 ]+", " ")
  w <- unlist(stringr::str_split(stringr::str_squish(s), " "))
  stops <- c("the", "and", "for", "with", "from", "that", "this", "into", "over", "upon", "a", "an", "of", "in", "on", "to", "by", "at", "as", "is", "are", "was", "were")
  w <- w[nzchar(w) & !(w %in% stops)]
  paste(utils::head(w, max_words), collapse = " ")
}

#' Enrich BibTeX seed papers via PubMed
#'
#' @param seeds_df Tibble returned by `parse_bibtex()`.
#' @param email Optional email for NCBI polite pool.
#' @param ncbi_key Optional NCBI API key.
#' @param cache_dir Cache directory path.
#' @param use_cache Logical.
#' @param pm_esearch_fn Optional dependency injection for `pubmed_esearch`.
#' @param pm_efetch_fn Optional dependency injection for `pubmed_efetch`.
#' @param progress_cb Optional progress callback function.
#' @return Enriched tibble with `full_title`, `full_abstract`, `full_text`,
#'   `pmid`, `in_pubmed`, `mesh_terms`, `all_keywords`.
enrich_seeds <- function(seeds_df,
                         email = NULL,
                         ncbi_key = NULL,
                         cache_dir = ".cache",
                         use_cache = TRUE,
                         pm_esearch_fn = pubmed_esearch,
                         pm_efetch_fn = pubmed_efetch,
                         progress_cb = NULL) {
  n <- nrow(seeds_df)
  pmids <- as.character(seeds_df$pmid %||% rep("", n))
  pmids[is.na(pmids)] <- ""

  # 1. Resolve missing PMIDs via DOI or Title (+ optional publication year)
  for (i in seq_len(n)) {
    if (is.function(progress_cb)) progress_cb(sprintf("Resolving seed %d/%d", i, n))
    if (nzchar(pmids[i])) next
    doi_i <- seeds_df$doi[i]
    title_i <- seeds_df$title[i]
    year_i <- as.character(seeds_df$year[i] %||% "")

    resolved <- ""
    if (nzchar(doi_i)) {
      res <- pm_esearch_fn(
        term = paste0('"', doi_i, '"[doi]'), retmax = 1L,
        email = email, api_key = ncbi_key, cache_dir = cache_dir, use_cache = use_cache
      )
      if (isTRUE(res$ok) && length(res$idlist) >= 1L) {
        resolved <- res$idlist[1L]
      }
    }

    if (!nzchar(resolved) && nzchar(title_i)) {
      q_exact <- paste0('"', clean_title_for_query(title_i), '"[title]')
      res_t <- pm_esearch_fn(
        term = q_exact, retmax = 3L,
        email = email, api_key = ncbi_key, cache_dir = cache_dir, use_cache = use_cache
      )
      cands <- if (isTRUE(res_t$ok)) res_t$idlist else character(0)

      sig <- title_sig_words(title_i, 12L)
      if (length(cands) == 0L && nzchar(sig)) {
        sig_q <- paste(paste0(strsplit(sig, " ")[[1]], "[title]"), collapse = " AND ")
        # Try with 4-digit year constraint first to avoid chronological truncation on common title words
        if (grepl("^[12][0-9]{3}$", year_i)) {
          res_ty <- pm_esearch_fn(
            term = paste0(sig_q, " AND ", year_i, "[dp]"), retmax = 5L,
            email = email, api_key = ncbi_key, cache_dir = cache_dir, use_cache = use_cache
          )
          if (isTRUE(res_ty$ok) && length(res_ty$idlist) > 0L) {
            cands <- res_ty$idlist
          }
        }
        if (length(cands) == 0L) {
          res_t2 <- pm_esearch_fn(
            term = sig_q, retmax = 5L,
            email = email, api_key = ncbi_key, cache_dir = cache_dir, use_cache = use_cache
          )
          cands <- if (isTRUE(res_t2$ok)) res_t2$idlist else character(0)
        }
      }

      if (length(cands) > 0L) {
        cand_meta <- pm_efetch_fn(
          pmids = cands, email = email, api_key = ncbi_key,
          cache_dir = cache_dir, use_cache = use_cache
        )
        if (nrow(cand_meta) > 0L) {
          sims <- vapply(cand_meta$pm_title, function(pt) title_jaccard(title_i, pt), numeric(1))
          best_idx <- which.max(sims)
          if (length(best_idx) == 1L && sims[best_idx] >= 0.75) {
            resolved <- cand_meta$pmid[best_idx]
          }
        }
      }
    }
    pmids[i] <- resolved
  }

  # 2. Batch efetch for all resolved PMIDs
  pm_meta <- pm_efetch_fn(
    pmids = pmids[nzchar(pmids)], email = email, api_key = ncbi_key,
    cache_dir = cache_dir, use_cache = use_cache
  )

  # 3. Merge metadata cleanly
  pm_map <- if (nrow(pm_meta) > 0L) {
    split(pm_meta, pm_meta$pmid)
  } else {
    list()
  }

  full_title <- character(n)
  full_abstract <- character(n)
  full_text <- character(n)
  mesh_list <- vector("list", n)
  kw_list <- vector("list", n)

  for (i in seq_len(n)) {
    p_row <- if (nzchar(pmids[i]) && !is.null(pm_map[[pmids[i]]])) pm_map[[pmids[i]]][1L, ] else NULL
    t_bib <- seeds_df$title[i]
    t_pm <- if (!is.null(p_row)) p_row$pm_title[1L] else ""
    full_title[i] <- dplyr::case_when(
      nzchar(t_bib) ~ t_bib,
      nzchar(t_pm) ~ t_pm,
      TRUE ~ ""
    )

    a_bib <- seeds_df$abstract[i]
    a_pm <- if (!is.null(p_row)) p_row$pm_abstract[1L] else ""
    # Prefer PubMed abstract when longer/complete or when BibTeX abstract is empty
    full_abstract[i] <- dplyr::case_when(
      nzchar(a_pm) && nchar(a_pm) >= nchar(a_bib) ~ a_pm,
      nzchar(a_bib) ~ a_bib,
      nzchar(a_pm) ~ a_pm,
      TRUE ~ ""
    )

    m_terms <- if (!is.null(p_row) && length(p_row$pm_mesh[[1L]]) > 0L) p_row$pm_mesh[[1L]] else character(0)
    k_bib <- split_keywords(seeds_df$keywords[i])
    k_pm <- if (!is.null(p_row) && length(p_row$pm_keywords[[1L]]) > 0L) p_row$pm_keywords[[1L]] else character(0)

    mesh_list[[i]] <- unique(m_terms[nzchar(m_terms)])
    kw_list[[i]] <- unique(c(k_bib, k_pm)[nzchar(c(k_bib, k_pm))])
    # Note: PubMed [tiab] indexes Title + Abstract + Author Keywords (k_pm)
    full_text[i] <- stringr::str_squish(paste(full_title[i], full_abstract[i], paste(k_pm, collapse = " ")))
  }

  out <- seeds_df |>
    dplyr::mutate(
      pmid = pmids,
      in_pubmed = nzchar(pmids),
      full_title = full_title,
      full_abstract = full_abstract,
      full_text = full_text,
      mesh_terms = mesh_list,
      all_keywords = kw_list
    )

  n_pm <- sum(out$in_pubmed)
  n_abs <- sum(nzchar(out$full_abstract))
  cli::cli_inform(c(
    "i" = "Seed resolution: {n} total | {n_pm} in PubMed | {n_abs} with abstracts."
  ))

  out
}

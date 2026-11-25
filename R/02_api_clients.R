# R/02_api_clients.R
# Disk-cached, rate-limited clients for PubMed (NCBI E-utilities).

.rate_state <- new.env(parent = emptyenv())
.rate_state$pubmed_last <- 0
.rate_state$pubmed_interval <- 0.45
.rate_state$pubmed_429_hits <- 0L

.respect_rate_limit <- function(min_interval = 0.45) {
  now <- proc.time()[["elapsed"]]
  last <- .rate_state$pubmed_last %||% 0
  wait <- min_interval - (now - last)
  if (wait > 0) Sys.sleep(wait)
  .rate_state$pubmed_last <- proc.time()[["elapsed"]]
}

.mark_request_done <- function() {
  .rate_state$pubmed_last <- proc.time()[["elapsed"]]
}

#' Build deterministic cache key
#'
#' @param prefix String prefix identifying the cache namespace.
#' @param payload Arbitrary R object hashed with MD5.
#' @return Character cache key.
cache_key <- function(prefix, payload) {
  paste0(prefix, "_", digest::digest(payload, algo = "md5"))
}

#' Read cached RDS value if available
#'
#' @param cache_dir Directory containing `.rds` cache files.
#' @param key Cache key string returned by [cache_key()].
#' @param use_cache Logical; whether to read from disk cache.
#' @return Cached R object or `NULL` on cache miss.
cache_get <- function(cache_dir, key, use_cache = TRUE) {
  if (!isTRUE(use_cache) || is.null(cache_dir) || !nzchar(cache_dir)) return(NULL)
  path <- file.path(cache_dir, paste0(key, ".rds"))
  if (file.exists(path)) {
    return(tryCatch(readRDS(path), error = function(e) NULL))
  }
  NULL
}

#' Write cached RDS value
#'
#' @param cache_dir Directory containing `.rds` cache files.
#' @param key Cache key string returned by [cache_key()].
#' @param value R object to serialize to disk.
#' @param use_cache Logical; whether to write to disk cache.
#' @return Invisibly returns `value`.
cache_set <- function(cache_dir, key, value, use_cache = TRUE) {
  if (!isTRUE(use_cache) || is.null(cache_dir) || !nzchar(cache_dir)) return(invisible(NULL))
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(cache_dir, paste0(key, ".rds"))
  saveRDS(value, path)
  invisible(value)
}

.is_short_transient <- function(resp) {
  st <- httr2::resp_status(resp)
  if (st == 429L) {
    ra <- suppressWarnings(httr2::resp_retry_after(resp))
    return(is.na(ra) || ra <= 5)
  }
  st %in% c(500L, 502L, 503L, 504L)
}

#' Query PubMed E-utilities esearch via HTTP POST
#'
#' @param term PubMed query string.
#' @param retmax Maximum number of PMIDs to return (0 for count only).
#' @param email Optional email for NCBI tool identification.
#' @param api_key Optional NCBI API key.
#' @param cache_dir Directory for RDS cache.
#' @param use_cache Logical; whether to read/write cache.
#' @return List with `count` (integer), `idlist` (character), `warnings` (character), `ok` (logical).
pubmed_esearch <- function(term, retmax = 0L, email = NULL, api_key = NULL,
                           cache_dir = ".cache", use_cache = TRUE) {
  retmax <- as.integer(retmax)
  ck <- cache_key("pm_esearch_v1", list(term = term, retmax = retmax))
  cached <- cache_get(cache_dir, ck, use_cache = use_cache)
  if (!is.null(cached)) return(cached)

  # If count-only (retmax == 0) is requested, reuse any cached full-idlist esearch for the same term
  if (retmax == 0L && isTRUE(use_cache)) {
    for (alt_max in c(2500L, 1200L, 600L)) {
      ck_alt <- cache_key("pm_esearch_v1", list(term = term, retmax = alt_max))
      alt_val <- cache_get(cache_dir, ck_alt, use_cache = TRUE)
      if (!is.null(alt_val) && isTRUE(alt_val$ok)) {
        return(list(count = alt_val$count, idlist = character(0), warnings = alt_val$warnings, ok = TRUE))
      }
    }
  }

  has_key <- !is.null(api_key) && nzchar(api_key)
  interval <- if (has_key) 0.14 else (.rate_state$pubmed_interval %||% 0.45)
  .respect_rate_limit(interval)

  body <- list(
    db = "pubmed",
    term = term,
    retmax = as.character(retmax),
    retmode = "json",
    tool = "r_query_optimizer"
  )
  if (!is.null(email) && nzchar(email)) body$email <- email
  if (has_key) body$api_key <- api_key

  out <- tryCatch({
    on.exit(.mark_request_done(), add = TRUE)
    t_start <- proc.time()[["elapsed"]]
    req <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi") |>
      httr2::req_body_form(!!!body) |>
      httr2::req_retry(max_tries = 4, is_transient = .is_short_transient, max_seconds = 20, backoff = ~ 1.2 * (2 ^ .x)) |>
      httr2::req_timeout(30)

    resp <- httr2::req_perform(req)
    elapsed <- proc.time()[["elapsed"]] - t_start
    if (!has_key && (httr2::resp_status(resp) == 429L || elapsed > 2.5)) {
      .rate_state$pubmed_interval <- min((.rate_state$pubmed_interval %||% 0.45) + 0.15, 1.0)
      .rate_state$pubmed_429_hits <- (.rate_state$pubmed_429_hits %||% 0L) + 1L
      Sys.sleep(1.0)
    }
    parsed <- httr2::resp_body_json(resp, simplifyVector = TRUE)
    esr <- parsed$esearchresult
    if (is.null(esr) || !is.null(esr$ERROR) || is.null(esr$count)) {
      err_msg <- as.character(esr$ERROR %||% parsed$error %||% "Missing esearchresult count")
      return(list(count = NA_integer_, idlist = character(0), warnings = err_msg, ok = FALSE))
    }

    cnt <- suppressWarnings(as.integer(esr$count))
    if (is.na(cnt)) {
      return(list(count = NA_integer_, idlist = character(0), warnings = "Non-numeric esearch count", ok = FALSE))
    }
    ids <- as.character(esr$idlist %||% character(0))
    warns <- c(
      unlist(esr$warninglist$phrasesignored %||% character(0)),
      unlist(esr$warninglist$quotedphrasesnotfound %||% character(0)),
      unlist(esr$errorlist$phrasesnotfound %||% character(0))
    )
    list(count = cnt, idlist = ids, warnings = as.character(warns), ok = TRUE)
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (!has_key && grepl("429|rate limit", msg, ignore.case = TRUE)) {
      .rate_state$pubmed_interval <- min((.rate_state$pubmed_interval %||% 0.45) + 0.15, 1.0)
      .rate_state$pubmed_429_hits <- (.rate_state$pubmed_429_hits %||% 0L) + 1L
      Sys.sleep(1.0)
    }
    list(count = NA_integer_, idlist = character(0), warnings = msg, ok = FALSE)
  })

  if (isTRUE(out$ok)) cache_set(cache_dir, ck, out, use_cache = use_cache)
  out
}

#' Fetch article metadata (title, abstract, MeSH, keywords) from PubMed efetch XML
#'
#' Uses both batch-level and per-PMID disk caching so previously fetched articles
#' never need to be re-downloaded when seed lists or noise samples change.
#'
#' @param pmids Character vector of PMIDs.
#' @param email Optional email for NCBI tool identification.
#' @param api_key Optional NCBI API key.
#' @param cache_dir Directory for RDS cache.
#' @param use_cache Logical; whether to read/write cache.
#' @return Tibble with `pmid`, `pm_title`, `pm_abstract`, `pm_mesh`, `pm_keywords`.
pubmed_efetch <- function(pmids, email = NULL, api_key = NULL,
                          cache_dir = ".cache", use_cache = TRUE) {
  pmids <- unique(as.character(pmids[!is.na(pmids) & nzchar(pmids)]))
  empty_tbl <- tibble::tibble(
    pmid = character(0),
    pm_title = character(0),
    pm_abstract = character(0),
    pm_mesh = list(),
    pm_keywords = list()
  )
  if (length(pmids) == 0L) return(empty_tbl)

  # 1. Check exact batch cache first (fast path)
  has_key <- !is.null(api_key) && nzchar(api_key)
  batches <- split(pmids, ceiling(seq_along(pmids) / 100L))
  rows <- list()
  uncached_pmids <- character(0)

  for (batch in batches) {
    ck <- cache_key("pm_efetch_v1", sort(batch))
    cached <- cache_get(cache_dir, ck, use_cache = use_cache)
    if (!is.null(cached) && nrow(cached) > 0L) {
      rows[[length(rows) + 1L]] <- cached
    } else {
      uncached_pmids <- c(uncached_pmids, batch)
    }
  }

  # 2. Check per-PMID cache (and harvest any existing pm_efetch_v1_* files once if needed)
  if (length(uncached_pmids) > 0L && isTRUE(use_cache) && !is.null(cache_dir) && nzchar(cache_dir)) {
    still_missing <- character(0)
    for (id in uncached_pmids) {
      art_row <- cache_get(cache_dir, paste0("pm_art_v1_", id), use_cache = TRUE)
      if (!is.null(art_row) && nrow(art_row) > 0L) {
        rows[[length(rows) + 1L]] <- art_row
      } else {
        still_missing <- c(still_missing, id)
      }
    }
    # If some are still missing, scan existing batch cache files in cache_dir to see if they were already fetched in another batch
    if (length(still_missing) > 0L && dir.exists(cache_dir)) {
      batch_files <- list.files(cache_dir, pattern = "^pm_efetch_v1_.*\\.rds$", full.names = TRUE)
      for (bf in batch_files) {
        if (length(still_missing) == 0L) break
        b_df <- tryCatch(readRDS(bf), error = function(e) NULL)
        if (!is.null(b_df) && nrow(b_df) > 0L) {
          hit_rows <- b_df[b_df$pmid %in% still_missing, , drop = FALSE]
          if (nrow(hit_rows) > 0L) {
            rows[[length(rows) + 1L]] <- hit_rows
            for (ri in seq_len(nrow(hit_rows))) {
              one_r <- hit_rows[ri, , drop = FALSE]
              cache_set(cache_dir, paste0("pm_art_v1_", one_r$pmid[1L]), one_r, use_cache = TRUE)
            }
            still_missing <- setdiff(still_missing, hit_rows$pmid)
          }
        }
      }
    }
    uncached_pmids <- still_missing
  }

  # 3. Fetch any remaining uncached PMIDs from NCBI E-utilities in batches of 150
  if (length(uncached_pmids) > 0L) {
    net_batches <- split(uncached_pmids, ceiling(seq_along(uncached_pmids) / 150L))
    for (batch in net_batches) {
      ck <- cache_key("pm_efetch_v1", sort(batch))
      interval <- if (has_key) 0.14 else (.rate_state$pubmed_interval %||% 0.45)
      .respect_rate_limit(interval)
      body <- list(
        db = "pubmed",
        id = paste(batch, collapse = ","),
        retmode = "xml",
        tool = "r_query_optimizer"
      )
      if (!is.null(email) && nzchar(email)) body$email <- email
      if (has_key) body$api_key <- api_key

      batch_tbl <- tryCatch({
        on.exit(.mark_request_done(), add = TRUE)
        t_start <- proc.time()[["elapsed"]]
        req <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi") |>
          httr2::req_body_form(!!!body) |>
          httr2::req_retry(max_tries = 4, is_transient = .is_short_transient, max_seconds = 20, backoff = ~ 1.2 * (2 ^ .x)) |>
          httr2::req_timeout(45)
        resp <- httr2::req_perform(req)
        elapsed <- proc.time()[["elapsed"]] - t_start
        if (!has_key && (httr2::resp_status(resp) == 429L || elapsed > 2.5)) {
          .rate_state$pubmed_interval <- min((.rate_state$pubmed_interval %||% 0.45) + 0.15, 1.0)
          .rate_state$pubmed_429_hits <- (.rate_state$pubmed_429_hits %||% 0L) + 1L
          Sys.sleep(1.0)
        }
        xml_txt <- httr2::resp_body_string(resp)
        doc <- xml2::read_xml(xml_txt)
        articles <- xml2::xml_find_all(doc, ".//PubmedArticle")

        if (length(articles) == 0L) {
          empty_tbl
        } else {
          tibble::tibble(
            pmid = vapply(articles, function(art) {
              xml2::xml_text(xml2::xml_find_first(art, ".//MedlineCitation/PMID"))
            }, character(1)),
            pm_title = vapply(articles, function(art) {
              t_node <- xml2::xml_find_first(art, ".//ArticleTitle")
              if (inherits(t_node, "xml_missing")) "" else stringr::str_squish(xml2::xml_text(t_node))
            }, character(1)),
            pm_abstract = vapply(articles, function(art) {
              abs_nodes <- xml2::xml_find_all(art, ".//Abstract/AbstractText")
              if (length(abs_nodes) == 0L) "" else stringr::str_squish(paste(xml2::xml_text(abs_nodes), collapse = " "))
            }, character(1)),
            pm_mesh = lapply(articles, function(art) {
              nodes <- xml2::xml_find_all(art, ".//MeshHeadingList/MeshHeading/DescriptorName")
              unique(tolower(stringr::str_squish(xml2::xml_text(nodes))))
            }),
            pm_keywords = lapply(articles, function(art) {
              nodes <- xml2::xml_find_all(art, ".//KeywordList/Keyword")
              unique(tolower(stringr::str_squish(xml2::xml_text(nodes))))
            })
          )
        }
      }, error = function(e) {
        warning("PubMed efetch failed for batch: ", conditionMessage(e), call. = FALSE)
        empty_tbl
      })

      if (nrow(batch_tbl) > 0L) {
        cache_set(cache_dir, ck, batch_tbl, use_cache = use_cache)
        for (ri in seq_len(nrow(batch_tbl))) {
          one_r <- batch_tbl[ri, , drop = FALSE]
          cache_set(cache_dir, paste0("pm_art_v1_", one_r$pmid[1L]), one_r, use_cache = use_cache)
        }
      }
      rows[[length(rows) + 1L]] <- batch_tbl
    }
  }

  out <- dplyr::bind_rows(rows)
  if (nrow(out) > 0L) {
    out <- out[!duplicated(out$pmid), , drop = FALSE]
  }
  out
}

#' Evaluate a PubMed query: total hit count + which seed PMIDs are retrieved
#'
#' Requests up to 2500 PMIDs in the initial `esearch` call so that queries with
#' <= 2500 hits determine seed recall in a single HTTP request instead of two.
#'
#' @param pm_query PubMed query string.
#' @param seed_pmids Character vector of seed PMIDs.
#' @param email Optional email for NCBI tool identification.
#' @param api_key Optional NCBI API key.
#' @param cache_dir Directory for RDS cache.
#' @param use_cache Logical; whether to read/write cache.
#' @return List with `pm_hits`, `pm_found`, `pm_recall`, `ok`.
pubmed_eval_query <- function(pm_query, seed_pmids, email = NULL, api_key = NULL,
                              cache_dir = ".cache", use_cache = TRUE) {
  seed_pmids <- unique(as.character(seed_pmids[!is.na(seed_pmids) & nzchar(seed_pmids)]))

  # Check if legacy retmax=0L cache exists first to avoid unnecessary network calls on cached runs
  ck0 <- cache_key("pm_esearch_v1", list(term = pm_query, retmax = 0L))
  cached0 <- cache_get(cache_dir, ck0, use_cache = use_cache)
  res1 <- if (!is.null(cached0) && isTRUE(cached0$ok)) {
    cached0
  } else {
    pubmed_esearch(
      term = pm_query, retmax = 2500L, email = email, api_key = api_key,
      cache_dir = cache_dir, use_cache = use_cache
    )
  }
  if (!isTRUE(res1$ok)) {
    return(list(pm_hits = NA_integer_, pm_found = character(0), pm_recall = NA_real_, ok = FALSE))
  }
  total_hits <- res1$count
  found_pmids <- character(0)

  if (length(seed_pmids) > 0L && total_hits > 0L) {
    if (length(res1$idlist) >= total_hits) {
      # Fast single-request path: all matching PMIDs were returned in res1$idlist
      found_pmids <- intersect(seed_pmids, res1$idlist)
    } else {
      batches <- split(seed_pmids, ceiling(seq_along(seed_pmids) / 200L))
      for (batch in batches) {
        uid_clause <- paste0(batch, "[uid]", collapse = " OR ")
        check_q <- paste0("(", pm_query, ") AND (", uid_clause, ")")
        res2 <- pubmed_esearch(
          term = check_q, retmax = length(batch), email = email, api_key = api_key,
          cache_dir = cache_dir, use_cache = use_cache
        )
        if (!isTRUE(res2$ok)) {
          return(list(pm_hits = total_hits, pm_found = character(0), pm_recall = NA_real_, ok = FALSE))
        }
        found_pmids <- union(found_pmids, res2$idlist)
      }
    }
  }

  recall <- if (length(seed_pmids) > 0L) length(found_pmids) / length(seed_pmids) else NA_real_
  list(pm_hits = total_hits, pm_found = found_pmids, pm_recall = recall, ok = TRUE)
}

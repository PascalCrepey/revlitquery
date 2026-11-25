# R/06_optimizer.R
# Local candidate generation, live PubMed evaluation, concept-grouped exclusion
# refinement, greedy Pareto frontier search, and trade-off selection.

#' Estimate local seed coverage and rough PubMed hit count for a candidate query
#'
#' @param q Canonical `search_query`.
#' @param terms_df Scored candidate terms tibble.
#' @param seeds_df Enriched seeds tibble.
#' @return List with `local_cov`, `local_found`, `local_recall`, `est_hits`.
estimate_local_query <- function(q, terms_df, seeds_df) {
  n_seeds <- nrow(seeds_df)
  raw_texts <- as.character(seeds_df$full_text %||% paste(seeds_df$title, seeds_df$abstract))
  raw_titles <- as.character(seeds_df$full_title %||% seeds_df$title)
  resolvable_mask <- if (any(seeds_df$in_pubmed %||% FALSE)) as.logical(seeds_df$in_pubmed) else rep(TRUE, n_seeds)
  n_res <- max(sum(resolvable_mask), 1L)

  if (length(q$blocks) == 0L) {
    return(list(local_cov = rep(FALSE, n_seeds), local_found = 0L, local_recall = 0, est_hits = Inf))
  }

  n_pm_corpus <- 38000000
  overall_cov <- rep(TRUE, n_seeds)
  block_probs <- numeric(length(q$blocks))
  has_domain_ctx <- any(vapply(q$blocks, function(b) identical(b$name, "domain_context"), logical(1)))

  for (bi in seq_along(q$blocks)) {
    blk <- q$blocks[[bi]]
    blk_cov <- rep(FALSE, n_seeds)
    blk_bg <- 0

    for (t_str in blk$terms) {
      idx <- match(t_str, terms_df$term)
      if (!is.na(idx)) {
        t_cov <- if (blk$field == "title") terms_df$cov_title[[idx]] else terms_df$cov_tiab[[idx]]
        t_bg <- as.numeric(terms_df$bg_hits[idx])
      } else {
        t_cov <- term_matches_text(t_str, if (blk$field == "title") raw_titles else raw_texts)
        t_bg <- if (grepl(" ", t_str)) 20000 else 200000
      }
      # Short 2-char alphanumeric codes (e.g. "h5") have non-domain background hits unless paired with domain_context
      if (bi == 1L && nchar(t_str) <= 3L && grepl("^[a-z]+[0-9]+$", t_str)) {
        t_bg <- if (has_domain_ctx) t_bg * 0.12 else t_bg * 0.85
      }
      if (blk$field == "title") t_bg <- t_bg * 0.18
      blk_cov <- blk_cov | t_cov
      blk_bg <- blk_bg + t_bg
    }

    overall_cov <- overall_cov & blk_cov
    raw_prob <- min(max(blk_bg, 1) / n_pm_corpus, 1)
    if (bi == 1L) {
      block_probs[bi] <- raw_prob
    } else if (identical(blk$name, "domain_context")) {
      block_probs[bi] <- 0.90
    } else if (identical(blk$name, "population_host")) {
      block_probs[bi] <- max(raw_prob, if (blk$field == "title") 0.18 else 0.45)
    } else if (blk$name %in% c("outcome_severity", "outcome_or_epi")) {
      block_probs[bi] <- min(max(raw_prob * 8.5, if (blk$field == "title") 0.012 else 0.035), 0.35)
    } else {
      block_probs[bi] <- min(max(raw_prob * 6.0, 0.03), 0.40)
    }
  }

  excl_mult <- 1.0
  if (length(q$exclude_terms) > 0L) {
    for (ex in q$exclude_terms) {
      p_ex <- .parse_exclude_term(ex)
      tgt <- if (identical(p_ex$field, "title")) raw_titles else raw_texts
      overall_cov <- overall_cov & !term_matches_text(p_ex$term, tgt)
      excl_mult <- excl_mult * (if (identical(p_ex$field, "title")) 0.95 else 0.89)
    }
  }

  cov_res <- overall_cov & resolvable_mask
  est_hits <- n_pm_corpus * prod(block_probs) * max(excl_mult, 0.25)
  list(
    local_cov = overall_cov,
    local_found = sum(cov_res),
    local_recall = sum(cov_res) / n_res,
    est_hits = est_hits
  )
}

# Build specificity tiers for a single concept block
.block_tiers <- function(blk, terms_df, max_terms = 5L) {
  terms <- blk$terms
  if (length(terms) == 0L) {
    return(list())
  }
  idx <- match(terms, terms_df$term)
  bgs <- ifelse(is.na(idx), 50000, terms_df$bg_hits[idx])
  dfs <- ifelse(is.na(idx), 1, terms_df$df_tiab[idx])

  # `terms` is already in greedy cumulative coverage order (anchor first, then best complements)
  ord_narrow <- terms[order(bgs, -dfs)]
  narrow_attr <- attr(blk, "narrow_term") %||% character(0)

  tiers <- list(
    terms[1L],
    if (length(narrow_attr) > 0L) narrow_attr[1L] else ord_narrow[1L],
    ord_narrow[1L],
    utils::head(terms, min(2L, length(terms))),
    utils::head(ord_narrow, min(2L, length(ord_narrow))),
    utils::head(terms, min(3L, length(terms))),
    utils::head(terms, min(4L, length(terms))),
    utils::head(terms, min(max_terms, length(terms)))
  )
  sigs <- vapply(tiers, function(v) paste(sort(v), collapse = "|"), character(1))
  tiers[!duplicated(sigs)]
}

#' Generate diverse initial Stage 1 query candidates
#'
#' @param blocks List of discovered concept blocks.
#' @param terms_df Scored candidate terms tibble.
#' @param seeds_df Enriched seeds tibble.
#' @param config Optional config list.
#' @param min_recall Target minimum recall used to guide diversity selection.
#' @param max_hits Optional maximum PubMed hit count cap.
#' @param max_gs_chars Maximum Google Scholar query characters (default 256).
#' @param max_candidates Maximum initial candidates to return.
#' @return List of `search_query` objects.
generate_initial_candidates <- function(blocks,
                                        terms_df,
                                        seeds_df,
                                        config = list(),
                                        min_recall = 0.60,
                                        max_hits = 2000,
                                        max_gs_chars = 256L,
                                        max_candidates = 26L) {
  excl <- as.character(config$exclude_terms %||% character(0))
  req_terms <- tolower(stringr::str_squish(as.character(config$required_terms %||% character(0))))
  req_terms <- req_terms[nzchar(req_terms)]

  n_blocks <- length(blocks)
  if (n_blocks == 0L) {
    return(list())
  }

  block_tiers <- lapply(blocks, .block_tiers,
    terms_df = terms_df,
    max_terms = as.integer(config$max_terms_per_block %||% 5L)
  )

  cands <- list()
  add_cand <- function(q) {
    if (!is_valid_gs_query(q, max_gs_chars = max_gs_chars)) {
      return()
    }
    if (length(req_terms) > 0L) {
      q_terms <- unlist(lapply(q$blocks, `[[`, "terms"))
      if (!all(req_terms %in% q_terms)) {
        return()
      }
    }
    cands[[length(cands) + 1L]] <<- q
  }

  # 1-block queries only if only 1 block exists
  if (n_blocks == 1L) {
    for (t1 in block_tiers[[1L]]) {
      for (f1 in c("tiab", "title")) {
        add_cand(make_query(list(list(name = blocks[[1L]]$name, field = f1, terms = t1)), excl))
      }
    }
  }

  # 2-block combinations across all pairs from available blocks
  if (n_blocks >= 2L) {
    pairs <- utils::combn(seq_len(n_blocks), 2L, simplify = FALSE)
    field_2 <- list(c("tiab", "tiab"), c("title", "tiab"), c("tiab", "title"), c("title", "title"))
    for (pr in pairs) {
      i <- pr[1L]
      j <- pr[2L]
      for (t_i in block_tiers[[i]]) {
        for (t_j in block_tiers[[j]]) {
          for (fc in field_2) {
            add_cand(make_query(list(
              list(name = blocks[[i]]$name, field = fc[1L], terms = t_i),
              list(name = blocks[[j]]$name, field = fc[2L], terms = t_j)
            ), excl))
          }
        }
      }
    }
  }

  # 3-block combinations if >= 3 blocks exist
  if (n_blocks >= 3L) {
    triples <- utils::combn(seq_len(n_blocks), 3L, simplify = FALSE)
    field_3 <- list(c("tiab", "tiab", "tiab"), c("title", "tiab", "tiab"), c("tiab", "title", "tiab"))
    for (tr in triples) {
      i <- tr[1L]
      j <- tr[2L]
      m <- tr[3L]
      for (t_i in utils::head(block_tiers[[i]], 4L)) {
        for (t_j in utils::head(block_tiers[[j]], 4L)) {
          for (t_m in utils::head(block_tiers[[m]], 3L)) {
            for (fc in field_3) {
              add_cand(make_query(list(
                list(name = blocks[[i]]$name, field = fc[1L], terms = t_i),
                list(name = blocks[[j]]$name, field = fc[2L], terms = t_j),
                list(name = blocks[[m]]$name, field = fc[3L], terms = t_m)
              ), excl))
            }
          }
        }
      }
    }
  }

  # Domain-context qualified combinations: (primary_topic) AND (domain_context) AND (outcome / PICO blocks)
  domain_q <- attr(blocks, "domain_qualifier") %||% character(0)
  if (length(domain_q) > 0L && n_blocks >= 2L) {
    dom_blk <- list(name = "domain_context", field = "tiab", terms = domain_q)
    for (j in 2:n_blocks) {
      for (t_1 in block_tiers[[1L]]) {
        for (t_j in block_tiers[[j]]) {
          for (f1 in c("tiab", "title")) {
            add_cand(make_query(list(
              list(name = blocks[[1L]]$name, field = f1, terms = t_1),
              dom_blk,
              list(name = blocks[[j]]$name, field = "tiab", terms = t_j)
            ), excl))
          }
        }
      }
    }
    if (n_blocks >= 3L) {
      for (t_1 in utils::head(block_tiers[[1L]], 3L)) {
        for (t_2 in utils::head(block_tiers[[2L]], 3L)) {
          for (t_3 in utils::head(block_tiers[[3L]], 2L)) {
            add_cand(make_query(list(
              list(name = blocks[[1L]]$name, field = "tiab", terms = t_1),
              dom_blk,
              list(name = blocks[[2L]]$name, field = "tiab", terms = t_2),
              list(name = blocks[[3L]]$name, field = "tiab", terms = t_3)
            ), excl))
          }
        }
      }
    }
  }

  # Ensure at least a fallback candidate using required_terms if strict filter emptied pool
  if (length(cands) == 0L && length(req_terms) > 0L) {
    add_cand(make_query(lapply(seq_along(req_terms), function(i) {
      list(name = paste0("req_", i), field = "tiab", terms = req_terms[i])
    }), excl))
  }
  if (length(cands) == 0L) {
    return(list())
  }

  # Deduplicate by signature
  sigs <- vapply(cands, query_signature, character(1))
  cands <- cands[!duplicated(sigs)]

  # Score locally and select diverse frontier across recall buckets
  loc <- lapply(cands, estimate_local_query, terms_df = terms_df, seeds_df = seeds_df)
  found_cnt <- vapply(loc, `[[`, integer(1), "local_found")
  est_hits <- vapply(loc, `[[`, numeric(1), "est_hits")

  n_res <- if (any(seeds_df$in_pubmed %||% FALSE)) sum(seeds_df$in_pubmed) else nrow(seeds_df)
  min_cnt <- max(1L, floor((min_recall - 0.15) * n_res))
  avail_cnts <- sort(unique(found_cnt[found_cnt >= min_cnt]), decreasing = TRUE)
  if (length(avail_cnts) == 0L) avail_cnts <- sort(unique(found_cnt), decreasing = TRUE)

  is_all_tiab <- vapply(cands, function(q) {
    all(vapply(q$blocks, function(b) identical(b$field, "tiab"), logical(1)))
  }, logical(1))

  selected_idx <- integer(0)
  for (cnt in avail_cnts) {
    bucket <- which(found_cnt == cnt)
    over_cap <- if (is.finite(max_hits)) est_hits[bucket] > max_hits * 1.5 else rep(FALSE, length(bucket))
    bucket <- bucket[order(over_cap, est_hits[bucket])]
    tiab_in_bucket <- bucket[is_all_tiab[bucket]]
    picks <- unique(c(
      bucket[1L],
      if (length(tiab_in_bucket) > 0L) tiab_in_bucket[1L] else integer(0),
      bucket
    ))
    selected_idx <- c(selected_idx, utils::head(picks, 3L))
  }
  selected_idx <- unique(selected_idx)

  if (length(selected_idx) < max_candidates) {
    over_all <- if (is.finite(max_hits)) est_hits > max_hits * 1.5 else rep(FALSE, length(est_hits))
    rem <- setdiff(order(over_all, -found_cnt, est_hits), selected_idx)
    selected_idx <- c(selected_idx, utils::head(rem, max_candidates - length(selected_idx)))
  } else if (length(selected_idx) > max_candidates) {
    selected_idx <- utils::head(selected_idx, max_candidates)
  }

  cands[selected_idx]
}

#' Compute non-dominated Pareto frontier indicator
#'
#' Maximizes `combined_recall`, minimizes `pm_hits`. When `max_hits` is finite,
#' restricts the frontier to queries with `pm_hits <= max_hits` (falling back to
#' all valid queries only if none satisfy the cap).
#'
#' @param eval_df Data frame / tibble of evaluated queries.
#' @param max_hits Optional maximum PubMed hit count cap.
#' @return Logical vector of length `nrow(eval_df)`.
compute_pareto_frontier <- function(eval_df, max_hits = Inf) {
  n <- nrow(eval_df)
  if (n == 0L) {
    return(logical(0))
  }
  is_p <- rep(FALSE, n)

  valid <- !is.na(eval_df$combined_recall) & !is.na(eval_df$pm_hits)
  if (any(valid & eval_df$combined_recall > 0)) {
    valid <- valid & eval_df$combined_recall > 0
  }
  if (is.finite(max_hits)) {
    within_idx <- which(valid & eval_df$pm_hits <= max_hits)
    v_idx <- if (length(within_idx) > 0L) within_idx else which(valid)
  } else {
    v_idx <- which(valid)
  }
  if (length(v_idx) == 0L) {
    return(is_p)
  }

  for (i in v_idx) {
    r_i <- eval_df$combined_recall[i]
    h_i <- eval_df$pm_hits[i]
    dominated <- FALSE
    for (j in v_idx) {
      if (i == j) next
      r_j <- eval_df$combined_recall[j]
      h_j <- eval_df$pm_hits[j]
      if (r_j >= r_i && h_j <= h_i && (r_j > r_i || h_j < h_i)) {
        dominated <- TRUE
        break
      }
    }
    if (!dominated) is_p[i] <- TRUE
  }

  # Break exact ties on (combined_recall, pm_hits)
  p_idx <- which(is_p)
  if (length(p_idx) > 1L) {
    tie_key <- paste(round(eval_df$combined_recall[p_idx], 8), eval_df$pm_hits[p_idx], sep = "|")
    groups <- split(p_idx, tie_key)
    for (grp in groups) {
      if (length(grp) <= 1L) next
      ord <- order(eval_df$gs_chars[grp], eval_df$n_terms[grp])
      is_p[grp[-ord[1L]]] <- FALSE
    }
  }
  is_p
}

# Generate Stage 2 greedy neighbors from current Pareto queries
.generate_neighbors <- function(pareto_queries,
                                pareto_missed_idx,
                                blocks,
                                terms_df,
                                seeds_df,
                                config,
                                seen_sigs,
                                max_gs_chars = 256L,
                                setcover_ex = character(0)) {
  req_terms <- tolower(stringr::str_squish(as.character(config$required_terms %||% character(0))))
  req_terms <- req_terms[nzchar(req_terms)]
  max_terms_blk <- as.integer(config$max_terms_per_block %||% 5L)

  neighbors <- list()
  add_nbr <- function(q) {
    if (!is_valid_gs_query(q, max_gs_chars = max_gs_chars)) {
      return()
    }
    if (length(req_terms) > 0L) {
      q_terms <- unlist(lapply(q$blocks, `[[`, "terms"))
      if (!all(req_terms %in% q_terms)) {
        return()
      }
    }
    sig <- query_signature(q)
    if (!(sig %in% seen_sigs)) {
      neighbors[[sig]] <<- q
    }
  }

  for (qi in seq_along(pareto_queries)) {
    q <- pareto_queries[[qi]]
    missed <- pareto_missed_idx[[qi]]
    nb <- length(q$blocks)

    for (bi in seq_len(nb)) {
      blk <- q$blocks[[bi]]

      # 1. Prune 1 term if block has >= 2 terms
      if (length(blk$terms) >= 2L) {
        for (ti in seq_along(blk$terms)) {
          new_blks <- q$blocks
          new_blks[[bi]]$terms <- blk$terms[-ti]
          add_nbr(make_query(new_blks, q$exclude_terms))
        }
      }

      # 2. Tighten field (tiab -> title) or Loosen field (title -> tiab)
      new_blks_f <- q$blocks
      new_blks_f[[bi]]$field <- if (blk$field == "tiab") "title" else "tiab"
      add_nbr(make_query(new_blks_f, q$exclude_terms))

      # 3. Narrowing swap: replace a broad term with a narrower lexical variant
      for (ti in seq_along(blk$terms)) {
        cur_t <- blk$terms[ti]
        cur_row <- match(cur_t, terms_df$term)
        cur_bg <- if (!is.na(cur_row)) terms_df$bg_hits[cur_row] else Inf
        lex_cands <- terms_df$term[
          !(terms_df$term %in% blk$terms) &
            terms_df$bg_hits < cur_bg &
            vapply(terms_df$term, shares_stem, logical(1), t2 = cur_t)
        ]
        for (lc in utils::head(lex_cands, 2L)) {
          new_blks_s <- q$blocks
          new_blks_s[[bi]]$terms[ti] <- lc
          add_nbr(make_query(new_blks_s, q$exclude_terms))
        }
      }

      # 4. Rescue missed seeds: add 1 high-lift term from the SAME semantic family covering missed seeds
      if (length(missed) > 0L && length(blk$terms) < max_terms_blk) {
        other_terms <- unlist(lapply(q$blocks[-bi], `[[`, "terms"))
        b1_coupled <- attr(blocks, "block1_coupled") %||% character(0)
        term_fams <- attr(blocks, "term_families") %||% assign_term_families(terms_df$term, b1_coupled)
        blk_fam <- if (blk$name %in% c("primary_topic", "outcome_severity", "population_host", "epi_methods")) {
          blk$name
        } else if (identical(blk$name, "outcome_or_epi")) {
          "outcome_severity"
        } else {
          unname(term_fams[blk$terms[1L]]) %||% "other"
        }
        cand_rows <- which(
          !(terms_df$term %in% blk$terms) &
            !(terms_df$term %in% other_terms) &
            terms_df$lift >= 5
        )
        cand_rows <- Filter(function(r) {
          trm <- terms_df$term[r]
          tf <- unname(term_fams[trm]) %||% assign_term_family(trm, b1_coupled)
          if (blk_fam == "other") {
            !(tf %in% c("primary_topic", "outcome_severity", "population_host", "epi_methods")) &&
              !(trm %in% b1_coupled)
          } else {
            identical(tf, blk_fam) || any(vapply(blk$terms, shares_stem, logical(1), t2 = trm))
          }
        }, cand_rows)
        if (bi >= 2L && length(b1_coupled) > 0L) {
          cand_rows <- cand_rows[!(terms_df$term[cand_rows] %in% b1_coupled)]
        }
        if (length(other_terms) > 0L) {
          cand_rows <- Filter(function(r) {
            !any(vapply(other_terms, shares_stem, logical(1), t2 = terms_df$term[r]))
          }, cand_rows)
        }
        best_r <- integer(0)
        if (length(cand_rows) > 0L) {
          rescue_gain <- vapply(cand_rows, function(r) {
            cov_m <- sum(terms_df$cov_tiab[[r]][missed])
            if (cov_m == 0L) {
              return(-1)
            }
            cov_m * log10(pmax(terms_df$lift[r], 1.1)) / log10(terms_df$bg_hits[r] + 10)
          }, numeric(1))
          best_r <- cand_rows[order(rescue_gain, decreasing = TRUE)]
          best_r <- best_r[rescue_gain[order(rescue_gain, decreasing = TRUE)] > 0]
          for (ri in utils::head(best_r, 2L)) {
            new_blks_r <- q$blocks
            new_blks_r[[bi]]$terms <- c(blk$terms, terms_df$term[ri])
            add_nbr(make_query(new_blks_r, q$exclude_terms))
          }
        }

        # 4b. If outcome block still misses seeds not covered by any outcome_severity term,
        # allow 1 high-lift epi_methods term (e.g. "zoonotic") + optional remaining outcome rescue in "outcome_or_epi"
        if (blk$name %in% c("outcome_severity", "outcome_or_epi")) {
          covered_by_same <- if (length(best_r) > 0L) {
            missed[Reduce(`|`, lapply(best_r, function(r) terms_df$cov_tiab[[r]][missed]))]
          } else {
            integer(0)
          }
          still_missed <- setdiff(missed, covered_by_same)
          if (length(still_missed) > 0L) {
            epi_rows <- which(
              !(terms_df$term %in% blk$terms) &
                !(terms_df$term %in% other_terms) &
                !(terms_df$term %in% b1_coupled) &
                unname(term_fams[terms_df$term]) == "epi_methods" &
                terms_df$lift >= 25 &
                terms_df$bg_hits <= 100000
            )
            if (length(epi_rows) > 0L) {
              epi_gain <- vapply(epi_rows, function(r) {
                cov_m <- sum(terms_df$cov_tiab[[r]][still_missed])
                if (cov_m == 0L) return(-1)
                cov_m * log10(pmax(terms_df$lift[r], 1.1)) / log10(terms_df$bg_hits[r] + 10)
              }, numeric(1))
              best_epi <- epi_rows[order(epi_gain, decreasing = TRUE)]
              best_epi <- best_epi[epi_gain[order(epi_gain, decreasing = TRUE)] > 0]
              for (ri in utils::head(best_epi, 1L)) {
                new_blks_e <- q$blocks
                new_blks_e[[bi]]$name <- "outcome_or_epi"
                new_blks_e[[bi]]$terms <- c(blk$terms, terms_df$term[ri])
                add_nbr(make_query(new_blks_e, q$exclude_terms))
                # Also generate 2-step full rescue (epi term + best same-family outcome term) to reach 100% in 1 round
                if (length(best_r) > 0L && (length(blk$terms) + 2L) <= max_terms_blk) {
                  new_blks_e2 <- q$blocks
                  new_blks_e2[[bi]]$name <- "outcome_or_epi"
                  new_blks_e2[[bi]]$terms <- c(blk$terms, terms_df$term[ri], terms_df$term[best_r[1L]])
                  add_nbr(make_query(new_blks_e2, q$exclude_terms))
                }
              }
            }
          }
        }
      }
    }

    # 5. Add a 3rd narrowing block (or domain_context qualifier) if not already present
    used_names <- vapply(q$blocks, `[[`, character(1), "name")
    domain_q <- attr(blocks, "domain_qualifier") %||% character(0)
    if (length(domain_q) > 0L && !("domain_context" %in% used_names) && nb >= 1L) {
      dom_blk <- list(name = "domain_context", field = "tiab", terms = domain_q)
      new_blks_d <- c(q$blocks[1L], list(dom_blk), if (nb >= 2L) q$blocks[2:nb] else list())
      add_nbr(make_query(new_blks_d, q$exclude_terms))
    }
    if (nb == 2L && length(blocks) >= 3L) {
      avail_blks <- Filter(function(b) !(b$name %in% used_names), blocks)
      for (ab in avail_blks) {
        for (n_t in seq_len(min(2L, length(ab$terms)))) {
          new_blks_a <- c(q$blocks, list(list(name = ab$name, field = "tiab", terms = ab$terms[seq_len(n_t)])))
          add_nbr(make_query(new_blks_a, q$exclude_terms))
        }
      }
    }

    # 6. Extend exclusion set on high-recall queries using set-cover exclusions
    if (length(setcover_ex) > 0L) {
      add_ex <- setdiff(setcover_ex, q$exclude_terms)
      if (length(add_ex) >= 2L) {
        add_nbr(make_query(q$blocks, unique(c(q$exclude_terms, utils::head(add_ex, 3L)))))
      }
    }
  }

  unname(neighbors)
}

#' Run Stage 1 + Stage 2 Pareto query optimization
#'
#' @param seeds_df Enriched seeds tibble.
#' @param terms_df Scored candidate terms tibble.
#' @param blocks Concept blocks list.
#' @param config Optional config list.
#' @param min_recall Minimum desired seed recall (0..1).
#' @param max_hits Optional maximum PubMed hit count threshold.
#' @param max_evals Maximum live API query evaluations.
#' @param max_gs_chars Maximum Google Scholar character limit (default 256).
#' @param discover_exclusions Logical; whether to mine non-seed PubMed results for exclusion hints.
#' @param email Optional email for NCBI polite pool.
#' @param ncbi_key Optional NCBI API key.
#' @param cache_dir Cache directory path.
#' @param use_cache Logical; whether to read/write cache.
#' @param eval_pubmed_fn Function `(pm_query, seed_pmids, ...)` returning `pm_hits`, `pm_found`, `pm_recall`, `ok`.
#' @param pm_esearch_fn Optional dependency injection for [pubmed_esearch()].
#' @param pm_efetch_fn Optional dependency injection for [pubmed_efetch()].
#' @param progress_cb Optional progress callback function.
#' @return Optimization result list containing `eval_df`, `pareto_df`, `exclusion_hints`, `recommended_sig`, `seeds_df`, `blocks`.
optimize_queries <- function(seeds_df,
                             terms_df,
                             blocks,
                             config = list(),
                             min_recall = 0.60,
                             max_hits = 2000,
                             max_evals = 60L,
                             max_gs_chars = 256L,
                             discover_exclusions = TRUE,
                             email = NULL,
                             ncbi_key = NULL,
                             cache_dir = ".cache",
                             use_cache = TRUE,
                             eval_pubmed_fn = pubmed_eval_query,
                             pm_esearch_fn = pubmed_esearch,
                             pm_efetch_fn = pubmed_efetch,
                             progress_cb = NULL) {
  max_evals <- as.integer(max_evals)
  if (is.null(max_hits) || length(max_hits) == 0L || is.na(max_hits) || max_hits <= 0) {
    max_hits <- Inf
  }
  seed_pmids <- seeds_df$pmid[seeds_df$in_pubmed]
  resolvable_idx <- which(seeds_df$in_pubmed)
  n_resolvable <- length(resolvable_idx)

  stage1_budget <- max(5L, floor(max_evals * 0.45))
  init_cands <- generate_initial_candidates(
    blocks = blocks, terms_df = terms_df, seeds_df = seeds_df,
    config = config, min_recall = min_recall, max_hits = max_hits,
    max_gs_chars = max_gs_chars, max_candidates = stage1_budget
  )
  if (length(init_cands) == 0L) {
    stop("No valid candidate queries could be generated under the constraints (max_gs_chars=", max_gs_chars, ").", call. = FALSE)
  }

  eval_records <- list()
  seen_sigs <- character(0)

  evaluate_one <- function(q, stage_label) {
    if (is.function(progress_cb)) {
      progress_cb(sprintf("Evaluating query %d/%d (%s)", length(eval_records) + 1L, max_evals, stage_label))
    }
    sig <- query_signature(q)
    pm_q <- to_pubmed_query(q)
    gs_q <- to_google_scholar_query(q)

    pm_res <- eval_pubmed_fn(
      pm_q, seed_pmids,
      email = email, api_key = ncbi_key, cache_dir = cache_dir, use_cache = use_cache
    )

    pm_hit_idx <- if (length(pm_res$pm_found) > 0L) which(seeds_df$pmid %in% pm_res$pm_found) else integer(0)

    loc <- estimate_local_query(q, terms_df, seeds_df)
    if (n_resolvable > 0L) {
      comb_in_resolvable <- intersect(pm_hit_idx, resolvable_idx)
      comb_recall <- length(comb_in_resolvable) / n_resolvable
      seeds_found_n <- length(comb_in_resolvable)
      missed_idx <- setdiff(resolvable_idx, comb_in_resolvable)
    } else {
      comb_recall <- loc$local_recall
      seeds_found_n <- loc$local_found
      missed_idx <- which(!loc$local_cov)
    }

    missed_keys <- paste(seeds_df$citekey[missed_idx], collapse = "; ")
    n_terms <- sum(vapply(q$blocks, function(b) length(b$terms), integer(1)))

    list(
      query_obj = q,
      signature = sig,
      stage = stage_label,
      is_custom = FALSE,
      custom_label = "",
      pubmed_query = pm_q,
      google_scholar_query = gs_q,
      gs_chars = nchar(gs_q, type = "chars"),
      n_blocks = length(q$blocks),
      n_terms = n_terms,
      pm_hits = as.integer(pm_res$pm_hits),
      pm_recall = as.numeric(pm_res$pm_recall),
      combined_recall = as.numeric(comb_recall),
      seeds_found = as.integer(seeds_found_n),
      seeds_resolvable = as.integer(if (n_resolvable > 0L) n_resolvable else nrow(seeds_df)),
      found_idx = list(if (n_resolvable > 0L) intersect(pm_hit_idx, resolvable_idx) else which(loc$local_cov)),
      missed_idx = list(missed_idx),
      missed_citekeys = missed_keys,
      pubmed_url = pubmed_url(pm_q),
      google_scholar_url = google_scholar_url(gs_q)
    )
  }

  # Stage 1 evaluation
  cli::cli_inform(c("i" = "Stage 1: Evaluating {length(init_cands)} initial candidate queries..."))
  for (q in init_cands) {
    if (length(eval_records) >= max_evals) break
    sig <- query_signature(q)
    if (sig %in% seen_sigs) next
    seen_sigs <- c(seen_sigs, sig)
    eval_records[[length(eval_records) + 1L]] <- evaluate_one(q, "stage1")
  }

  to_eval_df <- function(recs) {
    df <- tibble::tibble(
      eval_id = seq_along(recs),
      signature = vapply(recs, `[[`, character(1), "signature"),
      stage = vapply(recs, `[[`, character(1), "stage"),
      is_custom = vapply(recs, `[[`, logical(1), "is_custom"),
      custom_label = vapply(recs, `[[`, character(1), "custom_label"),
      pubmed_query = vapply(recs, `[[`, character(1), "pubmed_query"),
      google_scholar_query = vapply(recs, `[[`, character(1), "google_scholar_query"),
      gs_chars = vapply(recs, `[[`, integer(1), "gs_chars"),
      n_blocks = vapply(recs, `[[`, integer(1), "n_blocks"),
      n_terms = vapply(recs, `[[`, integer(1), "n_terms"),
      pm_hits = vapply(recs, `[[`, integer(1), "pm_hits"),
      pm_recall = vapply(recs, `[[`, numeric(1), "pm_recall"),
      combined_recall = vapply(recs, `[[`, numeric(1), "combined_recall"),
      seeds_found = vapply(recs, `[[`, integer(1), "seeds_found"),
      seeds_resolvable = vapply(recs, `[[`, integer(1), "seeds_resolvable"),
      missed_citekeys = vapply(recs, `[[`, character(1), "missed_citekeys"),
      pubmed_url = vapply(recs, `[[`, character(1), "pubmed_url"),
      google_scholar_url = vapply(recs, `[[`, character(1), "google_scholar_url"),
      found_idx = lapply(recs, function(r) r$found_idx[[1L]]),
      missed_idx = lapply(recs, function(r) r$missed_idx[[1L]]),
      query_obj = lapply(recs, `[[`, "query_obj")
    )
    df$within_max_hits <- is.na(df$pm_hits) | (df$pm_hits <= max_hits)
    df$is_pareto <- compute_pareto_frontier(df, max_hits = max_hits)
    df
  }

  # Stage 2a: Discover concept-grouped exclusion hints & apply set-cover / concept exclusion bundles
  exclusion_hints_df <- .empty_exclusion_hints()
  if (isTRUE(discover_exclusions) && length(eval_records) < max_evals && exists("discover_exclusion_hints", mode = "function")) {
    s1_df <- to_eval_df(eval_records)
    hi_bound <- if (is.finite(max_hits)) max(max_hits * 3, 6000) else 8000
    max_s1_rec <- max(s1_df$combined_recall, na.rm = TRUE)
    s1_under_cap <- s1_df |>
      dplyr::filter(
        !is.na(.data$pm_hits) &
          .data$pm_hits >= 100L &
          .data$pm_hits <= max_hits &
          .data$combined_recall >= max(min_recall, max_s1_rec - 0.10)
      ) |>
      dplyr::arrange(dplyr::desc(.data$combined_recall), .data$pm_hits)
    s1_focused <- s1_df |>
      dplyr::filter(!is.na(.data$pm_hits) & .data$pm_hits >= 100L & .data$pm_hits <= hi_bound) |>
      dplyr::arrange(dplyr::desc(.data$combined_recall), .data$pm_hits)
    s1_noisy <- s1_df |>
      dplyr::filter(!is.na(.data$pm_hits) & .data$pm_hits >= 100L) |>
      dplyr::arrange(dplyr::desc(.data$combined_recall), .data$pm_hits)
    base_row <- if (nrow(s1_under_cap) > 0L) {
      s1_under_cap[1L, ]
    } else if (nrow(s1_focused) > 0L) {
      s1_focused[1L, ]
    } else if (nrow(s1_noisy) > 0L) {
      s1_noisy[1L, ]
    } else {
      s1_df[which.max(s1_df$combined_recall), ]
    }

    if (nrow(base_row) == 1L && !is.na(base_row$pm_hits) && base_row$pm_hits > 30L) {
      cli::cli_inform(c("i" = "Stage 2a: Mining non-seed PubMed results for concept-grouped exclusion hints..."))
      exclusion_hints_df <- tryCatch(
        discover_exclusion_hints(
          base_query = base_row$query_obj[[1L]],
          seeds_df = seeds_df,
          base_pm_hits = base_row$pm_hits,
          sample_n = 160L,
          top_verify = 12L,
          max_hints = 35L,
          email = email,
          ncbi_key = ncbi_key,
          cache_dir = cache_dir,
          use_cache = use_cache,
          pm_esearch_fn = pm_esearch_fn,
          pm_efetch_fn = pm_efetch_fn,
          progress_cb = progress_cb
        ),
        error = function(e) .empty_exclusion_hints()
      )

      if (nrow(exclusion_hints_df) > 0L) {
        safe_ex <- exclusion_hints_df$exclude_spec[exclusion_hints_df$seeds_lost_n == 0L]
        setcover_ex <- attr(exclusion_hints_df, "setcover_ex") %||% safe_ex
        code_bundle <- attr(exclusion_hints_df, "code_bundle") %||% character(0)
        excl_blocks <- attr(exclusion_hints_df, "exclusion_blocks") %||% list()
        non_code_safe <- setdiff(setcover_ex, attr(exclusion_hints_df, "all_code_ex") %||% code_bundle)
        lowrisk_ex <- exclusion_hints_df$exclude_spec[
          exclusion_hints_df$seeds_lost_n <= 1L & exclusion_hints_df$pct_hits_reduced >= 15
        ]

        # Select best Stage 1 queries across the top 4 distinct high-recall levels
        s1_ranked <- s1_df |>
          dplyr::filter(!is.na(.data$pm_hits) & .data$combined_recall >= min_recall) |>
          dplyr::arrange(dplyr::desc(.data$seeds_found), .data$pm_hits)
        if (nrow(s1_ranked) == 0L) {
          s1_ranked <- s1_df |> dplyr::arrange(dplyr::desc(.data$combined_recall), .data$pm_hits)
        }
        top_levels <- utils::head(unique(s1_ranked$seeds_found), 4L)
        top_s1_idx <- unlist(lapply(top_levels, function(lv) {
          rows_lv <- which(s1_ranked$seeds_found == lv)
          utils::head(rows_lv, 1L)
        }))
        top_s1 <- s1_ranked[unique(c(top_s1_idx, seq_len(min(3L, nrow(s1_ranked))))), ]
        top_s1 <- utils::head(top_s1, 4L)

        excl_cands <- list()
        add_ex_cand <- function(bq, ex_vec) {
          if (length(ex_vec) == 0L) return()
          q_new <- make_query(bq$blocks, unique(c(bq$exclude_terms, ex_vec)))
          if (is_valid_gs_query(q_new, max_gs_chars)) excl_cands[[length(excl_cands) + 1L]] <<- q_new
        }

        for (ri in seq_len(nrow(top_s1))) {
          bq <- top_s1$query_obj[[ri]]
          # 1. Pure sibling subtype/code exclusion bundle (e.g. NOT (h7n9 OR h1n1 OR h1n1pdm09))
          if (length(code_bundle) > 0L) {
            add_ex_cand(bq, code_bundle)
            if (length(non_code_safe) >= 2L) add_ex_cand(bq, c(code_bundle, utils::head(non_code_safe, 2L)))
            if (length(non_code_safe) >= 5L) add_ex_cand(bq, c(code_bundle, utils::head(non_code_safe, 5L)))
          }
          # 2. Concept-block exclusion bundles
          if (length(excl_blocks) >= 1L) {
            add_ex_cand(bq, excl_blocks[[1L]]$terms)
          }
          if (length(excl_blocks) >= 2L) {
            add_ex_cand(bq, unique(c(excl_blocks[[1L]]$terms, excl_blocks[[2L]]$terms)))
          }
          # 3. Progressive global greedy set-cover exclusion bundles
          for (k_ex in c(2L, 4L, 7L, 10L)) {
            if (length(setcover_ex) >= k_ex) {
              add_ex_cand(bq, setcover_ex[seq_len(k_ex)])
            }
          }
          if (length(lowrisk_ex) >= 2L) {
            add_ex_cand(bq, utils::head(lowrisk_ex, 3L))
          }
        }

        # Leave at least ~8 evals for Stage 2b greedy refinement
        s2a_limit <- max(0L, min(length(excl_cands), max_evals - length(eval_records) - 8L))
        if (s2a_limit > 0L) {
          cli::cli_inform(c("i" = "Stage 2a: Evaluating {s2a_limit} exclusion-refined candidate queries..."))
          eval_n_2a <- 0L
          for (q_ex in excl_cands) {
            if (eval_n_2a >= s2a_limit || length(eval_records) >= max_evals) break
            sig <- query_signature(q_ex)
            if (sig %in% seen_sigs) next
            seen_sigs <- c(seen_sigs, sig)
            eval_records[[length(eval_records) + 1L]] <- evaluate_one(q_ex, "stage2_excl")
            eval_n_2a <- eval_n_2a + 1L
          }
        }
      }
    }
  }

  # Stage 2b greedy refinement loop
  round_num <- 0L
  setcover_ex <- attr(exclusion_hints_df, "setcover_ex") %||% character(0)
  while (length(eval_records) < max_evals) {
    round_num <- round_num + 1L
    cur_df <- to_eval_df(eval_records)
    p_rows <- which(cur_df$is_pareto)
    over_hi <- which(!cur_df$within_max_hits & cur_df$combined_recall >= min_recall)
    if (length(over_hi) > 0L) {
      over_hi <- utils::head(over_hi[order(-cur_df$combined_recall[over_hi], cur_df$pm_hits[over_hi])], 2L)
    }
    seed_rows <- unique(c(p_rows, over_hi))
    if (length(seed_rows) == 0L) break

    prev_pareto_sigs <- cur_df$signature[p_rows]
    nbrs <- .generate_neighbors(
      pareto_queries = cur_df$query_obj[seed_rows],
      pareto_missed_idx = cur_df$missed_idx[seed_rows],
      blocks = blocks,
      terms_df = terms_df,
      seeds_df = seeds_df,
      config = config,
      seen_sigs = seen_sigs,
      max_gs_chars = max_gs_chars,
      setcover_ex = setcover_ex
    )
    if (length(nbrs) == 0L) break

    loc_n <- lapply(nbrs, estimate_local_query, terms_df = terms_df, seeds_df = seeds_df)
    n_found <- vapply(loc_n, `[[`, integer(1), "local_found")
    n_hits <- vapply(loc_n, `[[`, numeric(1), "est_hits")

    rem_budget <- max_evals - length(eval_records)
    batch_size <- min(8L, rem_budget)
    pick_idx <- integer(0)
    for (cnt in sort(unique(n_found), decreasing = TRUE)) {
      b_idx <- which(n_found == cnt)
      over_cap <- if (is.finite(max_hits)) n_hits[b_idx] > max_hits * 1.25 else rep(FALSE, length(b_idx))
      b_idx <- b_idx[order(over_cap, n_hits[b_idx])]
      pick_idx <- c(pick_idx, utils::head(b_idx, 2L))
    }
    pick_idx <- utils::head(unique(pick_idx), batch_size)
    if (length(pick_idx) == 0L) break

    cli::cli_inform(c("i" = "Stage 2b (round {round_num}): Evaluating {length(pick_idx)} refined neighbor queries..."))
    for (idx in pick_idx) {
      if (length(eval_records) >= max_evals) break
      q_n <- nbrs[[idx]]
      sig <- query_signature(q_n)
      if (sig %in% seen_sigs) next
      seen_sigs <- c(seen_sigs, sig)
      eval_records[[length(eval_records) + 1L]] <- evaluate_one(q_n, paste0("stage2_r", round_num))
    }

    new_df <- to_eval_df(eval_records)
    new_pareto_sigs <- new_df$signature[new_df$is_pareto]
    if (length(setdiff(new_pareto_sigs, prev_pareto_sigs)) == 0L) {
      break
    }
  }

  eval_df <- to_eval_df(eval_records)
  pareto_df <- eval_df |>
    dplyr::filter(.data$is_pareto) |>
    dplyr::arrange(dplyr::desc(.data$combined_recall), .data$pm_hits, .data$gs_chars) |>
    dplyr::mutate(rank = dplyr::row_number())

  meets_both <- pareto_df |>
    dplyr::filter(.data$combined_recall >= min_recall & .data$pm_hits <= max_hits)

  if (nrow(meets_both) > 0L) {
    rec_sig <- meets_both$signature[1L]
    rec_status <- "meets_constraints"
  } else {
    under_cap <- pareto_df |> dplyr::filter(.data$pm_hits <= max_hits)
    if (nrow(under_cap) > 0L) {
      rec_sig <- under_cap$signature[1L]
      rec_status <- "below_min_recall"
    } else {
      meets_rec <- pareto_df |> dplyr::filter(.data$combined_recall >= min_recall)
      if (nrow(meets_rec) > 0L) {
        rec_sig <- meets_rec$signature[which.min(meets_rec$pm_hits)]
        rec_status <- "exceeds_max_hits"
      } else {
        rec_sig <- pareto_df$signature[1L]
        rec_status <- "exceeds_max_hits"
      }
    }
  }

  pareto_df <- pareto_df |>
    dplyr::mutate(is_recommended = .data$signature == rec_sig)

  out_res <- list(
    eval_df = eval_df,
    pareto_df = pareto_df,
    exclusion_hints = exclusion_hints_df,
    recommended_sig = rec_sig,
    recommended_status = rec_status,
    min_recall = min_recall,
    max_hits = max_hits,
    max_gs_chars = max_gs_chars,
    seeds_df = seeds_df,
    terms_df = terms_df,
    blocks = blocks
  )

  # Evaluate any user-supplied custom benchmark queries (outside max_evals budget)
  cq_list <- .normalize_custom_queries_config(config$custom_queries)
  if (length(cq_list) > 0L) {
    cli::cli_inform(c("i" = "Evaluating {length(cq_list)} custom benchmark quer{?y/ies}..."))
    for (cq in cq_list) {
      c_row <- evaluate_custom_query(
        pm_query = cq$pubmed_query,
        seeds_df = seeds_df,
        label = cq$label,
        gs_query = cq$google_scholar_query,
        email = email,
        ncbi_key = ncbi_key,
        cache_dir = cache_dir,
        use_cache = use_cache,
        eval_pubmed_fn = eval_pubmed_fn
      )
      out_res <- add_custom_query_to_result(out_res, c_row)
    }
  }

  if (is.finite(max_hits) && nrow(out_res$pareto_df) > 0L && all(out_res$pareto_df$pm_hits > max_hits)) {
    warning("No evaluated query satisfied max_hits <= ", max_hits, ".", call. = FALSE)
  }

  out_res
}

# Normalize config$custom_queries (character vector or list of named lists)
.normalize_custom_queries_config <- function(cq) {
  if (is.null(cq) || length(cq) == 0L) {
    return(list())
  }
  out <- list()
  if (is.character(cq)) {
    for (i in seq_along(cq)) {
      pq <- stringr::str_squish(cq[i])
      if (nzchar(pq)) {
        out[[length(out) + 1L]] <- list(
          label = paste0("Custom ", i),
          pubmed_query = pq,
          google_scholar_query = NULL
        )
      }
    }
  } else if (is.list(cq)) {
    for (i in seq_along(cq)) {
      item <- cq[[i]]
      if (is.character(item) && length(item) == 1L) {
        pq <- stringr::str_squish(item)
        if (nzchar(pq)) {
          out[[length(out) + 1L]] <- list(
            label = names(cq)[i] %||% paste0("Custom ", i),
            pubmed_query = pq,
            google_scholar_query = NULL
          )
        }
      } else if (is.list(item)) {
        pq <- stringr::str_squish(as.character(item$pubmed_query %||% item$pubmed %||% item$query %||% ""))
        if (nzchar(pq)) {
          lbl <- stringr::str_squish(as.character(item$label %||% item$name %||% paste0("Custom ", i)))
          if (!nzchar(lbl)) lbl <- paste0("Custom ", i)
          gsq <- item$google_scholar_query %||% item$gs_query %||% NULL
          out[[length(out) + 1L]] <- list(
            label = lbl,
            pubmed_query = pq,
            google_scholar_query = gsq
          )
        }
      }
    }
  }
  out
}

#' Evaluate a single user-supplied PubMed query against enriched seeds
#'
#' @param pm_query Raw PubMed query string.
#' @param seeds_df Enriched seeds tibble.
#' @param label Display label for the custom query.
#' @param gs_query Optional Google Scholar query string (auto-approximated if NULL/blank).
#' @param email Optional email for NCBI polite pool.
#' @param ncbi_key Optional NCBI API key.
#' @param cache_dir Cache directory path.
#' @param use_cache Logical; whether to read/write cache.
#' @param eval_pubmed_fn Function `(pm_query, seed_pmids, ...)` returning `pm_hits`, `pm_found`, `pm_recall`, `ok`.
#' @return 1-row tibble matching `eval_df` schema.
evaluate_custom_query <- function(pm_query,
                                  seeds_df,
                                  label = "Custom",
                                  gs_query = NULL,
                                  email = NULL,
                                  ncbi_key = NULL,
                                  cache_dir = ".cache",
                                  use_cache = TRUE,
                                  eval_pubmed_fn = pubmed_eval_query) {
  pm_q <- stringr::str_squish(as.character(pm_query %||% ""))
  if (!nzchar(pm_q)) stop("Custom PubMed query cannot be empty.", call. = FALSE)
  lbl <- stringr::str_squish(as.character(label %||% ""))
  if (!nzchar(lbl)) lbl <- "Custom"

  gs_q <- if (!is.null(gs_query) && nzchar(stringr::str_squish(as.character(gs_query)))) {
    stringr::str_squish(as.character(gs_query))
  } else {
    pubmed_to_google_scholar_approx(pm_q)
  }

  seed_pmids <- seeds_df$pmid[seeds_df$in_pubmed]
  resolvable_idx <- which(seeds_df$in_pubmed)
  n_resolvable <- length(resolvable_idx)

  pm_res <- eval_pubmed_fn(
    pm_q, seed_pmids,
    email = email, api_key = ncbi_key, cache_dir = cache_dir, use_cache = use_cache
  )
  if (!isTRUE(pm_res$ok)) {
    stop("PubMed evaluation failed for custom query.", call. = FALSE)
  }

  pm_hit_idx <- if (length(pm_res$pm_found) > 0L) which(seeds_df$pmid %in% pm_res$pm_found) else integer(0)

  if (n_resolvable > 0L) {
    comb_in_res <- intersect(pm_hit_idx, resolvable_idx)
    comb_recall <- length(comb_in_res) / n_resolvable
    seeds_found_n <- length(comb_in_res)
    missed_idx <- setdiff(resolvable_idx, comb_in_res)
    found_final <- comb_in_res
  } else {
    comb_recall <- 0
    seeds_found_n <- 0L
    missed_idx <- seq_len(nrow(seeds_df))
    found_final <- integer(0)
  }

  sig <- paste0("custom::", digest::digest(pm_q, algo = "md5"))
  tibble::tibble(
    eval_id = 1L,
    signature = sig,
    stage = "custom",
    is_custom = TRUE,
    custom_label = lbl,
    pubmed_query = pm_q,
    google_scholar_query = gs_q,
    gs_chars = nchar(gs_q, type = "chars"),
    n_blocks = NA_integer_,
    n_terms = NA_integer_,
    pm_hits = as.integer(pm_res$pm_hits),
    pm_recall = as.numeric(pm_res$pm_recall),
    combined_recall = as.numeric(comb_recall),
    seeds_found = as.integer(seeds_found_n),
    seeds_resolvable = as.integer(if (n_resolvable > 0L) n_resolvable else nrow(seeds_df)),
    missed_citekeys = paste(seeds_df$citekey[missed_idx], collapse = "; "),
    pubmed_url = pubmed_url(pm_q),
    google_scholar_url = google_scholar_url(gs_q),
    found_idx = list(found_final),
    missed_idx = list(missed_idx),
    query_obj = list(NULL),
    within_max_hits = TRUE,
    is_pareto = FALSE
  )
}

#' Merge an evaluated custom query into an existing or minimal optimization result
#'
#' @param result Optimization result list from `optimize_queries()` (or minimal shell).
#' @param custom_row 1-row tibble returned by `evaluate_custom_query()`.
#' @return Updated optimization result list with `eval_df` and `pareto_df` recomputed.
add_custom_query_to_result <- function(result, custom_row) {
  eval_df <- result$eval_df %||% tibble::tibble()
  max_hits <- result$max_hits %||% Inf
  min_recall <- result$min_recall %||% 0.60

  if (nrow(eval_df) == 0L) {
    custom_row$eval_id <- 1L
    new_eval <- custom_row
  } else {
    if (!("is_custom" %in% names(eval_df))) eval_df$is_custom <- FALSE
    if (!("custom_label" %in% names(eval_df))) eval_df$custom_label <- ""
    dup_idx <- which(eval_df$signature == custom_row$signature)
    if (length(dup_idx) > 0L) {
      custom_row$eval_id <- eval_df$eval_id[dup_idx[1L]]
      eval_df[dup_idx[1L], ] <- custom_row
      new_eval <- eval_df
    } else {
      custom_row$eval_id <- max(eval_df$eval_id, na.rm = TRUE) + 1L
      new_eval <- dplyr::bind_rows(eval_df, custom_row)
    }
  }

  new_eval$within_max_hits <- is.na(new_eval$pm_hits) | (new_eval$pm_hits <= max_hits)
  new_eval$is_pareto <- compute_pareto_frontier(new_eval, max_hits = max_hits)

  pareto_df <- new_eval |>
    dplyr::filter(.data$is_pareto) |>
    dplyr::arrange(dplyr::desc(.data$combined_recall), .data$pm_hits, .data$gs_chars) |>
    dplyr::mutate(rank = dplyr::row_number())

  rec <- recompute_recommendation(pareto_df, min_recall = min_recall, max_hits = max_hits)
  result$eval_df <- new_eval
  result$pareto_df <- rec$pareto_df
  result$recommended_sig <- rec$recommended_sig
  result$recommended_status <- rec$recommended_status
  result$last_added_eval_id <- custom_row$eval_id[1L]
  result
}

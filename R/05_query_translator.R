# R/05_query_translator.R
# Canonical query object and semantically equivalent translations to PubMed
# and Google Scholar (<= 256 chars) + URL builders.

#' Construct a canonical query object
#'
#' @param blocks List of blocks, each `list(name = ..., field = "tiab"|"title", terms = character())`.
#' @param exclude_terms Optional character vector of NOT terms.
#' @return Canonical query list of class `"search_query"`.
#' @examples
#' q <- make_query(
#'   blocks = list(
#'     list(name = "topic", field = "tiab", terms = c("h5n1", "h5")),
#'     list(name = "outcome", field = "tiab", terms = "fatal*")
#'   ),
#'   exclude_terms = c("h7n9", "vaccin*")
#' )
#' to_pubmed_query(q)
#' to_google_scholar_query(q)
make_query <- function(blocks, exclude_terms = character(0)) {
  clean_blocks <- lapply(seq_along(blocks), function(i) {
    b <- blocks[[i]]
    fld <- match.arg(b$field %||% "tiab", c("tiab", "title"))
    t_vec <- unique(tolower(stringr::str_squish(as.character(b$terms %||% character(0)))))
    t_vec <- t_vec[nzchar(t_vec)]
    # Drop terms within the same OR block that are strict superstrings or wildcard matches of another term in that block
    if (length(t_vec) > 1L) {
      keep_t <- vapply(seq_along(t_vec), function(idx) {
        t_cur <- t_vec[idx]
        t_cur_norm <- stringr::str_squish(stringr::str_replace_all(t_cur, "[-/]+", " "))
        others <- t_vec[-idx]
        !any(vapply(others, function(o) {
          if (identical(o, t_cur)) {
            return(FALSE)
          }
          o_norm <- stringr::str_squish(stringr::str_replace_all(o, "[-/]+", " "))
          has_star <- grepl("\\*", o_norm)
          if (!has_star && nchar(o_norm) >= nchar(t_cur_norm)) {
            return(FALSE)
          }
          pat <- escape_regex(o_norm)
          pat <- stringr::str_replace_all(pat, "\\\\\\*", "[a-z0-9]*")
          pat <- stringr::str_replace_all(pat, " ", "[\\\\s-]+")
          grepl(paste0("(?<![a-z0-9])", pat, "(?![a-z0-9])"), t_cur_norm, perl = TRUE)
        }, logical(1)))
      }, logical(1))
      t_vec <- t_vec[keep_t]
    }
    list(
      name = b$name %||% paste0("concept_", i),
      field = fld,
      terms = t_vec
    )
  })
  clean_blocks <- Filter(function(b) length(b$terms) > 0L, clean_blocks)
  excl <- unique(tolower(stringr::str_squish(as.character(exclude_terms %||% character(0)))))
  excl <- excl[nzchar(excl)]

  structure(
    list(blocks = clean_blocks, exclude_terms = excl),
    class = "search_query"
  )
}

#' Deterministic canonical signature for deduplicating queries
#'
#' @param q Canonical `search_query` object.
#' @return Character scalar signature.
query_signature <- function(q) {
  b_sigs <- vapply(q$blocks, function(b) {
    paste0(b$field, ":", paste(sort(b$terms), collapse = "|"))
  }, character(1))
  paste0(paste(sort(b_sigs), collapse = "&&"), "!!", paste(sort(q$exclude_terms), collapse = "|"))
}

# Parse optional "title:" prefix on an exclude term
.parse_exclude_term <- function(x) {
  x <- tolower(stringr::str_squish(as.character(x %||% "")))
  if (startsWith(x, "title:")) {
    list(field = "title", term = stringr::str_squish(sub("^title:", "", x)))
  } else if (startsWith(x, "ti:")) {
    list(field = "title", term = stringr::str_squish(sub("^ti:", "", x)))
  } else if (startsWith(x, "tiab:")) {
    list(field = "tiab", term = stringr::str_squish(sub("^tiab:", "", x)))
  } else {
    list(field = "tiab", term = x)
  }
}

# Normalize a prefix wildcard stem to a clean word for engines without prefix '*' (Google Scholar)
.expand_wildcard_stem <- function(term) {
  t_clean <- stringr::str_squish(stringr::str_replace_all(term, '["\',]', " "))
  if (!grepl("\\*$", t_clean)) {
    return(t_clean)
  }
  stem <- sub("\\*+$", "", t_clean)
  stem_map <- c(
    "mortalit" = "mortality",
    "morbidit" = "morbidity",
    "epidemiol" = "epidemiology",
    "epidemio" = "epidemiology",
    "seroprev" = "seroprevalence",
    "serolog" = "serological",
    "serol" = "serological",
    "pathogen" = "pathogenic",
    "transmi" = "transmission",
    "zoonot" = "zoonotic",
    "virulen" = "virulence",
    "hospital" = "hospitalization",
    "severit" = "severity",
    "sever" = "severe",
    "surviv" = "survival",
    "vaccin" = "vaccine",
    "immun" = "immune",
    "immuni" = "immunity",
    "antibod" = "antibody",
    "mutat" = "mutation",
    "phylo" = "phylogenetic",
    "phyloge" = "phylogenetic",
    "reassort" = "reassortment",
    "epithel" = "epithelial",
    "epitheli" = "epithelial",
    "structur" = "structural",
    "express" = "expression",
    "nucle" = "nucleoprotein",
    "necro" = "necrosis",
    "interfer" = "interferon",
    "cytokin" = "cytokine",
    "bovin" = "bovine"
  )
  if (stem %in% names(stem_map)) unname(stem_map[[stem]]) else stem
}

# Format a term for Google Scholar (quote only if multi-word or contains non-alphanumeric chars)
.gs_format_term <- function(term, in_title = FALSE) {
  t_clean <- .expand_wildcard_stem(term)
  needs_quotes <- grepl("[^a-zA-Z0-9]", t_clean)
  core <- if (needs_quotes) paste0('"', t_clean, '"') else t_clean
  if (isTRUE(in_title)) paste0("intitle:", core) else core
}

#' Translate canonical query to PubMed E-utilities query string
#'
#' @param q Canonical `search_query` object.
#' @return Character scalar.
to_pubmed_query <- function(q) {
  if (length(q$blocks) == 0L) {
    return("")
  }
  multi <- length(q$blocks) > 1L || length(q$exclude_terms) > 0L
  b_strs <- vapply(q$blocks, function(b) {
    tag <- if (b$field == "title") "[ti]" else "[tiab]"
    parts <- sprintf('"%s"%s', b$terms, tag)
    joined <- paste(parts, collapse = " OR ")
    if (length(parts) > 1L || multi) paste0("(", joined, ")") else joined
  }, character(1))

  out <- paste(b_strs, collapse = " AND ")
  if (length(q$exclude_terms) > 0L) {
    not_parts <- vapply(q$exclude_terms, function(ex) {
      p <- .parse_exclude_term(ex)
      tag <- if (p$field == "title") "[ti]" else "[tiab]"
      sprintf('"%s"%s', p$term, tag)
    }, character(1))
    not_expr <- paste(not_parts, collapse = " OR ")
    out <- paste0("(", out, ") NOT (", not_expr, ")")
  }
  out
}

#' Translate canonical query to Google Scholar query string
#'
#' Uses implicit AND (space) between concept blocks, uppercase OR, `intitle:` for
#' title-restricted blocks, and unquoted clean single tokens to stay within 256 chars.
#'
#' @param q Canonical `search_query` object.
#' @return Character scalar.
to_google_scholar_query <- function(q) {
  if (length(q$blocks) == 0L) {
    return("")
  }
  b_strs <- vapply(q$blocks, function(b) {
    in_title <- identical(b$field, "title")
    parts <- vapply(b$terms, .gs_format_term, character(1), in_title = in_title)
    joined <- paste(parts, collapse = " OR ")
    if (length(parts) > 1L) paste0("(", joined, ")") else joined
  }, character(1))

  out <- paste(b_strs, collapse = " ")
  if (length(q$exclude_terms) > 0L) {
    excl_parts <- vapply(q$exclude_terms, function(t_str) {
      p <- .parse_exclude_term(t_str)
      paste0("-", .gs_format_term(p$term, in_title = identical(p$field, "title")))
    }, character(1))
    out <- paste(c(out, excl_parts), collapse = " ")
  }
  stringr::str_squish(out)
}

#' Count characters in translated Google Scholar query
#'
#' @param q Canonical `search_query` object.
#' @return Integer character count.
gs_char_count <- function(q) {
  nchar(to_google_scholar_query(q), type = "chars")
}

#' Check Google Scholar <= 256 character constraint
#'
#' @param q Canonical `search_query` object.
#' @param max_gs_chars Maximum allowed character count (default 256).
#' @return Logical scalar.
is_valid_gs_query <- function(q, max_gs_chars = 256L) {
  gs <- to_google_scholar_query(q)
  nzchar(gs) && nchar(gs, type = "chars") <= as.integer(max_gs_chars)
}

#' Clickable PubMed search URL
#'
#' @param pm_query PubMed query string.
#' @return URL string.
pubmed_url <- function(pm_query) {
  paste0("https://pubmed.ncbi.nlm.nih.gov/?term=", utils::URLencode(pm_query, reserved = TRUE))
}

#' Clickable Google Scholar search URL
#'
#' @param gs_query Google Scholar query string.
#' @return URL string.
google_scholar_url <- function(gs_query) {
  paste0("https://scholar.google.com/scholar?q=", utils::URLencode(gs_query, reserved = TRUE))
}

#' Convert a raw PubMed query string into an approximate Google Scholar query
#'
#' Preserves boolean OR/AND/NOT semantics, converts `[ti]`/`[Title]` tags to
#' `intitle:`, strips `[tiab]`/`[MeSH]` tags and `*` wildcards, converts `AND` to
#' implicit spaces, and formats `NOT (...)` exclusions as `-term`.
#'
#' @param pm_query Raw PubMed query string.
#' @return Approximate Google Scholar query string.
#' @examples
#' pubmed_to_google_scholar_approx(
#'   '(("H5N1"[tiab] OR "H5"[tiab]) AND "fatal*"[tiab]) NOT "H7N9"[tiab]'
#' )
pubmed_to_google_scholar_approx <- function(pm_query) {
  q <- stringr::str_squish(as.character(pm_query %||% ""))
  if (!nzchar(q)) {
    return("")
  }

  # Convert title tags before stripping general tags (vectorized over all matches in `m`)
  q <- stringr::str_replace_all(q, '(?i)"([^"]+)"\\s*\\[(?:ti|title)\\]', function(m) {
    inner <- stringr::str_match(m, '(?i)"([^"]+)"')[, 2]
    inner <- sub("\\*$", "", inner)
    ifelse(grepl("\\s", inner), paste0('intitle:"', inner, '"'), paste0("intitle:", inner))
  })
  q <- stringr::str_replace_all(q, "(?i)([A-Za-z0-9_*-]+)\\s*\\[(?:ti|title)\\]", function(m) {
    inner <- stringr::str_match(m, "(?i)([A-Za-z0-9_*-]+)\\s*\\[")[, 2]
    paste0("intitle:", sub("\\*$", "", inner))
  })

  # Strip remaining PubMed field tags like [Title/Abstract], [tiab], [MeSH], [tw], etc.
  q <- stringr::str_replace_all(q, "\\[[^\\]]+\\]", "")
  # Strip truncation wildcards
  q <- stringr::str_replace_all(q, "([A-Za-z0-9_-]+)\\*", "\\1")
  # Unquote single clean alphanumeric tokens (save GS characters)
  q <- stringr::str_replace_all(q, '"([A-Za-z0-9_-]+)"', "\\1")
  q <- stringr::str_squish(q)

  # Split on top-level NOT if present
  parts <- stringr::str_split(q, "(?i)\\s+NOT\\s+", n = 2L)[[1]]
  pos_part <- stringr::str_squish(parts[1L])
  neg_part <- if (length(parts) >= 2L) stringr::str_squish(parts[2L]) else ""

  # Convert AND -> space, normalize OR spacing
  pos_part <- stringr::str_replace_all(pos_part, "(?i)\\s+AND\\s+", " ")
  pos_part <- stringr::str_replace_all(pos_part, "(?i)\\s+OR\\s+", " OR ")
  pos_part <- stringr::str_replace_all(pos_part, "\\(\\s+", "(")
  pos_part <- stringr::str_replace_all(pos_part, "\\s+\\)", ")")
  # Unwrap single redundant outer parens around positive part when NOT was stripped
  if (nzchar(neg_part) && startsWith(pos_part, "(") && endsWith(pos_part, ")")) {
    inner <- substr(pos_part, 2L, nchar(pos_part) - 1L)
    chars <- strsplit(inner, "", fixed = TRUE)[[1]]
    bal <- cumsum(ifelse(chars == "(", 1L, ifelse(chars == ")", -1L, 0L)))
    if (length(bal) > 0L && all(bal >= 0L)) pos_part <- inner
  }

  excl_str <- ""
  if (nzchar(neg_part)) {
    neg_clean <- stringr::str_replace_all(neg_part, "[()]", " ")
    neg_clean <- stringr::str_replace_all(neg_clean, "(?i)\\b(?:OR|AND|NOT)\\b", " ")
    toks <- stringr::str_extract_all(neg_clean, '(?:intitle:)?(?:"[^"]+"|[^\\s]+)')[[1]]
    toks <- toks[nzchar(toks)]
    if (length(toks) > 0L) {
      excl_str <- paste(paste0("-", toks), collapse = " ")
    }
  }

  stringr::str_squish(paste(pos_part, excl_str))
}

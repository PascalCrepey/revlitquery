# R/01_bibtex_parser.R
# Character-by-character BibTeX parser with LaTeX cleaning, DOI normalization,
# and PMID extraction.

#' Clean LaTeX markup and normalize whitespace in BibTeX field text
#'
#' @param x Character scalar or vector.
#' @return Cleaned character vector.
clean_latex_text <- function(x) {
  if (is.null(x) || length(x) == 0L) return(character(0))
  out <- as.character(x)
  out[is.na(out)] <- ""

  # Common LaTeX accent sequences (with or without braces)
  accents <- list(
    c("\\\\'\\{?a\\}?", "a"), c("\\\\'\\{?e\\}?", "e"), c("\\\\'\\{?i\\}?", "i"),
    c("\\\\'\\{?o\\}?", "o"), c("\\\\'\\{?u\\}?", "u"), c("\\\\'\\{?y\\}?", "y"),
    c("\\\\`\\{?a\\}?", "a"), c("\\\\`\\{?e\\}?", "e"), c("\\\\`\\{?i\\}?", "i"),
    c("\\\\`\\{?o\\}?", "o"), c("\\\\`\\{?u\\}?", "u"),
    c('\\\\"\\{?a\\}?', "a"), c('\\\\"\\{?e\\}?', "e"), c('\\\\"\\{?i\\}?', "i"),
    c('\\\\"\\{?o\\}?', "o"), c('\\\\"\\{?u\\}?', "u"), c('\\\\"\\{?y\\}?', "y"),
    c("\\\\\\^\\{?a\\}?", "a"), c("\\\\\\^\\{?e\\}?", "e"), c("\\\\\\^\\{?i\\}?", "i"),
    c("\\\\\\^\\{?o\\}?", "o"), c("\\\\\\^\\{?u\\}?", "u"),
    c("\\\\~\\{?a\\}?", "a"), c("\\\\~\\{?n\\}?", "n"), c("\\\\~\\{?o\\}?", "o"),
    c("\\\\c\\{?c\\}?", "c"), c("\\\\c\\{?C\\}?", "C"),
    c("\\\\'\\{?A\\}?", "A"), c("\\\\'\\{?E\\}?", "E"), c("\\\\'\\{?I\\}?", "I"),
    c("\\\\'\\{?O\\}?", "O"), c("\\\\'\\{?U\\}?", "U"),
    c('\\\\"\\{?A\\}?', "A"), c('\\\\"\\{?O\\}?', "O"), c('\\\\"\\{?U\\}?', "U"),
    c("\\\\ss\\b", "ss"), c("\\\\ae\\b", "ae"), c("\\\\oe\\b", "oe"),
    c("\\\\o\\b", "o"), c("\\\\O\\b", "O"), c("\\\\l\\b", "l"), c("\\\\L\\b", "L")
  )
  for (pair in accents) {
    out <- stringr::str_replace_all(out, pair[[1]], pair[[2]])
  }

  # Common LaTeX escapes and dashes
  out <- stringr::str_replace_all(out, "\\\\&", "&")
  out <- stringr::str_replace_all(out, "\\\\%", "%")
  out <- stringr::str_replace_all(out, "\\\\_", "_")
  out <- stringr::str_replace_all(out, "\\\\#", "#")
  out <- stringr::str_replace_all(out, "\\\\textendash\\b|---?", "-")
  out <- stringr::str_replace_all(out, "\\\\textemdash\\b", "-")
  out <- stringr::str_replace_all(out, "\\\\textit\\{([^}]*)\\}", "\\1")
  out <- stringr::str_replace_all(out, "\\\\textbf\\{([^}]*)\\}", "\\1")
  out <- stringr::str_replace_all(out, "\\\\emph\\{([^}]*)\\}", "\\1")
  out <- stringr::str_replace_all(out, "\\\\textsuperscript\\{([^}]*)\\}", "\\1")
  out <- stringr::str_replace_all(out, "\\\\textsubscript\\{([^}]*)\\}", "\\1")
  out <- stringr::str_replace_all(out, "~", " ")

  # Strip remaining curly braces and stray inline math delimiters
  out <- stringr::str_replace_all(out, "[{}]", "")
  out <- stringr::str_replace_all(out, "\\s+", " ")
  stringr::str_trim(out)
}

#' Normalize a DOI string to lowercase bare identifier
#'
#' @param x Character scalar or vector.
#' @return Normalized DOI string (empty string if missing).
normalize_doi <- function(x) {
  if (is.null(x) || length(x) == 0L) return(character(0))
  out <- as.character(x)
  out[is.na(out)] <- ""
  out <- stringr::str_trim(out)
  out <- stringr::str_replace(out, "^https?://(dx\\.)?doi\\.org/", "")
  out <- stringr::str_replace(out, "^doi:\\s*", "")
  out <- stringr::str_replace_all(out, "[{}]", "")
  tolower(stringr::str_trim(out))
}

#' Extract a PubMed ID (PMID) from explicit fields or URL/note text
#'
#' @param pmid_field Explicit pmid/pubmedid field.
#' @param url_field Optional url field.
#' @param note_field Optional note field.
#' @return Character scalar (digits or "").
extract_pmid_value <- function(pmid_field = "", url_field = "", note_field = "") {
  cand <- stringr::str_trim(as.character(pmid_field %||% ""))
  if (nzchar(cand)) {
    m <- stringr::str_match(cand, "(\\d{5,10})")
    if (!is.na(m[1, 2])) return(m[1, 2])
  }
  for (txt in c(url_field %||% "", note_field %||% "")) {
    txt <- as.character(txt)
    m <- stringr::str_match(txt, "pubmed\\.ncbi\\.nlm\\.nih\\.gov/(\\d{5,10})")
    if (!is.na(m[1, 2])) return(m[1, 2])
    m2 <- stringr::str_match(txt, "(?i)\\bPMID[:\\s]+(\\d{5,10})\\b")
    if (!is.na(m2[1, 2])) return(m2[1, 2])
  }
  ""
}

#' Split BibTeX keywords string into cleaned lowercase character vector
#'
#' @param kw_str Raw keywords string.
#' @return Character vector of unique keywords.
split_keywords <- function(kw_str) {
  if (is.null(kw_str) || is.na(kw_str) || !nzchar(stringr::str_trim(kw_str))) {
    return(character(0))
  }
  parts <- unlist(stringr::str_split(clean_latex_text(kw_str), "\\s*[;,]\\s*"))
  parts <- tolower(stringr::str_trim(parts))
  unique(parts[nzchar(parts)])
}

# Null-coalescing helper
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L) b else a
}

# Internal: skip whitespace and optional BibTeX % comments outside values
.skip_ws <- function(chars, pos, n) {
  while (pos <= n && grepl("\\s", chars[pos])) {
    pos <- pos + 1L
  }
  pos
}

# Internal: read a BibTeX field value starting at `pos` (handles {...}, "...", bare, and # concat)
.read_bib_value <- function(chars, pos, n) {
  parts <- character(0)
  repeat {
    pos <- .skip_ws(chars, pos, n)
    if (pos > n) break
    ch <- chars[pos]
    if (ch == "{") {
      depth <- 1L
      pos <- pos + 1L
      start <- pos
      while (pos <= n && depth > 0L) {
        c_cur <- chars[pos]
        if (c_cur == "\\" && pos < n) {
          pos <- pos + 2L
          next
        }
        if (c_cur == "{") {
          depth <- depth + 1L
        } else if (c_cur == "}") {
          depth <- depth - 1L
          if (depth == 0L) break
        }
        pos <- pos + 1L
      }
      val_part <- if (pos > start) paste0(chars[start:(pos - 1L)], collapse = "") else ""
      parts <- c(parts, val_part)
      if (pos <= n && chars[pos] == "}") pos <- pos + 1L
    } else if (ch == '"') {
      pos <- pos + 1L
      start <- pos
      depth <- 0L
      while (pos <= n) {
        c_cur <- chars[pos]
        if (c_cur == "\\" && pos < n) {
          pos <- pos + 2L
          next
        }
        if (c_cur == "{") depth <- depth + 1L
        else if (c_cur == "}" && depth > 0L) depth <- depth - 1L
        else if (c_cur == '"' && depth == 0L) break
        pos <- pos + 1L
      }
      val_part <- if (pos > start) paste0(chars[start:(pos - 1L)], collapse = "") else ""
      parts <- c(parts, val_part)
      if (pos <= n && chars[pos] == '"') pos <- pos + 1L
    } else {
      start <- pos
      while (pos <= n && !(chars[pos] %in% c(",", "}", ")", "#"))) {
        pos <- pos + 1L
      }
      val_part <- if (pos > start) stringr::str_trim(paste0(chars[start:(pos - 1L)], collapse = "")) else ""
      parts <- c(parts, val_part)
    }

    pos <- .skip_ws(chars, pos, n)
    if (pos <= n && chars[pos] == "#") {
      pos <- pos + 1L
      next
    }
    break
  }
  list(value = paste0(parts, collapse = ""), pos = pos)
}

# Internal: parse body of a single BibTeX entry (citekey, field = value, ...)
.parse_entry_body <- function(body_str) {
  chars <- strsplit(body_str, "", fixed = TRUE)[[1]]
  n <- length(chars)
  pos <- .skip_ws(chars, 1L, n)

  # Read citekey until comma
  key_start <- pos
  while (pos <= n && chars[pos] != ",") {
    pos <- pos + 1L
  }
  citekey <- if (pos > key_start) stringr::str_trim(paste0(chars[key_start:(pos - 1L)], collapse = "")) else ""
  if (pos <= n && chars[pos] == ",") pos <- pos + 1L

  fields <- list()
  while (pos <= n) {
    pos <- .skip_ws(chars, pos, n)
    if (pos > n) break
    if (chars[pos] == ",") {
      pos <- pos + 1L
      next
    }
    # Read field name
    name_start <- pos
    while (pos <= n && grepl("[A-Za-z0-9_:-]", chars[pos])) {
      pos <- pos + 1L
    }
    if (pos == name_start) {
      pos <- pos + 1L
      next
    }
    field_name <- tolower(paste0(chars[name_start:(pos - 1L)], collapse = ""))
    pos <- .skip_ws(chars, pos, n)
    if (pos > n || chars[pos] != "=") {
      # Skip until next comma if malformed
      while (pos <= n && chars[pos] != ",") pos <- pos + 1L
      next
    }
    pos <- pos + 1L # consume '='
    res <- .read_bib_value(chars, pos, n)
    fields[[field_name]] <- res$value
    pos <- res$pos
    pos <- .skip_ws(chars, pos, n)
    if (pos <= n && chars[pos] == ",") pos <- pos + 1L
  }
  list(citekey = citekey, fields = fields)
}

#' Parse a BibTeX file into a normalized tibble of seed papers
#'
#' @param path Path to a `.bib` file.
#' @return A tibble with columns `seed_id`, `citekey`, `entry_type`, `title`,
#'   `abstract`, `doi`, `pmid`, `keywords`, `year`, `author`, `journal`.
parse_bibtex <- function(path) {
  if (!file.exists(path)) {
    stop("BibTeX file not found: ", path, call. = FALSE)
  }
  raw <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Strip UTF-8 BOM if present
  raw <- sub("^\ufeff", "", raw)
  chars <- strsplit(raw, "", fixed = TRUE)[[1]]
  n <- length(chars)
  pos <- 1L

  entries <- list()
  while (pos <= n) {
    # Scan for '@'
    while (pos <= n && chars[pos] != "@") {
      pos <- pos + 1L
    }
    if (pos > n) break
    pos <- pos + 1L # consume '@'
    pos <- .skip_ws(chars, pos, n)

    type_start <- pos
    while (pos <= n && grepl("[A-Za-z0-9_]", chars[pos])) {
      pos <- pos + 1L
    }
    if (pos == type_start) next
    entry_type <- tolower(paste0(chars[type_start:(pos - 1L)], collapse = ""))
    pos <- .skip_ws(chars, pos, n)
    if (pos > n) break

    open_ch <- chars[pos]
    if (!(open_ch %in% c("{", "("))) next
    close_ch <- if (open_ch == "{") "}" else ")"
    pos <- pos + 1L
    body_start <- pos
    depth <- 1L
    in_quotes <- FALSE

    while (pos <= n && depth > 0L) {
      c_cur <- chars[pos]
      if (c_cur == "\\" && pos < n) {
        pos <- pos + 2L
        next
      }
      if (c_cur == '"' && depth == 1L) {
        in_quotes <- !in_quotes
      } else if (!in_quotes) {
        if (c_cur == open_ch) {
          depth <- depth + 1L
        } else if (c_cur == close_ch) {
          depth <- depth - 1L
          if (depth == 0L) break
        }
      }
      pos <- pos + 1L
    }

    body_str <- if (pos > body_start) paste0(chars[body_start:(pos - 1L)], collapse = "") else ""
    if (pos <= n && chars[pos] == close_ch) pos <- pos + 1L

    if (entry_type %in% c("comment", "preamble", "string")) next

    parsed <- .parse_entry_body(body_str)
    f <- parsed$fields

    title <- clean_latex_text(f$title %||% "")
    abstract <- clean_latex_text(f$abstract %||% f$summary %||% "")
    doi <- normalize_doi(f$doi %||% "")
    pmid <- extract_pmid_value(
      pmid_field = f$pmid %||% f$pubmedid %||% "",
      url_field = f$url %||% "",
      note_field = f$note %||% ""
    )
    kw_vec <- split_keywords(f$keywords %||% f$keyword %||% "")
    year <- clean_latex_text(f$year %||% "")
    author <- clean_latex_text(f$author %||% "")
    journal <- clean_latex_text(f$journal %||% f$booktitle %||% f$journaltitle %||% "")

    entries[[length(entries) + 1L]] <- list(
      citekey = if (nzchar(parsed$citekey)) parsed$citekey else paste0("seed_", length(entries) + 1L),
      entry_type = entry_type,
      title = title,
      abstract = abstract,
      doi = doi,
      pmid = pmid,
      keywords = paste(kw_vec, collapse = "; "),
      year = year,
      author = author,
      journal = journal
    )
  }

  if (length(entries) == 0L) {
    stop("No valid BibTeX entries found in: ", path, call. = FALSE)
  }

  tibble::tibble(
    seed_id = seq_along(entries),
    citekey = vapply(entries, `[[`, character(1), "citekey"),
    entry_type = vapply(entries, `[[`, character(1), "entry_type"),
    title = vapply(entries, `[[`, character(1), "title"),
    abstract = vapply(entries, `[[`, character(1), "abstract"),
    doi = vapply(entries, `[[`, character(1), "doi"),
    pmid = vapply(entries, `[[`, character(1), "pmid"),
    keywords = vapply(entries, `[[`, character(1), "keywords"),
    year = vapply(entries, `[[`, character(1), "year"),
    author = vapply(entries, `[[`, character(1), "author"),
    journal = vapply(entries, `[[`, character(1), "journal")
  )
}

#!/usr/bin/env Rscript
# run_query_optimizer.R
# Command-line wrapper for revlitquery::optimize_search_queries().

.get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE)))
  }
  "."
}

.load_optimizer_modules <- function(base_dir = .get_script_dir()) {
  r_dir <- file.path(base_dir, "R")
  if (!dir.exists(r_dir) && dir.exists("R")) r_dir <- "R"
  files <- sort(list.files(r_dir, pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
  if (length(files) == 0L) {
    stop("Could not find R/*.R modules in: ", r_dir, call. = FALSE)
  }
  for (f in files) source(f, local = FALSE)
}

.load_optimizer_modules()

.print_usage <- function() {
  cat(
    "Usage: Rscript run_query_optimizer.R --bib <file.bib> [options]\n\n",
    "Options:\n",
    "  --bib <path>           Path to BibTeX input file (required)\n",
    "  --out <dir>            Output directory (default: output)\n",
    "  --config <path>        Optional YAML configuration file\n",
    "  --min-recall <num>     Minimum seed recall 0..1 (default: 0.60)\n",
    "  --max-hits <num>       Maximum PubMed hits cap (default: 2000; Inf for no cap)\n",
    "  --max-evals <int>      Maximum live query evaluations (default: 60)\n",
    "  --max-gs-chars <int>   Max Google Scholar query characters (default: 256)\n",
    "  --custom-query <str>   Custom PubMed query to evaluate (repeatable)\n",
    "  --custom-label <str>   Optional label for next/last --custom-query\n",
    "  --no-cache             Disable disk API response caching\n",
    "  --email <str>          Email for NCBI polite pool\n",
    "  --ncbi-key <str>       NCBI API key (or set NCBI_API_KEY env var)\n",
    "  --help, -h             Show this help message and exit\n",
    sep = ""
  )
}

.parse_cli_args <- function(args) {
  opts <- list(
    bib = NULL,
    out = "output",
    config = NULL,
    min_recall = 0.60,
    max_hits = 2000,
    max_evals = 60L,
    max_gs_chars = 256L,
    use_cache = TRUE,
    custom_queries = list(),
    pending_label = NULL,
    email = Sys.getenv("NCBI_EMAIL", ""),
    ncbi_key = Sys.getenv("NCBI_API_KEY", ""),
    help = FALSE
  )

  i <- 1L
  n <- length(args)
  while (i <= n) {
    arg <- args[i]
    if (arg %in% c("--help", "-h")) {
      opts$help <- TRUE
      i <- i + 1L
    } else if (arg == "--no-cache") {
      opts$use_cache <- FALSE
      i <- i + 1L
    } else if (startsWith(arg, "--")) {
      if (grepl("=", arg, fixed = TRUE)) {
        parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
        key <- parts[1L]
        val <- paste(parts[-1L], collapse = "=")
        i <- i + 1L
      } else {
        key <- sub("^--", "", arg)
        if (i == n) stop("Missing value for flag: ", arg, call. = FALSE)
        val <- args[i + 1L]
        i <- i + 2L
      }
      if (key == "bib") opts$bib <- val
      else if (key == "out") opts$out <- val
      else if (key == "config") opts$config <- val
      else if (key == "min-recall") opts$min_recall <- as.numeric(val)
      else if (key == "max-hits") opts$max_hits <- as.numeric(val)
      else if (key == "max-evals") opts$max_evals <- as.integer(val)
      else if (key == "max-gs-chars") opts$max_gs_chars <- as.integer(val)
      else if (key == "custom-label") {
        if (length(opts$custom_queries) > 0L && is.null(opts$pending_label)) {
          last_i <- length(opts$custom_queries)
          if (identical(opts$custom_queries[[last_i]]$label, paste0("Custom ", last_i))) {
            opts$custom_queries[[last_i]]$label <- val
          } else {
            opts$pending_label <- val
          }
        } else {
          opts$pending_label <- val
        }
      } else if (key == "custom-query") {
        idx_q <- length(opts$custom_queries) + 1L
        lbl <- opts$pending_label %||% paste0("Custom ", idx_q)
        opts$pending_label <- NULL
        opts$custom_queries[[idx_q]] <- list(label = lbl, pubmed_query = val)
      }
      else if (key == "email") opts$email <- val
      else if (key == "ncbi-key") opts$ncbi_key <- val
      else stop("Unknown CLI option: --", key, call. = FALSE)
    } else {
      stop("Unexpected positional argument: ", arg, call. = FALSE)
    }
  }
  opts
}

if (sys.nframe() == 0L && !interactive()) {
  cli_args <- commandArgs(trailingOnly = TRUE)
  parsed <- .parse_cli_args(cli_args)
  if (isTRUE(parsed$help) || is.null(parsed$bib)) {
    .print_usage()
    if (!isTRUE(parsed$help) && is.null(parsed$bib)) {
      quit(save = "no", status = 1L)
    }
    quit(save = "no", status = 0L)
  }

  optimize_search_queries(
    bib_path = parsed$bib,
    output_dir = parsed$out,
    config_path = parsed$config,
    min_recall = parsed$min_recall,
    max_hits = parsed$max_hits,
    max_evals = parsed$max_evals,
    max_gs_chars = parsed$max_gs_chars,
    use_cache = parsed$use_cache,
    custom_queries = parsed$custom_queries,
    email = parsed$email,
    ncbi_key = parsed$ncbi_key
  )
}

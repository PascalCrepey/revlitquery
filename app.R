# app.R
# Interactive Shiny launcher for revlitquery (can also be launched via revlitquery::run_revlitquery_app())

if (!requireNamespace("revlitquery", quietly = TRUE) || !exists("revlitquery_app", mode = "function")) {
  r_dir <- if (dir.exists("R")) "R" else file.path("..", "R")
  files <- sort(list.files(r_dir, pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
  if (length(files) == 0L) stop("Could not locate R/*.R modules from ", getwd(), call. = FALSE)
  for (f in files) source(f, local = FALSE)
}

revlitquery_app()

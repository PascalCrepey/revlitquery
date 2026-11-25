library(testthat)

if (file.exists("R/01_bibtex_parser.R") || file.exists("../R/01_bibtex_parser.R")) {
  proj_root <- if (file.exists("R/01_bibtex_parser.R")) "." else ".."
  r_files <- sort(list.files(file.path(proj_root, "R"), pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
  for (f in r_files) source(f, local = FALSE)
  test_dir(file.path(proj_root, "tests", "testthat"), stop_on_failure = TRUE)
} else {
  library(revlitquery)
  test_check("revlitquery")
}

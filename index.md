# revlitquery

Given a `.bib` file of target/seed papers from a systematic or
literature review, **`revlitquery`** reverse-engineers semantically
equivalent **PubMed** and **Google Scholar** (`<= 256` characters)
boolean search queries and optimizes their **sensitivity (seed recall)**
vs. **specificity (total hit count)** Pareto frontier using NCBI
E-utilities:

- **PubMed (NCBI E-utilities)**: live total hit counts, single-request
  seed PMID retrieval checks, year-assisted title resolution, and
  per-PMID disk caching.
- **Google Scholar**: exact semantic translation adhering to Google
  Scholar syntax (`intitle:`, implicit `AND`, uppercase `OR`,
  `-exclude`, `<= 256` characters) with corpus-derived wildcard stem
  completion and clickable verification links.
- **Unsupervised Inclusion & Exclusion Concept Discovery**:
  domain-agnostic distributional context (`SVD`/`LSA` + MeSH/Keyword
  projection + seed complementarity) Louvain graph clustering with
  orthogonal set-cover for inclusion `AND`/`OR` blocks, paired with
  Normalized Pointwise Mutual Information (`NPMI`) co-occurrence topic
  clustering on non-seed PubMed results for zero-loss `NOT` / `-`
  exclusion blocks across **any scientific discipline**.

------------------------------------------------------------------------

## Installation

Install `revlitquery` from GitHub using `remotes` or `pak` (pass
`build_vignettes = TRUE` to build the Shiny interface and CLI/R workflow
vignettes locally):

``` r

# install.packages("remotes")
remotes::install_github("PascalCrepey/revlitquery", build_vignettes = TRUE)
```

Or install locally from a cloned copy of the repository:

``` bash
R CMD build .
R CMD INSTALL revlitquery_0.1.0.tar.gz
```

## Website & Vignettes

Browse the full documentation website at
**<https://pascalcrepey.github.io/revlitquery/>**, or view the included
vignettes directly from R:

``` r

# 1. Interactive Shiny Interface Guide (with step-by-step screenshots)
vignette("shiny-interface", package = "revlitquery")

# 2. Command-Line (CLI) and Programmatic R Workflow
vignette("cli-and-r-workflow", package = "revlitquery")
```

------------------------------------------------------------------------

## Quick Start in R

### 1. Launch the Interactive Shiny Dashboard

``` r

library(revlitquery)

run_revlitquery_app()
```

Or directly from a local clone without installing:

``` r

shiny::runApp("app.R", launch.browser = TRUE)
```

### 2. Run the Programmatic Optimizer

``` r

library(revlitquery)

res <- optimize_search_queries(
  bib_path   = "path/to/seed_articles.bib",
  output_dir = "output",
  min_recall = 0.90,
  max_hits   = 1000,
  max_evals  = 60
)
```

### 3. Command-Line Interface (CLI)

``` bash
Rscript run_query_optimizer.R --bib path/to/seed_articles.bib --min-recall 0.90 --max-hits 1000 --out output
```

------------------------------------------------------------------------

## Controlling Specificity vs. Sensitivity

- By default, the optimizer enforces a ceiling of `max_hits = 2000`
  PubMed hits on the Pareto frontier and automatically mines non-seed
  PubMed results for concept-organized safe exclusion (`NOT` / `-`)
  terms.
- Inspect `output/exclusion_hints.csv` or the **Exclusion Hints** tab in
  the Shiny app to see up to 35 ranked candidate exclusions (grouped by
  concept, supporting `[tiab]`, `title:`-scoped `-intitle:`, and prefix
  wildcards `*`).
- Benchmark your own PubMed query alongside the Pareto frontier via
  `custom_queries = "..."` in
  [`optimize_search_queries()`](https://pascalcrepey.github.io/revlitquery/reference/optimize_search_queries.md),
  `--custom-query "..."` on the CLI, or in the **Pareto Explorer** tab
  of the Shiny app.

## Optional Environment Variables

- `NCBI_EMAIL`: email sent to the NCBI E-utilities polite pool.
- `NCBI_API_KEY`: raises NCBI E-utilities rate limit from ~3 req/s to
  ~10 req/s.

## Outputs (`output/`)

| File | Description |
|----|----|
| `pareto_queries.csv` | Non-dominated PubMed + Google Scholar queries and URLs |
| `all_evaluated_queries.csv` | All queries evaluated during Stage 1 & Stage 2 search |
| `seed_resolution.csv` | Resolution status of each BibTeX seed (PMID, abstract) |
| `seed_coverage_matrix.csv` | Per-paper `TRUE/FALSE` coverage across each Pareto query |
| `exclusion_hints.csv` | Concept-grouped exclusion terms and prefix wildcards with hit reductions |
| `report.md` | Human-readable summary report with copy-pasteable queries & links |
| `pareto_frontier.png` | Log-scale total hits vs. seed recall plot |

## Running Tests

``` bash
Rscript tests/testthat.R
```

## License

MIT © Pascal Crépey

# Command-Line (CLI) and Programmatic R Workflow

## Overview

`revlitquery` reverse-engineers boolean literature search queries for
**PubMed** and **Google Scholar** (`<= 256` characters) starting from a
seed corpus of known relevant articles in BibTeX (`.bib`) format.

In addition to the interactive Shiny interface
([`run_revlitquery_app()`](https://pascalcrepey.github.io/revlitquery/reference/run_revlitquery_app.md)),
`revlitquery` can be driven completely from:

1.  **The Command-Line Interface (CLI)** via `run_query_optimizer.R` or
    `Rscript -e "revlitquery::optimize_search_queries(...)"`.
2.  **Programmatic R scripts** using
    [`optimize_search_queries()`](https://pascalcrepey.github.io/revlitquery/reference/optimize_search_queries.md)
    or the modular pipeline functions
    ([`parse_bibtex()`](https://pascalcrepey.github.io/revlitquery/reference/parse_bibtex.md),
    [`enrich_seeds()`](https://pascalcrepey.github.io/revlitquery/reference/enrich_seeds.md),
    [`extract_candidate_terms()`](https://pascalcrepey.github.io/revlitquery/reference/extract_candidate_terms.md),
    [`build_concept_blocks()`](https://pascalcrepey.github.io/revlitquery/reference/build_concept_blocks.md),
    [`optimize_queries()`](https://pascalcrepey.github.io/revlitquery/reference/optimize_queries.md),
    [`discover_exclusion_hints()`](https://pascalcrepey.github.io/revlitquery/reference/discover_exclusion_hints.md),
    and
    [`write_optimization_report()`](https://pascalcrepey.github.io/revlitquery/reference/write_optimization_report.md)).

------------------------------------------------------------------------

## 1. Command-Line Interface (CLI)

### Basic CLI Usage

From a terminal in the repository root, run `run_query_optimizer.R` with
`--bib` pointing to your BibTeX file of seed references:

``` bash
Rscript run_query_optimizer.R \
  --bib path/to/seed_articles.bib \
  --out output
```

Or directly via the installed `revlitquery` package from anywhere on
your system:

``` bash
Rscript -e 'revlitquery::optimize_search_queries(
  bib_path   = "path/to/seed_articles.bib",
  output_dir = "output",
  min_recall = 0.90,
  max_hits   = 1000
)'
```

### Complete CLI Options Reference

Run `Rscript run_query_optimizer.R --help` to display all supported
command-line flags:

``` text
Usage: Rscript run_query_optimizer.R --bib <file.bib> [options]

Options:
  --bib <path>           Path to BibTeX input file (required)
  --out <dir>            Output directory (default: output)
  --config <path>        Optional YAML configuration file
  --min-recall <num>     Minimum seed recall 0..1 (default: 0.60)
  --max-hits <num>       Maximum PubMed hits cap (default: 2000; Inf for no cap)
  --max-evals <int>      Maximum live query evaluations (default: 60)
  --max-gs-chars <int>   Max Google Scholar query characters (default: 256)
  --custom-query <str>   Custom PubMed query to evaluate (repeatable)
  --custom-label <str>   Optional label for next/last --custom-query
  --no-cache             Disable disk API response caching
  --email <str>          Email for NCBI polite pool
  --ncbi-key <str>       NCBI API key (or set NCBI_API_KEY env var)
  --help, -h             Show this help message and exit
```

### Common CLI Recipes

#### Enforcing High Recall with a Strict Specificity Cap (`<= 1000` Hits)

To require at least 95% seed recall while keeping total PubMed hits
below 1,000:

``` bash
Rscript run_query_optimizer.R \
  --bib path/to/seed_articles.bib \
  --min-recall 0.95 \
  --max-hits 1000 \
  --max-evals 60 \
  --out output_strict
```

#### Benchmarking Your Own Manual Query Against the Pareto Frontier

Use `--custom-label` and `--custom-query` (both repeatable) to evaluate
one or more manual queries alongside the automatically discovered Pareto
frontier:

``` bash
Rscript run_query_optimizer.R \
  --bib path/to/seed_articles.bib \
  --min-recall 0.90 \
  --max-hits 1000 \
  --custom-label "Manuscript query" \
  --custom-query '((H5N1[Title/Abstract] OR H5[Title/Abstract]) AND (Influenza[Title/Abstract]) AND (fatal*[Title/Abstract])) NOT (H7N9[Title/Abstract] OR H1N1[Title/Abstract] OR H1N1pdm09[Title/Abstract])' \
  --out output_benchmark
```

#### Using an NCBI API Key for Faster Execution

By default, NCBI E-utilities allows ~3 requests per second. Setting
`NCBI_API_KEY` (and `NCBI_EMAIL`) increases the rate limit to ~10
requests per second:

``` bash
export NCBI_EMAIL="researcher@university.edu"
export NCBI_API_KEY="your_ncbi_api_key_here"

Rscript run_query_optimizer.R \
  --bib path/to/seed_articles.bib \
  --out output
```

------------------------------------------------------------------------

## 2. YAML Configuration File (`--config`)

For reproducible pipelines, you can store term overrides, custom concept
blocks, and benchmark queries in a YAML file and pass
`--config config.yml` (or `config_path = "config.yml"` in R):

``` yaml
# Maximum Google Scholar query character limit (hard ceiling: 256)
max_gs_chars: 256

# Maximum number of AND-connected concept blocks (1..4)
max_blocks: 3

# Maximum OR-synonyms per concept block
max_terms_per_block: 5

# Terms that MUST appear in the primary concept block
required_terms:
  - h5n1

# Candidate terms to ignore during automatic extraction
banned_terms:
  - world health organization
  - united states

# Exclusion terms applied via PubMed NOT (...) and Google Scholar -term / -"phrase"
# Supports prefix wildcards (*) and title:-scoped exclusions (-intitle: in Google Scholar)
exclude_terms:
  - h7n9
  - h1n1*
  - vaccin*
  - title:ferret*

# Optional manual concept blocks (leave empty/commented for automatic PICO discovery)
# custom_blocks:
#   pathogen:
#     - h5n1
#     - h5
#   outcome_severity:
#     - fatal*
#     - mortality
#     - death*
#   population_host:
#     - cases
#     - human*

# Optional custom PubMed queries to benchmark alongside the Pareto frontier
custom_queries:
  - label: "Manuscript query"
    pubmed_query: '((H5N1[Title/Abstract] OR H5[Title/Abstract]) AND (Influenza[Title/Abstract]) AND (fatal*[Title/Abstract])) NOT (H7N9[Title/Abstract] OR H1N1[Title/Abstract] OR H1N1pdm09[Title/Abstract])'
```

------------------------------------------------------------------------

## 3. Programmatic Usage in R

### High-Level Entry Point: `optimize_search_queries()`

``` r

library(revlitquery)

res <- optimize_search_queries(
  bib_path   = "path/to/seed_articles.bib",
  output_dir = "output",
  min_recall = 0.90,
  max_hits   = 1000,
  max_evals  = 60,
  custom_queries = list(
    list(
      label = "Manuscript query",
      pubmed_query = paste0(
        "((H5N1[Title/Abstract] OR H5[Title/Abstract]) AND ",
        "(Influenza[Title/Abstract]) AND (fatal*[Title/Abstract])) ",
        "NOT (H7N9[Title/Abstract] OR H1N1[Title/Abstract] OR H1N1pdm09[Title/Abstract])"
      )
    )
  )
)

# Inspect the Pareto-optimal query table
res$pareto_df[, c("rank", "combined_recall", "pm_hits", "gs_chars", "pubmed_query", "google_scholar_query")]

# Inspect concept-grouped exclusion hints
head(res$exclusion_hints, 10)
```

### Step-by-Step Modular Pipeline (Offline Demo)

You can also call each stage of the `revlitquery` pipeline individually.
Below is a self-contained offline demonstration using a small synthetic
BibTeX corpus:

``` r

library(revlitquery)

# 1. Write and parse a 6-article demo BibTeX file
demo_bib <- file.path(tempdir(), "demo_seeds.bib")
revlitquery:::.write_demo_bibtex(demo_bib)
#> [1] "/tmp/RtmplhtVxH/demo_seeds.bib"

seeds_df <- parse_bibtex(demo_bib)
seeds_df[, c("citekey", "year", "pmid", "title")]
#> # A tibble: 6 × 4
#>   citekey                   year  pmid     title                                
#>   <chr>                     <chr> <chr>    <chr>                                
#> 1 li2008estimates           2008  18439345 Estimates of the case fatality rate …
#> 2 wang2012serological       2012  22362880 Serological evidence of human infect…
#> 3 fiebig2011epidemiological 2011  21679679 Epidemiological and clinical charact…
#> 4 kandun2008factors         2008  18959499 Factors associated with case fatalit…
#> 5 lai2016global             2016  27211899 Global epidemiology of avian influen…
#> 6 cowling2013comparative    2013  23803488 Comparative epidemiology of human in…
```

``` r

# 2. Prepare local seed text fields (offline; use enrich_seeds() for live PubMed enrichment)
seeds_local <- seeds_df
seeds_local$in_pubmed <- nzchar(seeds_local$pmid)
seeds_local$full_title <- seeds_local$title
seeds_local$full_abstract <- seeds_local$abstract
seeds_local$full_text <- tolower(paste(seeds_local$title, seeds_local$abstract))

# 3. Extract candidate n-grams and prefix wildcards (offline mode: probe_pubmed = FALSE)
terms_df <- extract_candidate_terms(
  seeds_df     = seeds_local,
  config       = list(),
  probe_pubmed = FALSE
)
#> Warning: Unknown or uninitialised column: `all_keywords`.
#> Warning: Unknown or uninitialised column: `mesh_terms`.
head(terms_df[, c("term", "n_words", "df_tiab", "df_title", "lift", "term_score")], 10)
#> # A tibble: 10 × 6
#>    term               n_words df_tiab df_title  lift term_score
#>    <chr>                <int>   <int>    <int> <dbl>      <dbl>
#>  1 avian influenza          2       5        5 1267.       4.01
#>  2 human h5n1               2       5        2 1267.       3.90
#>  3 h5n1                     1       6        6  190.       3.54
#>  4 case fatality            2       5        2 1267.       3.39
#>  5 human                    1       6        5  190.       2.94
#>  6 influenza                1       6        5  190.       2.94
#>  7 case fatality rate       3       3        1 9495.       2.55
#>  8 fatality                 1       6        2  190.       2.54
#>  9 avian                    1       5        5  158.       2.47
#> 10 case                     1       5        3  158.       2.22
```

``` r

# 4. Discover PICO / semantic concept blocks
blocks <- build_concept_blocks(
  terms_df = terms_df,
  seeds_df = seeds_local,
  config   = list(max_blocks = 3, max_terms_per_block = 4)
)
str(blocks)
#> List of 3
#>  $ :List of 2
#>   ..$ name : chr "primary_topic"
#>   ..$ terms: chr [1:2] "avian influenza" "influenza"
#>  $ :List of 2
#>   ..$ name : chr "outcome_severity"
#>   ..$ terms: chr [1:2] "case fatality" "fatality risk"
#>   ..- attr(*, "narrow_term")= chr "case fatality rate"
#>  $ :List of 2
#>   ..$ name : chr "population_host"
#>   ..$ terms: chr [1:2] "human" "confirmed human"
#>   ..- attr(*, "narrow_term")= chr "confirmed human"
#>  - attr(*, "block1_coupled")= chr [1:4] "avian influenza" "h5n1" "influenza" "avian"
#>  - attr(*, "joint_hits")= Named num(0) 
#>   ..- attr(*, "names")= chr(0) 
#>  - attr(*, "term_families")= Named chr [1:45] "primary_topic" "primary_topic" "primary_topic" "outcome_severity" ...
#>   ..- attr(*, "names")= chr [1:45] "avian influenza" "human h5n1" "h5n1" "case fatality" ...
#>  - attr(*, "domain_qualifier")= chr "h5n1"
```

``` r

# 5. Translate a structured query specification into PubMed & Google Scholar syntax
q <- make_query(
  blocks = list(
    list(name = "primary_topic", terms = c("h5n1", "h5"), field = "tiab"),
    list(name = "domain_context", terms = "influenza", field = "tiab"),
    list(name = "outcome_severity", terms = c("fatal*", "case fatality"), field = "tiab")
  ),
  exclude_terms = c("h7n9", "h1n1*", "vaccin*")
)

cat("PubMed query:\n", to_pubmed_query(q), "\n\n")
#> PubMed query:
#>  (("h5n1"[tiab] OR "h5"[tiab]) AND ("influenza"[tiab]) AND ("fatal*"[tiab])) NOT ("h7n9"[tiab] OR "h1n1*"[tiab] OR "vaccin*"[tiab])
cat("Google Scholar query:\n", to_google_scholar_query(q), "\n")
#> Google Scholar query:
#>  (h5n1 OR h5) influenza fatal -h7n9 -h1n1 -vaccine
```

------------------------------------------------------------------------

## 4. Output Files Generated in `output_dir`

When
[`optimize_search_queries()`](https://pascalcrepey.github.io/revlitquery/reference/optimize_search_queries.md)
finishes, it writes a self-contained report bundle to `output_dir`:

| File | Description |
|----|----|
| `pareto_queries.csv` | Non-dominated PubMed and Google Scholar query pairs, recall, hit counts, character lengths, missed citekeys, and clickable search URLs. |
| `all_evaluated_queries.csv` | Complete log of all candidate and custom queries evaluated across Stage 1 (grid) and Stage 2 (local + exclusion refinement). |
| `seed_resolution.csv` | Per-article PubMed resolution summary (`citekey`, `pmid`, `doi`, `in_pubmed`, abstract availability). |
| `seed_coverage_matrix.csv` | Paper-by-query matrix showing `TRUE`/`FALSE` retrieval of every seed article across each Pareto query (`Q1`, `Q2`, …). |
| `exclusion_hints.csv` | Up to 35 ranked `NOT` / `-` exclusion candidates (including prefix wildcards and `title:`-scoped terms) organized by exclusion concept (`competing_subtypes`, `experimental_lab`, `vaccines_immunology`, `veterinary_ecology`, `general_noise`) with verified PubMed hit reductions and seed loss counts. |
| `report.md` | Formatted Markdown summary report ready to share with co-authors or include in a systematic review appendix. |
| `pareto_frontier.png` | Publication-ready plot of the sensitivity (seed recall %) vs. specificity (total PubMed hits, log scale) Pareto frontier. |

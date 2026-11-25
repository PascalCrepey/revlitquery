# Evaluate a PubMed query: total hit count + which seed PMIDs are retrieved

Requests up to 2500 PMIDs in the initial \`esearch\` call so that
queries with \<= 2500 hits determine seed recall in a single HTTP
request instead of two.

## Usage

``` r
pubmed_eval_query(
  pm_query,
  seed_pmids,
  email = NULL,
  api_key = NULL,
  cache_dir = ".cache",
  use_cache = TRUE
)
```

## Arguments

- pm_query:

  PubMed query string.

- seed_pmids:

  Character vector of seed PMIDs.

- email:

  Optional email for NCBI tool identification.

- api_key:

  Optional NCBI API key.

- cache_dir:

  Directory for RDS cache.

- use_cache:

  Logical; whether to read/write cache.

## Value

List with \`pm_hits\`, \`pm_found\`, \`pm_recall\`, \`ok\`.

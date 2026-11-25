# Evaluate a single user-supplied PubMed query against enriched seeds

Evaluate a single user-supplied PubMed query against enriched seeds

## Usage

``` r
evaluate_custom_query(
  pm_query,
  seeds_df,
  label = "Custom",
  gs_query = NULL,
  email = NULL,
  ncbi_key = NULL,
  cache_dir = ".cache",
  use_cache = TRUE,
  eval_pubmed_fn = pubmed_eval_query
)
```

## Arguments

- pm_query:

  Raw PubMed query string.

- seeds_df:

  Enriched seeds tibble.

- label:

  Display label for the custom query.

- gs_query:

  Optional Google Scholar query string (auto-approximated if
  NULL/blank).

- email:

  Optional email for NCBI polite pool.

- ncbi_key:

  Optional NCBI API key.

- cache_dir:

  Cache directory path.

- use_cache:

  Logical; whether to read/write cache.

- eval_pubmed_fn:

  Function \`(pm_query, seed_pmids, ...)\` returning \`pm_hits\`,
  \`pm_found\`, \`pm_recall\`, \`ok\`.

## Value

1-row tibble matching \`eval_df\` schema.

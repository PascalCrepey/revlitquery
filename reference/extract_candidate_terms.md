# Extract and score candidate search terms from enriched seed papers

Extract and score candidate search terms from enriched seed papers

## Usage

``` r
extract_candidate_terms(
  seeds_df,
  config = list(),
  max_candidates = 50L,
  probe_pubmed = TRUE,
  eval_bg_fn = NULL,
  email = NULL,
  ncbi_key = NULL,
  cache_dir = ".cache",
  use_cache = TRUE,
  progress_cb = NULL
)
```

## Arguments

- seeds_df:

  Enriched tibble from \`enrich_seeds()\`.

- config:

  Optional config list (\`banned_terms\`, \`required_terms\`,
  \`custom_blocks\`).

- max_candidates:

  Max terms to probe against PubMed background frequency.

- probe_pubmed:

  Logical; whether to probe PubMed for background \`bg_hits\`.

- eval_bg_fn:

  Optional function \`(term) -\> integer\` for testing without network.

- email:

  Optional email for NCBI polite pool.

- ncbi_key:

  Optional NCBI API key.

- cache_dir:

  Cache directory path.

- use_cache:

  Logical; whether to read/write cache.

- progress_cb:

  Optional progress callback function.

## Value

Tibble of scored candidate terms with \`cov_tiab\` and \`cov_title\`
list-columns.

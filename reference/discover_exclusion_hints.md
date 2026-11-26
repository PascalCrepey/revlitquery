# Discover unsupervised, concept-clustered exclusion (\`NOT\` / \`-\`) terms from non-seed PubMed results

Samples articles returned by \`base_query\` that are NOT in \`seeds_df\`
across the full publication timeline, mines frequent 1- and 2-grams plus
morphological prefix wildcards (\`\*\`), checks seed loss locally
against PubMed-resolved seeds (including PubMed Author Keywords),
partitions candidates into unsupervised noise topic clusters via
\[cluster_exclusion_terms()\], builds greedy set-cover exclusion blocks,
and returns up to \`max_hints\` ranked exclusion hints.

## Usage

``` r
discover_exclusion_hints(
  base_query,
  seeds_df,
  base_pm_hits = NA_integer_,
  sample_n = 160L,
  top_verify = 12L,
  max_hints = 35L,
  email = NULL,
  ncbi_key = NULL,
  cache_dir = ".cache",
  use_cache = TRUE,
  pm_esearch_fn = pubmed_esearch,
  pm_efetch_fn = pubmed_efetch,
  progress_cb = NULL
)
```

## Arguments

- base_query:

  Canonical \`search_query\` object to analyze.

- seeds_df:

  Enriched seeds tibble.

- base_pm_hits:

  Integer total PubMed hits of \`base_query\`.

- sample_n:

  Maximum number of non-seed PubMed articles to sample.

- top_verify:

  Number of top exclusion candidates to verify live via \`esearch\`.

- max_hints:

  Maximum number of exclusion hints to return and display (default 35).

- email:

  Optional email for NCBI polite pool.

- ncbi_key:

  Optional NCBI API key.

- cache_dir:

  Cache directory path.

- use_cache:

  Logical; whether to read/write cache.

- pm_esearch_fn:

  Optional dependency injection for \[pubmed_esearch()\].

- pm_efetch_fn:

  Optional dependency injection for \[pubmed_efetch()\].

- progress_cb:

  Optional progress callback function.

## Value

Tibble of ranked exclusion hints with \`concept\` column and
\`exclusion_blocks\` attribute.

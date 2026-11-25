# Run Stage 1 + Stage 2 Pareto query optimization

Run Stage 1 + Stage 2 Pareto query optimization

## Usage

``` r
optimize_queries(
  seeds_df,
  terms_df,
  blocks,
  config = list(),
  min_recall = 0.6,
  max_hits = 2000,
  max_evals = 60L,
  max_gs_chars = 256L,
  discover_exclusions = TRUE,
  email = NULL,
  ncbi_key = NULL,
  cache_dir = ".cache",
  use_cache = TRUE,
  eval_pubmed_fn = pubmed_eval_query,
  pm_esearch_fn = pubmed_esearch,
  pm_efetch_fn = pubmed_efetch,
  progress_cb = NULL
)
```

## Arguments

- seeds_df:

  Enriched seeds tibble.

- terms_df:

  Scored candidate terms tibble.

- blocks:

  Concept blocks list.

- config:

  Optional config list.

- min_recall:

  Minimum desired seed recall (0..1).

- max_hits:

  Optional maximum PubMed hit count threshold.

- max_evals:

  Maximum live API query evaluations.

- max_gs_chars:

  Maximum Google Scholar character limit (default 256).

- discover_exclusions:

  Logical; whether to mine non-seed PubMed results for exclusion hints.

- eval_pubmed_fn:

  Function \`(pm_query, seed_pmids, ...)\` returning \`pm_hits\`,
  \`pm_found\`, \`pm_recall\`, \`ok\`.

## Value

Optimization result list containing \`eval_df\`, \`pareto_df\`,
\`exclusion_hints\`, \`recommended_sig\`, \`seeds_df\`, \`blocks\`.

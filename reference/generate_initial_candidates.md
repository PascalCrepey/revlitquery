# Generate diverse initial Stage 1 query candidates

Generate diverse initial Stage 1 query candidates

## Usage

``` r
generate_initial_candidates(
  blocks,
  terms_df,
  seeds_df,
  config = list(),
  min_recall = 0.6,
  max_hits = 2000,
  max_gs_chars = 256L,
  max_candidates = 26L
)
```

## Arguments

- blocks:

  List of discovered concept blocks.

- terms_df:

  Scored candidate terms tibble.

- seeds_df:

  Enriched seeds tibble.

- config:

  Optional config list.

- min_recall:

  Target minimum recall used to guide diversity selection.

- max_hits:

  Optional maximum PubMed hit count cap.

- max_gs_chars:

  Maximum Google Scholar query characters (default 256).

- max_candidates:

  Maximum initial candidates to return.

## Value

List of \`search_query\` objects.

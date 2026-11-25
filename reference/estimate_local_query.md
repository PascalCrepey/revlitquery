# Estimate local seed coverage and rough PubMed hit count for a candidate query

Estimate local seed coverage and rough PubMed hit count for a candidate
query

## Usage

``` r
estimate_local_query(q, terms_df, seeds_df)
```

## Arguments

- q:

  Canonical \`search_query\`.

- terms_df:

  Scored candidate terms tibble.

- seeds_df:

  Enriched seeds tibble.

## Value

List with \`local_cov\`, \`local_found\`, \`local_recall\`,
\`est_hits\`.

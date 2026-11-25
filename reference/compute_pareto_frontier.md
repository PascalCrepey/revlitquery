# Compute non-dominated Pareto frontier indicator

Maximizes \`combined_recall\`, minimizes \`pm_hits\`. When \`max_hits\`
is finite, restricts the frontier to queries with \`pm_hits \<=
max_hits\` (falling back to all valid queries only if none satisfy the
cap).

## Usage

``` r
compute_pareto_frontier(eval_df, max_hits = Inf)
```

## Arguments

- eval_df:

  Data frame / tibble of evaluated queries.

- max_hits:

  Optional maximum PubMed hit count cap.

## Value

Logical vector of length \`nrow(eval_df)\`.

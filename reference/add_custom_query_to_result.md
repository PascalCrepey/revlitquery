# Merge an evaluated custom query into an existing or minimal optimization result

Merge an evaluated custom query into an existing or minimal optimization
result

## Usage

``` r
add_custom_query_to_result(result, custom_row)
```

## Arguments

- result:

  Optimization result list from \`optimize_queries()\` (or minimal
  shell).

- custom_row:

  1-row tibble returned by \`evaluate_custom_query()\`.

## Value

Updated optimization result list with \`eval_df\` and \`pareto_df\`
recomputed.

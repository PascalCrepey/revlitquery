# Build display table for the Seed Coverage Matrix tab

Build display table for the Seed Coverage Matrix tab

## Usage

``` r
build_coverage_display_df(seeds_df, pareto_df, eval_df = NULL)
```

## Arguments

- seeds_df:

  Enriched seeds tibble.

- pareto_df:

  Pareto frontier tibble.

- eval_df:

  Optional all evaluated queries tibble (for custom queries).

## Value

Tibble ready for \`DT::datatable()\`.

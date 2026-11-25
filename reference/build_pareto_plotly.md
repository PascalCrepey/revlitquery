# Build interactive Plotly Pareto frontier chart

Build interactive Plotly Pareto frontier chart

## Usage

``` r
build_pareto_plotly(
  eval_df,
  pareto_df,
  selected_rank = NULL,
  min_recall = 0.6,
  max_hits = Inf,
  selected_eval_id = NULL
)
```

## Arguments

- eval_df:

  All evaluated queries tibble.

- pareto_df:

  Pareto-optimal queries tibble.

- selected_rank:

  Currently inspected Pareto rank (\`Q1\`, \`Q2\`, ...).

- min_recall:

  Target minimum recall (\`0..1\`).

- max_hits:

  Optional max PubMed hit threshold (\`Inf\` if none).

- selected_eval_id:

  Currently inspected \`eval_id\`.

## Value

A \`plotly\` htmlwidget.

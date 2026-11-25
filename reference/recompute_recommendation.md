# Recompute recommended query on an existing Pareto frontier when sliders move

Recompute recommended query on an existing Pareto frontier when sliders
move

## Usage

``` r
recompute_recommendation(pareto_df, min_recall = 0.6, max_hits = Inf)
```

## Arguments

- pareto_df:

  Pareto frontier tibble.

- min_recall:

  Minimum seed recall in \`0..1\`.

- max_hits:

  Maximum PubMed hit threshold (\`Inf\` or numeric).

## Value

List with updated \`pareto_df\`, \`recommended_sig\`, and
\`recommended_status\`.

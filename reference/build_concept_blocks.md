# Build stem-disjoint, conditionally selective concept blocks from scored candidate terms

Build stem-disjoint, conditionally selective concept blocks from scored
candidate terms

## Usage

``` r
build_concept_blocks(
  terms_df,
  seeds_df,
  config = list(),
  max_blocks = NULL,
  max_terms_per_block = NULL,
  probe_selectivity = NULL,
  eval_joint_fn = NULL
)
```

## Arguments

- terms_df:

  Scored candidate tibble from \`extract_candidate_terms()\`.

- seeds_df:

  Enriched seeds tibble.

- config:

  Optional config list (\`custom_blocks\`, \`max_blocks\`,
  \`max_terms_per_block\`).

- max_blocks:

  Maximum number of concept blocks.

- max_terms_per_block:

  Maximum terms per block.

- probe_selectivity:

  Logical or NULL (defaults to TRUE when \`terms_df\` was probed live).

- eval_joint_fn:

  Optional function \`(t1, t2) -\> integer\` returning PubMed joint hit
  count.

## Value

List of block lists \`list(name = ..., terms = ...)\`.

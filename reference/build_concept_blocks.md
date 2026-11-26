# Build orthogonal, unsupervised concept blocks from scored candidate terms

Discovers Boolean \`AND\` concept blocks and their \`OR\` synonyms
directly from the seed papers and PubMed co-occurrence statistics
without hardcoded word lists or domain-specific ontologies: 1.
Identifies the primary high-lift anchor concept (\`primary_topic\`) and
expands stem/code synonyms covering remaining seeds. 2. Separates any
universal domain-coupled context term (\`domain_context\`). 3.
Partitions remaining selective terms via \[discover_term_clusters()\]
(SVD local context + MeSH/keyword co-assignment + coverage
complementarity Louvain graph) and grows each concept block via greedy
minimum-hit seed set-cover.

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

  Scored candidate tibble from \[extract_candidate_terms()\].

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

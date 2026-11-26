# Unsupervised concept cluster discovery across candidate search terms

Partitions candidate terms into coherent conceptual facets using a
hybrid affinity graph combining: 1. Second-order distributional context
& MeSH/keyword SVD similarity (\`LSA\`), 2. Seed coverage
complementarity (linking substitutable \`OR\` synonyms that cover each
other's missed seed papers while keeping high-frequency orthogonal
\`AND\` pillars separate), and 3. Lexical stem, code-family, and
character n-gram similarity, followed by Louvain community detection
(\`igraph::cluster_louvain()\`).

## Usage

``` r
discover_term_clusters(terms_df, seeds_df, exclude_terms = character(0))
```

## Arguments

- terms_df:

  Scored candidate terms tibble from \[extract_candidate_terms()\].

- seeds_df:

  Enriched seeds tibble.

- exclude_terms:

  Optional character vector of terms (e.g., Block 1 coupled terms) to
  assign directly to \`"primary_topic"\` before clustering remaining
  facets.

## Value

Named character vector mapping each term in \`terms_df\$term\` to its
automatically discovered cluster label (e.g., \`"concept_fatal"\`,
\`"concept_cases"\`).

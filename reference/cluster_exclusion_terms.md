# Unsupervised co-occurrence graph clustering of candidate exclusion terms

Partitions candidate exclusion terms into data-driven noise topic
clusters based on their empirical co-occurrence across sampled non-seed
articles (combining Normalized Pointwise Mutual Information, Jaccard
overlap, shared MeSH/Author Keyword profiles, stem similarity, and
alphanumeric code pattern affinity), using Louvain modularity
maximization (\[igraph::cluster_louvain()\]).

## Usage

``` r
cluster_exclusion_terms(
  terms,
  noise_cov_list,
  rank_scores = NULL,
  noise_mesh = NULL,
  noise_keywords = NULL
)
```

## Arguments

- terms:

  Character vector of candidate exclusion terms.

- noise_cov_list:

  List of logical vectors (each of length \`n_noise\`) indicating which
  sampled non-seed articles contain each candidate term.

- rank_scores:

  Numeric vector of candidate importance scores (same length as
  \`terms\`).

- noise_mesh:

  Optional list of character vectors of MeSH headings per noise article.

- noise_keywords:

  Optional list of character vectors of author keywords per noise
  article.

## Value

Character vector of automatically generated concept cluster labels (same
length as \`terms\`, e.g., \`"topic_mice_receptor"\`,
\`"topic_h1n1_h7n9"\`).

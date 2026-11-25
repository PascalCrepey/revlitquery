# Fetch article metadata (title, abstract, MeSH, keywords) from PubMed efetch XML

Uses both batch-level and per-PMID disk caching so previously fetched
articles never need to be re-downloaded when seed lists or noise samples
change.

## Usage

``` r
pubmed_efetch(
  pmids,
  email = NULL,
  api_key = NULL,
  cache_dir = ".cache",
  use_cache = TRUE
)
```

## Arguments

- pmids:

  Character vector of PMIDs.

## Value

Tibble with \`pmid\`, \`pm_title\`, \`pm_abstract\`, \`pm_mesh\`,
\`pm_keywords\`.

# Enrich BibTeX seed papers via PubMed

Enrich BibTeX seed papers via PubMed

## Usage

``` r
enrich_seeds(
  seeds_df,
  email = NULL,
  ncbi_key = NULL,
  cache_dir = ".cache",
  use_cache = TRUE,
  pm_esearch_fn = pubmed_esearch,
  pm_efetch_fn = pubmed_efetch,
  progress_cb = NULL
)
```

## Arguments

- seeds_df:

  Tibble returned by \`parse_bibtex()\`.

- email:

  Optional email for NCBI polite pool.

- ncbi_key:

  Optional NCBI API key.

- cache_dir:

  Cache directory path.

- use_cache:

  Logical.

- pm_esearch_fn:

  Optional dependency injection for \`pubmed_esearch\`.

- pm_efetch_fn:

  Optional dependency injection for \`pubmed_efetch\`.

- progress_cb:

  Optional progress callback function.

## Value

Enriched tibble with \`full_title\`, \`full_abstract\`, \`full_text\`,
\`pmid\`, \`in_pubmed\`, \`mesh_terms\`, \`all_keywords\`.

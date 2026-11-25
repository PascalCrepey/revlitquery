# Query PubMed E-utilities esearch via HTTP POST

Query PubMed E-utilities esearch via HTTP POST

## Usage

``` r
pubmed_esearch(
  term,
  retmax = 0L,
  email = NULL,
  api_key = NULL,
  cache_dir = ".cache",
  use_cache = TRUE
)
```

## Arguments

- term:

  PubMed query string.

- retmax:

  Maximum number of PMIDs to return (0 for count only).

- email:

  Optional email for NCBI tool identification.

- api_key:

  Optional NCBI API key.

- cache_dir:

  Directory for RDS cache.

- use_cache:

  Logical; whether to read/write cache.

## Value

List with \`count\` (integer), \`idlist\` (character), \`warnings\`
(character), \`ok\` (logical).

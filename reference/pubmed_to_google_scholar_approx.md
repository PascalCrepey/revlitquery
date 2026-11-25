# Convert a raw PubMed query string into an approximate Google Scholar query

Preserves boolean OR/AND/NOT semantics, converts
\`\[ti\]\`/\`\[Title\]\` tags to \`intitle:\`, strips
\`\[tiab\]\`/\`\[MeSH\]\` tags and \`\*\` wildcards, converts \`AND\` to
implicit spaces, and formats \`NOT (...)\` exclusions as \`-term\`.

## Usage

``` r
pubmed_to_google_scholar_approx(pm_query)
```

## Arguments

- pm_query:

  Raw PubMed query string.

## Value

Approximate Google Scholar query string.

## Examples

``` r
pubmed_to_google_scholar_approx(
  '(("H5N1"[tiab] OR "H5"[tiab]) AND "fatal*"[tiab]) NOT "H7N9"[tiab]'
)
#> [1] "(H5N1 OR H5) fatal -H7N9"
```

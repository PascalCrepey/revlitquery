# Parse a BibTeX file into a normalized tibble of seed papers

Parse a BibTeX file into a normalized tibble of seed papers

## Usage

``` r
parse_bibtex(path)
```

## Arguments

- path:

  Path to a \`.bib\` file.

## Value

A tibble with columns \`seed_id\`, \`citekey\`, \`entry_type\`,
\`title\`, \`abstract\`, \`doi\`, \`pmid\`, \`keywords\`, \`year\`,
\`author\`, \`journal\`.

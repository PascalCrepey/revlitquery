# Translate canonical query to Google Scholar query string

Uses implicit AND (space) between concept blocks, uppercase OR,
\`intitle:\` for title-restricted blocks, and unquoted clean single
tokens to stay within 256 chars.

## Usage

``` r
to_google_scholar_query(q)
```

## Arguments

- q:

  Canonical \`search_query\` object.

## Value

Character scalar.

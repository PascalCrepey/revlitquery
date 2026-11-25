# Construct a canonical query object

Construct a canonical query object

## Usage

``` r
make_query(blocks, exclude_terms = character(0))
```

## Arguments

- blocks:

  List of blocks, each \`list(name = ..., field = "tiab"\|"title", terms
  = character())\`.

- exclude_terms:

  Optional character vector of NOT terms.

## Value

Canonical query list of class \`"search_query"\`.

## Examples

``` r
q <- make_query(
  blocks = list(
    list(name = "topic", field = "tiab", terms = c("h5n1", "h5")),
    list(name = "outcome", field = "tiab", terms = "fatal*")
  ),
  exclude_terms = c("h7n9", "vaccin*")
)
to_pubmed_query(q)
#> [1] "((\"h5n1\"[tiab] OR \"h5\"[tiab]) AND (\"fatal*\"[tiab])) NOT (\"h7n9\"[tiab] OR \"vaccin*\"[tiab])"
to_google_scholar_query(q)
#> [1] "(h5n1 OR h5) fatal -h7n9 -vaccine"
```

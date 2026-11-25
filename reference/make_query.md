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

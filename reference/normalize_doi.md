# Normalize a DOI string to lowercase bare identifier

Normalize a DOI string to lowercase bare identifier

## Usage

``` r
normalize_doi(x)
```

## Arguments

- x:

  Character scalar or vector.

## Value

Normalized DOI string (empty string if missing).

## Examples

``` r
normalize_doi("https://doi.org/10.1056/NEJMra052211")
#> [1] "10.1056/nejmra052211"
```

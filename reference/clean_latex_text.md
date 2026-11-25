# Clean LaTeX markup and normalize whitespace in BibTeX field text

Clean LaTeX markup and normalize whitespace in BibTeX field text

## Usage

``` r
clean_latex_text(x)
```

## Arguments

- x:

  Character scalar or vector.

## Value

Cleaned character vector.

## Examples

``` r
clean_latex_text("Human {H5N1} avian \\emph{influenza} infections")
#> [1] "Human H5N1 avian influenza infections"
```

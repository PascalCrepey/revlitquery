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

## Examples

``` r
tmp_bib <- tempfile(fileext = ".bib")
writeLines(c(
  "@article{li2008,",
  "  title = {Case fatality rate of human H5N1 influenza},",
  "  author = {Li, F. C.},",
  "  year = {2008},",
  "  pmid = {18439345}",
  "}"
), tmp_bib)
parse_bibtex(tmp_bib)
#> # A tibble: 1 × 11
#>   seed_id citekey entry_type title    abstract doi   pmid  keywords year  author
#>     <int> <chr>   <chr>      <chr>    <chr>    <chr> <chr> <chr>    <chr> <chr> 
#> 1       1 li2008  article    Case fa… ""       ""    1843… ""       2008  Li, F…
#> # ℹ 1 more variable: journal <chr>
unlink(tmp_bib)
```

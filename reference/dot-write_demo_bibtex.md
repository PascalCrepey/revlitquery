# Write a built-in 6-entry synthetic H5N1 demo BibTeX file to a temporary path

Used by the interactive Shiny app's "Load sample H5N1 .bib" button and
unit tests so the package does not bundle any external \`.bib\` files or
cached data.

## Usage

``` r
.write_demo_bibtex(dest_path = file.path(tempdir(), "sample_h5n1_cfr.bib"))
```

## Arguments

- dest_path:

  Destination file path (defaults to \`file.path(tempdir(),
  "sample_h5n1_cfr.bib")\`).

## Value

Character path to the written \`.bib\` file.

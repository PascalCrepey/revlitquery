# Create the revlitquery Shiny application object

Create the revlitquery Shiny application object

## Usage

``` r
revlitquery_app(initial_result = NULL, initial_bib_label = NULL)
```

## Arguments

- initial_result:

  Optional precomputed result list from \[optimize_search_queries()\] to
  display immediately on startup.

- initial_bib_label:

  Optional display label for \`initial_result\`.

## Value

A \`shiny.appobj\` object that can be passed to \[shiny::runApp()\].

## Examples

``` r
app <- revlitquery_app()
```

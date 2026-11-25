# Launch the interactive revlitquery Shiny application

Launch the interactive revlitquery Shiny application

## Usage

``` r
run_revlitquery_app(
  initial_result = NULL,
  initial_bib_label = NULL,
  launch.browser = TRUE,
  ...
)
```

## Arguments

- initial_result:

  Optional precomputed result list from \[optimize_search_queries()\] to
  display immediately on startup.

- initial_bib_label:

  Optional display label for \`initial_result\`.

- launch.browser:

  Logical; open the application in the default browser (default
  \`TRUE\`).

- ...:

  Additional arguments passed to \[shiny::runApp()\].

## Value

Invisibly returns the result of \[shiny::runApp()\].

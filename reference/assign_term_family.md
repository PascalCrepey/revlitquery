# Classify a candidate search term into a coherent semantic PICO family

Classify a candidate search term into a coherent semantic PICO family

## Usage

``` r
assign_term_family(term, block1_coupled = character(0))
```

## Arguments

- term:

  Single candidate term.

- block1_coupled:

  Optional character vector of terms coupled to Block 1.

## Value

One of \`"primary_topic"\`, \`"outcome_severity"\`,
\`"population_host"\`, \`"epi_methods"\`, or \`"other"\`.

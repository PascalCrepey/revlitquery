# Validate and assemble optimizer config list from Shiny sidebar inputs

Validate and assemble optimizer config list from Shiny sidebar inputs

## Usage

``` r
build_config_from_inputs(
  required_txt = "",
  banned_txt = "",
  exclude_txt = "",
  custom_blocks_yaml = "",
  max_blocks = 3L,
  max_terms_per_block = 5L,
  max_gs_chars = 256L
)
```

## Value

List with \`ok\`, \`config\`, and \`error\`.

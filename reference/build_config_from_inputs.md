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

## Arguments

- required_txt:

  Comma/newline-separated required terms string.

- banned_txt:

  Comma/newline-separated banned terms string.

- exclude_txt:

  Comma/newline-separated exclusion terms string.

- custom_blocks_yaml:

  Optional YAML text defining custom concept blocks.

- max_blocks:

  Maximum number of concept blocks.

- max_terms_per_block:

  Maximum number of terms per concept block.

- max_gs_chars:

  Maximum Google Scholar query character limit.

## Value

List with \`ok\`, \`config\`, and \`error\`.

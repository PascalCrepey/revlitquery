# Optimize semantically equivalent PubMed and Google Scholar queries for a BibTeX file

End-to-end reverse literature review query optimizer: parses a \`.bib\`
file of seed papers, enriches metadata via PubMed E-utilities, extracts
high-lift domain terms and concept blocks, evaluates candidate query
combinations along the sensitivity (seed recall) vs. specificity (total
PubMed hits) Pareto frontier, mines non-seed PubMed results for
concept-grouped exclusion (\`NOT\` / \`-\`) terms and prefix wildcards,
and writes a report bundle.

## Usage

``` r
optimize_search_queries(
  bib_path,
  output_dir = "output",
  config_path = NULL,
  min_recall = 0.6,
  max_hits = 2000,
  max_evals = 60L,
  max_gs_chars = 256L,
  use_cache = TRUE,
  cache_dir = NULL,
  custom_queries = NULL,
  email = Sys.getenv("NCBI_EMAIL", ""),
  ncbi_key = Sys.getenv("NCBI_API_KEY", "")
)
```

## Arguments

- bib_path:

  Path to input \`.bib\` file of seed papers.

- output_dir:

  Output directory for reports, CSVs, and Pareto plot (or \`NULL\` to
  skip disk export).

- config_path:

  Optional YAML config file path.

- min_recall:

  Minimum desired seed recall (0..1), default 0.60.

- max_hits:

  Optional maximum PubMed hit count threshold, default 2000.

- max_evals:

  Maximum live API query evaluations, default 60.

- max_gs_chars:

  Maximum Google Scholar query length, default 256.

- use_cache:

  Logical; cache API responses under \`cache_dir\`, default TRUE.

- cache_dir:

  Optional cache directory path (defaults to \`file.path(output_dir,
  ".cache")\`).

- custom_queries:

  Optional character vector or list of custom PubMed queries.

- email:

  Optional email for NCBI polite pool.

- ncbi_key:

  Optional NCBI API key.

## Value

Optimization result list (invisibly).

## Examples

``` r
if (interactive()) {
  demo_bib <- tempfile(fileext = ".bib")
  writeLines(c(
    "@article{li2008,",
    "  title = {Case fatality rate of human H5N1 avian influenza},",
    "  year = {2008},",
    "  pmid = {18439345}",
    "}"
  ), demo_bib)
  res <- optimize_search_queries(demo_bib, output_dir = tempdir(), max_evals = 10)
}
```

# Changelog

## revlitquery 0.1.0

- Initial release of `revlitquery`: Reverse Literature Review Search
  Query Optimizer for PubMed and Google Scholar.
- Added BibTeX parser
  ([`parse_bibtex()`](https://pascalcrepey.github.io/revlitquery/reference/parse_bibtex.md))
  and PubMed seed metadata enrichment
  ([`enrich_seeds()`](https://pascalcrepey.github.io/revlitquery/reference/enrich_seeds.md))
  with year-assisted title matching, Author Keyword indexing, and
  per-PMID disk caching.
- Added candidate term extraction
  ([`extract_candidate_terms()`](https://pascalcrepey.github.io/revlitquery/reference/extract_candidate_terms.md)),
  unsupervised morphological prefix wildcard synthesis (`*`), and
  domain-agnostic Hybrid Contrastive Graph + Orthogonal Set-Cover
  concept block discovery
  ([`discover_term_clusters()`](https://pascalcrepey.github.io/revlitquery/reference/discover_term_clusters.md),
  [`build_concept_blocks()`](https://pascalcrepey.github.io/revlitquery/reference/build_concept_blocks.md)).
- Added canonical boolean query representation
  ([`make_query()`](https://pascalcrepey.github.io/revlitquery/reference/make_query.md))
  and translators for PubMed
  ([`to_pubmed_query()`](https://pascalcrepey.github.io/revlitquery/reference/to_pubmed_query.md))
  and Google Scholar (`<= 256` characters;
  [`to_google_scholar_query()`](https://pascalcrepey.github.io/revlitquery/reference/to_google_scholar_query.md),
  [`pubmed_to_google_scholar_approx()`](https://pascalcrepey.github.io/revlitquery/reference/pubmed_to_google_scholar_approx.md)).
- Added two-stage sensitivity-specificity Pareto frontier optimizer
  ([`optimize_queries()`](https://pascalcrepey.github.io/revlitquery/reference/optimize_queries.md),
  [`optimize_search_queries()`](https://pascalcrepey.github.io/revlitquery/reference/optimize_search_queries.md))
  with timeline-stratified, unsupervised NPMI co-occurrence noise topic
  clustering
  ([`cluster_exclusion_terms()`](https://pascalcrepey.github.io/revlitquery/reference/cluster_exclusion_terms.md),
  [`discover_exclusion_hints()`](https://pascalcrepey.github.io/revlitquery/reference/discover_exclusion_hints.md)).
- Added interactive Shiny dashboard
  ([`revlitquery_app()`](https://pascalcrepey.github.io/revlitquery/reference/revlitquery_app.md),
  [`run_revlitquery_app()`](https://pascalcrepey.github.io/revlitquery/reference/run_revlitquery_app.md))
  and two package vignettes (`shiny-interface` and
  `cli-and-r-workflow`).

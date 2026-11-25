# revlitquery 0.1.0

* Initial release of `revlitquery`: Reverse Literature Review Search Query Optimizer for PubMed and Google Scholar.
* Added BibTeX parser (`parse_bibtex()`) and PubMed seed metadata enrichment (`enrich_seeds()`) with year-assisted title matching, Author Keyword indexing, and per-PMID disk caching.
* Added candidate term extraction (`extract_candidate_terms()`), morphological prefix wildcard synthesis (`*`), and PICO concept block discovery (`build_concept_blocks()`).
* Added canonical boolean query representation (`make_query()`) and translators for PubMed (`to_pubmed_query()`) and Google Scholar (`<= 256` characters; `to_google_scholar_query()`, `pubmed_to_google_scholar_approx()`).
* Added two-stage sensitivity-specificity Pareto frontier optimizer (`optimize_queries()`, `optimize_search_queries()`) with timeline-stratified, concept-grouped exclusion mining (`discover_exclusion_hints()`).
* Added interactive Shiny dashboard (`revlitquery_app()`, `run_revlitquery_app()`) and two package vignettes (`shiny-interface` and `cli-and-r-workflow`).

# Package index

## High-Level Optimizer & Shiny Application

Primary entry points to run the query optimizer or launch the
interactive Shiny dashboard.

- [`optimize_search_queries()`](https://pascalcrepey.github.io/revlitquery/reference/optimize_search_queries.md)
  : Optimize semantically equivalent PubMed and Google Scholar queries
  for a BibTeX file
- [`run_revlitquery_app()`](https://pascalcrepey.github.io/revlitquery/reference/run_revlitquery_app.md)
  : Launch the interactive revlitquery Shiny application
- [`revlitquery_app()`](https://pascalcrepey.github.io/revlitquery/reference/revlitquery_app.md)
  : Create the revlitquery Shiny application object

## BibTeX Parsing & PubMed Seed Enrichment

Parse BibTeX files and resolve seed articles against NCBI PubMed
E-utilities.

- [`parse_bibtex()`](https://pascalcrepey.github.io/revlitquery/reference/parse_bibtex.md)
  : Parse a BibTeX file into a normalized tibble of seed papers
- [`enrich_seeds()`](https://pascalcrepey.github.io/revlitquery/reference/enrich_seeds.md)
  : Enrich BibTeX seed papers via PubMed
- [`clean_latex_text()`](https://pascalcrepey.github.io/revlitquery/reference/clean_latex_text.md)
  : Clean LaTeX markup and normalize whitespace in BibTeX field text
- [`normalize_doi()`](https://pascalcrepey.github.io/revlitquery/reference/normalize_doi.md)
  : Normalize a DOI string to lowercase bare identifier
- [`extract_pmid_value()`](https://pascalcrepey.github.io/revlitquery/reference/extract_pmid_value.md)
  : Extract a PubMed ID (PMID) from explicit fields or URL/note text
- [`split_keywords()`](https://pascalcrepey.github.io/revlitquery/reference/split_keywords.md)
  : Split BibTeX keywords string into cleaned lowercase character vector
- [`title_jaccard()`](https://pascalcrepey.github.io/revlitquery/reference/title_jaccard.md)
  : Compute alphanumeric lowercase token Jaccard similarity between two
  titles
- [`clean_title_for_query()`](https://pascalcrepey.github.io/revlitquery/reference/clean_title_for_query.md)
  : Clean title string for search queries (strip quotes/brackets/colons
  that break E-utilities)
- [`title_sig_words()`](https://pascalcrepey.github.io/revlitquery/reference/title_sig_words.md)
  : First N significant words from a title for fallback matching

## Candidate Term Extraction & Concept Blocks

Extract n-grams, synthesize prefix wildcards, probe PubMed frequencies,
and group terms into PICO concept blocks.

- [`extract_candidate_terms()`](https://pascalcrepey.github.io/revlitquery/reference/extract_candidate_terms.md)
  : Extract and score candidate search terms from enriched seed papers
- [`build_concept_blocks()`](https://pascalcrepey.github.io/revlitquery/reference/build_concept_blocks.md)
  : Build stem-disjoint, conditionally selective concept blocks from
  scored candidate terms
- [`assign_term_family()`](https://pascalcrepey.github.io/revlitquery/reference/assign_term_family.md)
  : Classify a candidate search term into a coherent semantic PICO
  family
- [`assign_term_families()`](https://pascalcrepey.github.io/revlitquery/reference/assign_term_families.md)
  : Vectorized wrapper for \`assign_term_family()\`
- [`term_matches_text()`](https://pascalcrepey.github.io/revlitquery/reference/term_matches_text.md)
  : Check whether a candidate term matches each text in \`texts\`
  (case-insensitive, hyphen-flexible)
- [`term_stems()`](https://pascalcrepey.github.io/revlitquery/reference/term_stems.md)
  : Simple token stemmer for lexical overlap detection across concept
  blocks
- [`shares_stem()`](https://pascalcrepey.github.io/revlitquery/reference/shares_stem.md)
  : Check whether two terms share any significant token stem
- [`normalize_search_text()`](https://pascalcrepey.github.io/revlitquery/reference/normalize_search_text.md)
  : Normalize text for n-gram extraction (preserves alphanumeric tokens
  like h5n1)
- [`escape_regex()`](https://pascalcrepey.github.io/revlitquery/reference/escape_regex.md)
  : Escape regex metacharacters in literal term
- [`get_stopwords()`](https://pascalcrepey.github.io/revlitquery/reference/get_stopwords.md)
  : Stopword list (standard English + academic/methods filler)

## Query Construction & Translation (PubMed & Google Scholar)

Canonical search query representation and translation to PubMed and
Google Scholar (\<= 256 chars).

- [`make_query()`](https://pascalcrepey.github.io/revlitquery/reference/make_query.md)
  : Construct a canonical query object
- [`query_signature()`](https://pascalcrepey.github.io/revlitquery/reference/query_signature.md)
  : Deterministic canonical signature for deduplicating queries
- [`to_pubmed_query()`](https://pascalcrepey.github.io/revlitquery/reference/to_pubmed_query.md)
  : Translate canonical query to PubMed E-utilities query string
- [`to_google_scholar_query()`](https://pascalcrepey.github.io/revlitquery/reference/to_google_scholar_query.md)
  : Translate canonical query to Google Scholar query string
- [`pubmed_to_google_scholar_approx()`](https://pascalcrepey.github.io/revlitquery/reference/pubmed_to_google_scholar_approx.md)
  : Convert a raw PubMed query string into an approximate Google Scholar
  query
- [`gs_char_count()`](https://pascalcrepey.github.io/revlitquery/reference/gs_char_count.md)
  : Count characters in translated Google Scholar query
- [`is_valid_gs_query()`](https://pascalcrepey.github.io/revlitquery/reference/is_valid_gs_query.md)
  : Check Google Scholar \<= 256 character constraint
- [`pubmed_url()`](https://pascalcrepey.github.io/revlitquery/reference/pubmed_url.md)
  : Clickable PubMed search URL
- [`google_scholar_url()`](https://pascalcrepey.github.io/revlitquery/reference/google_scholar_url.md)
  : Clickable Google Scholar search URL

## Pareto Frontier Optimization & Exclusion Mining

Two-stage Pareto frontier search, timeline-stratified exclusion concept
discovery, and custom query benchmarking.

- [`optimize_queries()`](https://pascalcrepey.github.io/revlitquery/reference/optimize_queries.md)
  : Run Stage 1 + Stage 2 Pareto query optimization
- [`discover_exclusion_hints()`](https://pascalcrepey.github.io/revlitquery/reference/discover_exclusion_hints.md)
  : Discover concept-grouped candidate exclusion (NOT / -) terms from
  non-seed PubMed results
- [`assign_exclusion_concept()`](https://pascalcrepey.github.io/revlitquery/reference/assign_exclusion_concept.md)
  : Classify a candidate exclusion term into a scientific exclusion
  concept
- [`compute_pareto_frontier()`](https://pascalcrepey.github.io/revlitquery/reference/compute_pareto_frontier.md)
  : Compute non-dominated Pareto frontier indicator
- [`generate_initial_candidates()`](https://pascalcrepey.github.io/revlitquery/reference/generate_initial_candidates.md)
  : Generate diverse initial Stage 1 query candidates
- [`estimate_local_query()`](https://pascalcrepey.github.io/revlitquery/reference/estimate_local_query.md)
  : Estimate local seed coverage and rough PubMed hit count for a
  candidate query
- [`evaluate_custom_query()`](https://pascalcrepey.github.io/revlitquery/reference/evaluate_custom_query.md)
  : Evaluate a single user-supplied PubMed query against enriched seeds
- [`add_custom_query_to_result()`](https://pascalcrepey.github.io/revlitquery/reference/add_custom_query_to_result.md)
  : Merge an evaluated custom query into an existing or minimal
  optimization result

## NCBI E-utilities & Caching

Low-level PubMed API clients and disk caching helpers.

- [`pubmed_esearch()`](https://pascalcrepey.github.io/revlitquery/reference/pubmed_esearch.md)
  : Query PubMed E-utilities esearch via HTTP POST
- [`pubmed_efetch()`](https://pascalcrepey.github.io/revlitquery/reference/pubmed_efetch.md)
  : Fetch article metadata (title, abstract, MeSH, keywords) from PubMed
  efetch XML
- [`pubmed_eval_query()`](https://pascalcrepey.github.io/revlitquery/reference/pubmed_eval_query.md)
  : Evaluate a PubMed query: total hit count + which seed PMIDs are
  retrieved
- [`cache_key()`](https://pascalcrepey.github.io/revlitquery/reference/cache_key.md)
  : Build deterministic cache key
- [`cache_get()`](https://pascalcrepey.github.io/revlitquery/reference/cache_get.md)
  : Read cached RDS value if available
- [`cache_set()`](https://pascalcrepey.github.io/revlitquery/reference/cache_set.md)
  : Write cached RDS value

## Reporting & Shiny UI Helpers

Export reports, build interactive Plotly Pareto plots, and manage Shiny
UI state.

- [`write_optimization_report()`](https://pascalcrepey.github.io/revlitquery/reference/write_optimization_report.md)
  : Write all optimization artifacts (CSVs, Markdown report, PNG plot)
  to \`output_dir\`
- [`build_pareto_plotly()`](https://pascalcrepey.github.io/revlitquery/reference/build_pareto_plotly.md)
  : Build interactive Plotly Pareto frontier chart
- [`build_coverage_display_df()`](https://pascalcrepey.github.io/revlitquery/reference/build_coverage_display_df.md)
  : Build display table for the Seed Coverage Matrix tab
- [`recompute_recommendation()`](https://pascalcrepey.github.io/revlitquery/reference/recompute_recommendation.md)
  : Recompute recommended query on an existing Pareto frontier when
  sliders move
- [`build_config_from_inputs()`](https://pascalcrepey.github.io/revlitquery/reference/build_config_from_inputs.md)
  : Validate and assemble optimizer config list from Shiny sidebar
  inputs
- [`parse_term_list()`](https://pascalcrepey.github.io/revlitquery/reference/parse_term_list.md)
  : Split a textarea string of terms on commas, semicolons, or newlines
- [`parse_custom_blocks_yaml()`](https://pascalcrepey.github.io/revlitquery/reference/parse_custom_blocks_yaml.md)
  : Parse optional custom concept blocks from YAML text
- [`blocks_to_yaml()`](https://pascalcrepey.github.io/revlitquery/reference/blocks_to_yaml.md)
  : Convert discovered concept blocks list into editable YAML text

## Internal Helpers

- [`add_custom_query_to_result()`](https://pascalcrepey.github.io/revlitquery/reference/add_custom_query_to_result.md)
  : Merge an evaluated custom query into an existing or minimal
  optimization result
- [`assign_exclusion_concept()`](https://pascalcrepey.github.io/revlitquery/reference/assign_exclusion_concept.md)
  : Classify a candidate exclusion term into a scientific exclusion
  concept
- [`assign_term_families()`](https://pascalcrepey.github.io/revlitquery/reference/assign_term_families.md)
  : Vectorized wrapper for \`assign_term_family()\`
- [`assign_term_family()`](https://pascalcrepey.github.io/revlitquery/reference/assign_term_family.md)
  : Classify a candidate search term into a coherent semantic PICO
  family
- [`blocks_to_yaml()`](https://pascalcrepey.github.io/revlitquery/reference/blocks_to_yaml.md)
  : Convert discovered concept blocks list into editable YAML text
- [`build_concept_blocks()`](https://pascalcrepey.github.io/revlitquery/reference/build_concept_blocks.md)
  : Build stem-disjoint, conditionally selective concept blocks from
  scored candidate terms
- [`build_config_from_inputs()`](https://pascalcrepey.github.io/revlitquery/reference/build_config_from_inputs.md)
  : Validate and assemble optimizer config list from Shiny sidebar
  inputs
- [`build_coverage_display_df()`](https://pascalcrepey.github.io/revlitquery/reference/build_coverage_display_df.md)
  : Build display table for the Seed Coverage Matrix tab
- [`build_pareto_plotly()`](https://pascalcrepey.github.io/revlitquery/reference/build_pareto_plotly.md)
  : Build interactive Plotly Pareto frontier chart
- [`cache_get()`](https://pascalcrepey.github.io/revlitquery/reference/cache_get.md)
  : Read cached RDS value if available
- [`cache_key()`](https://pascalcrepey.github.io/revlitquery/reference/cache_key.md)
  : Build deterministic cache key
- [`cache_set()`](https://pascalcrepey.github.io/revlitquery/reference/cache_set.md)
  : Write cached RDS value
- [`clean_latex_text()`](https://pascalcrepey.github.io/revlitquery/reference/clean_latex_text.md)
  : Clean LaTeX markup and normalize whitespace in BibTeX field text
- [`clean_title_for_query()`](https://pascalcrepey.github.io/revlitquery/reference/clean_title_for_query.md)
  : Clean title string for search queries (strip quotes/brackets/colons
  that break E-utilities)
- [`compute_pareto_frontier()`](https://pascalcrepey.github.io/revlitquery/reference/compute_pareto_frontier.md)
  : Compute non-dominated Pareto frontier indicator
- [`discover_exclusion_hints()`](https://pascalcrepey.github.io/revlitquery/reference/discover_exclusion_hints.md)
  : Discover concept-grouped candidate exclusion (NOT / -) terms from
  non-seed PubMed results
- [`.is_code_family()`](https://pascalcrepey.github.io/revlitquery/reference/dot-is_code_family.md)
  : Check whether two terms share an alphanumeric subtype/code prefix
  (e.g. h5n1 and h5)
- [`.term_subsumes()`](https://pascalcrepey.github.io/revlitquery/reference/dot-term_subsumes.md)
  : Check whether term \`t_broad\` lexically subsumes \`t_narrow\`
  inside an OR block
- [`enrich_seeds()`](https://pascalcrepey.github.io/revlitquery/reference/enrich_seeds.md)
  : Enrich BibTeX seed papers via PubMed
- [`escape_regex()`](https://pascalcrepey.github.io/revlitquery/reference/escape_regex.md)
  : Escape regex metacharacters in literal term
- [`estimate_local_query()`](https://pascalcrepey.github.io/revlitquery/reference/estimate_local_query.md)
  : Estimate local seed coverage and rough PubMed hit count for a
  candidate query
- [`evaluate_custom_query()`](https://pascalcrepey.github.io/revlitquery/reference/evaluate_custom_query.md)
  : Evaluate a single user-supplied PubMed query against enriched seeds
- [`extract_candidate_terms()`](https://pascalcrepey.github.io/revlitquery/reference/extract_candidate_terms.md)
  : Extract and score candidate search terms from enriched seed papers
- [`extract_pmid_value()`](https://pascalcrepey.github.io/revlitquery/reference/extract_pmid_value.md)
  : Extract a PubMed ID (PMID) from explicit fields or URL/note text
- [`generate_initial_candidates()`](https://pascalcrepey.github.io/revlitquery/reference/generate_initial_candidates.md)
  : Generate diverse initial Stage 1 query candidates
- [`get_stopwords()`](https://pascalcrepey.github.io/revlitquery/reference/get_stopwords.md)
  : Stopword list (standard English + academic/methods filler)
- [`google_scholar_url()`](https://pascalcrepey.github.io/revlitquery/reference/google_scholar_url.md)
  : Clickable Google Scholar search URL
- [`gs_char_count()`](https://pascalcrepey.github.io/revlitquery/reference/gs_char_count.md)
  : Count characters in translated Google Scholar query
- [`is_valid_gs_query()`](https://pascalcrepey.github.io/revlitquery/reference/is_valid_gs_query.md)
  : Check Google Scholar \<= 256 character constraint
- [`make_query()`](https://pascalcrepey.github.io/revlitquery/reference/make_query.md)
  : Construct a canonical query object
- [`normalize_doi()`](https://pascalcrepey.github.io/revlitquery/reference/normalize_doi.md)
  : Normalize a DOI string to lowercase bare identifier
- [`normalize_search_text()`](https://pascalcrepey.github.io/revlitquery/reference/normalize_search_text.md)
  : Normalize text for n-gram extraction (preserves alphanumeric tokens
  like h5n1)
- [`optimize_queries()`](https://pascalcrepey.github.io/revlitquery/reference/optimize_queries.md)
  : Run Stage 1 + Stage 2 Pareto query optimization
- [`optimize_search_queries()`](https://pascalcrepey.github.io/revlitquery/reference/optimize_search_queries.md)
  : Optimize semantically equivalent PubMed and Google Scholar queries
  for a BibTeX file
- [`parse_bibtex()`](https://pascalcrepey.github.io/revlitquery/reference/parse_bibtex.md)
  : Parse a BibTeX file into a normalized tibble of seed papers
- [`parse_custom_blocks_yaml()`](https://pascalcrepey.github.io/revlitquery/reference/parse_custom_blocks_yaml.md)
  : Parse optional custom concept blocks from YAML text
- [`parse_term_list()`](https://pascalcrepey.github.io/revlitquery/reference/parse_term_list.md)
  : Split a textarea string of terms on commas, semicolons, or newlines
- [`pubmed_efetch()`](https://pascalcrepey.github.io/revlitquery/reference/pubmed_efetch.md)
  : Fetch article metadata (title, abstract, MeSH, keywords) from PubMed
  efetch XML
- [`pubmed_esearch()`](https://pascalcrepey.github.io/revlitquery/reference/pubmed_esearch.md)
  : Query PubMed E-utilities esearch via HTTP POST
- [`pubmed_eval_query()`](https://pascalcrepey.github.io/revlitquery/reference/pubmed_eval_query.md)
  : Evaluate a PubMed query: total hit count + which seed PMIDs are
  retrieved
- [`pubmed_to_google_scholar_approx()`](https://pascalcrepey.github.io/revlitquery/reference/pubmed_to_google_scholar_approx.md)
  : Convert a raw PubMed query string into an approximate Google Scholar
  query
- [`pubmed_url()`](https://pascalcrepey.github.io/revlitquery/reference/pubmed_url.md)
  : Clickable PubMed search URL
- [`query_signature()`](https://pascalcrepey.github.io/revlitquery/reference/query_signature.md)
  : Deterministic canonical signature for deduplicating queries
- [`recompute_recommendation()`](https://pascalcrepey.github.io/revlitquery/reference/recompute_recommendation.md)
  : Recompute recommended query on an existing Pareto frontier when
  sliders move
- [`revlitquery_app()`](https://pascalcrepey.github.io/revlitquery/reference/revlitquery_app.md)
  : Create the revlitquery Shiny application object
- [`run_revlitquery_app()`](https://pascalcrepey.github.io/revlitquery/reference/run_revlitquery_app.md)
  : Launch the interactive revlitquery Shiny application
- [`shares_stem()`](https://pascalcrepey.github.io/revlitquery/reference/shares_stem.md)
  : Check whether two terms share any significant token stem
- [`split_keywords()`](https://pascalcrepey.github.io/revlitquery/reference/split_keywords.md)
  : Split BibTeX keywords string into cleaned lowercase character vector
- [`term_matches_text()`](https://pascalcrepey.github.io/revlitquery/reference/term_matches_text.md)
  : Check whether a candidate term matches each text in \`texts\`
  (case-insensitive, hyphen-flexible)
- [`term_stems()`](https://pascalcrepey.github.io/revlitquery/reference/term_stems.md)
  : Simple token stemmer for lexical overlap detection across concept
  blocks
- [`title_jaccard()`](https://pascalcrepey.github.io/revlitquery/reference/title_jaccard.md)
  : Compute alphanumeric lowercase token Jaccard similarity between two
  titles
- [`title_sig_words()`](https://pascalcrepey.github.io/revlitquery/reference/title_sig_words.md)
  : First N significant words from a title for fallback matching
- [`to_google_scholar_query()`](https://pascalcrepey.github.io/revlitquery/reference/to_google_scholar_query.md)
  : Translate canonical query to Google Scholar query string
- [`to_pubmed_query()`](https://pascalcrepey.github.io/revlitquery/reference/to_pubmed_query.md)
  : Translate canonical query to PubMed E-utilities query string
- [`write_optimization_report()`](https://pascalcrepey.github.io/revlitquery/reference/write_optimization_report.md)
  : Write all optimization artifacts (CSVs, Markdown report, PNG plot)
  to \`output_dir\`

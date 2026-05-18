# brainstorm: multi-query genome comparison

## Goal

Extend GeneScreen from one reference genome plus one query genome into one reference genome plus multiple query genomes, while preserving the existing Gene ID, Location, and Sequence input workflows. The report should compare the ref-derived input sequence against all selected query genomes, then support multi-track overview and detail visualizations across the matched regions.

The implementation should treat the current two-genome comparison as the `1 ref + 1 query` case of the same multi-query architecture, not as a separate legacy-only path.

## What I Already Know

* Current branch for this work: `multi-query-genomes`.
* Current task directory: `.trellis/tasks/05-18-multi-query-genome-comparison`.
* Reference genome is required and must have annotation for Gene ID mode.
* Multiple query genomes are selected; query annotations are optional.
* If a query genome has annotation, its annotation should render on that genome's corresponding LINKVIEW/detail track. If not, the track still renders without annotation.
* Input sequence can still come from the existing three sources:
  * Gene ID: extract ref gene from reference genome annotation.
  * Location: extract ref region from reference genome.
  * Sequence: use user-provided sequence as the ref-derived sequence/ref track for comparison.
* User's desired workflow example:
  * Ref `A`, queries `B`, `C`, `D`.
  * First align `A_gene` against `B/C/D`, producing candidate hits such as `B1/B2/B3`, `C1/C2`, `D1`.
  * Then append pairwise comparisons among candidate query regions, e.g. `B_gene vs C_gene`, `B vs D`, `C vs D`.
* Existing parameters should remain available:
  * upstream/downstream expansion
  * identity threshold
  * minimum alignment length
  * merge gap
* Overview chart needs a new multi-genome layout:
  * Tracks grouped by genome/chromosome.
  * Each genome's best match (`X1`) should be visually distinguished, likely green.
  * Overview selections remain interactive.
  * Selecting a candidate in the overview updates the detail view's current candidate for that genome.
* Detail view default should show the best hit from each genome, e.g. `A1-B1-C1-D1`.
* Detail view should support:
  * multi-track rendering
  * per-track dropdown/collapse controls to switch candidates, e.g. `A1-B2-C1-D1`
  * arbitrary track order via drag/reorder, e.g. `A1-C1-D1-B1`
  * switching among precomputed candidates/links inside the static report without recomputing alignments
  * showing only the currently selected candidate per genome track and the links/variants among those selected candidates
  * keeping non-selected candidates available in dropdown controls instead of rendering them as background links in the detail chart
  * preserving real genomic coordinates and strand/orientation in data while mapping them to the current display direction at render time
  * visually indicating reverse-strand candidates/links instead of rewriting their source coordinates
* Report's first two statistics cards should support `1 ref + N query` context.
  * The first card summarizes input/candidate context, including `1 ref + N query`, total candidate count, and per-query candidate counts.
  * The second card summarizes alignment/variant context for the current detail selection, including total links and total SNP/INS/DEL counts, with concise per-genome-pair breakdowns.
* GUI should support multiple query genome inputs in Gene ID, Location, and Sequence pages.
* GUI multi-query selection should use a shared control pattern:
  * Reference genome remains a single selection.
  * Query genomes are managed as an "added queries" list.
  * Users can add query genomes from a dropdown/search control.
  * Each selected query row shows genome name, annotation availability/status, and a remove action.
  * The same control is reused across Gene ID, Location, and Sequence pages.
* Output should remain organized by input item:
  * Multi Gene ID / Location / Sequence workflows create one subdirectory per input item.
  * Each input item directory contains the ref-derived sequence/annotation artifacts, per-query candidate artifacts, query-query pairwise artifacts, and report files.
  * Suggested structure:
    * `output/<input-id>/ref/`
    * `output/<input-id>/queries/<query-name>/`
      * `ref__<query-name>.coords`
      * `ref__<query-name>.snps`
      * `ref__<query-name>.hl`
    * `output/<input-id>/pairwise/<query-a>__<query-b>/`
      * `<query-a>__<query-b>.coords`
      * `<query-a>__<query-b>.snps`
      * `<query-a>__<query-b>.hl`
    * `output/<input-id>/report/`
      * `index.html`
      * `data.json`

## Current Code Constraints

* Analysis is currently centered on a single `BlastAligner.align(reference, query, prefix)` call.
  * Gene ID and Location align a ref-derived FASTA against one query genome FASTA.
  * Sequence mode aligns user FASTA against one reference genome FASTA.
* GUI analysis processors live in `GeneScreen_GUI/core/analysis.py`:
  * `GeneIDProcessor`
  * `LocationProcessor`
  * `SequenceProcessor`
* CLI analysis processors live in `GeneScreen_CLI/bin/GeneScreen.py` with similar names.
* Existing visualization/report flow is dual-track and centered on `AlignmentGroup` and `VisualizationResult` in:
  * `GeneScreen_GUI/core/visualizer.py`
  * `GeneScreen_CLI/bin/GeneScreenVisualizer.py`
* Current GUI selector `GenomePairSelector` only supports one ref and one query.
* Query GFF extraction currently assumes a single query GFF and one lower track.
* Recent behavior: multi-input output paths were changed to `out/<safe_id>/` and should remain compatible.
* Current single-query/two-genome code should be updated to use the shared multi-query data model and output/report conventions wherever practical.

## Assumptions

* MVP should integrate into the existing Gene ID, Location, and Sequence pages rather than creating a separate top-level mode.
* CLI should also support multiple query genomes for parity, likely by allowing `-qry` to accept multiple values or a file.
* Existing single-query behavior should remain compatible.
* Two-genome comparison is considered `N=1` query under the new architecture.
* Pairwise query-query comparisons should be derived from matched query regions rather than whole-genome comparisons.
* Sequence mode treats the user-provided sequence as the ref-derived sequence/ref track. The selected reference genome is context for naming/index/optional annotation, not a mandatory extraction source.

## Open Questions

* None currently. No additional must-have report interactions beyond candidate switching, overview selection, and track drag/reorder for the first implementation.

## Requirements

* Support selecting multiple query genomes in GUI for Gene ID, Location, and Sequence workflows.
* GUI uses a shared ref selector plus add/remove query list for Gene ID, Location, and Sequence workflows.
* Support multiple query genomes in CLI while preserving current single-query usage where practical.
* New CLI multi-query syntax should model each genome as a complete entry:
  * Preloaded reference entry: `-ref <ref-name>`
  * Explicit reference entry: `-ref <ref-name> <ref.fa> <ref.gff|ref.gff3>`
  * Preloaded query entry, repeated: `-qry <query-name>`
  * Explicit query entry, repeated: `-qry <query-name> <query.fa> [query.gff|query.gff3]`
  * Example: `-ref A A.fa A.gff3 -qry B B.fa B.gff3 -qry C C.fa`
  * Example with preloaded genomes: `-ref A -qry B -qry C`
  * Mixed entries are valid because each genome entry is resolved independently, e.g. `-ref A -qry B B.fa B.gff3`.
  * Query annotation remains optional per query.
  * Reference annotation is required for ref-derived Gene ID extraction, whether it comes from the preloaded genome database or explicit CLI path.
  * Single-query legacy syntax such as `-ref Nippon -qry ZS97` remains valid through the preloaded-entry form.
* Extract/prepare the ref-derived input sequence once per input item.
* In Sequence mode, the user-provided sequence is the ref-derived input sequence even if it was not extracted from the selected reference genome.
* Align the ref-derived sequence against every query genome.
* Preserve per-query optional annotation rendering.
* Generate a unified multi-genome report per input item.
* Preserve the current multi-input organization by creating one output subdirectory per input item.
* Multi-query artifacts inside each input item directory are grouped by ref, per-query candidates, pairwise query-query comparisons, and report files.
* Existing `.coords/.snps/.hl` compatibility artifacts are still generated per comparison pair.
  * Ref-query artifacts live under `queries/<query-name>/`.
  * Query-query artifacts live under `pairwise/<query-a>__<query-b>/`.
* Static report data is written as `report/data.json`; `report/index.html` loads and renders that data.
* Generate overview and detail visualizations that can represent multiple genomes and multiple candidate hits per genome.
* Detail view displays the currently selected candidate combination only; non-selected candidates are switchable but not rendered as background context in the detail chart.
* Overview candidate selection updates the detail view candidate selection for the clicked genome while preserving the current selections and order of other tracks.
* Candidate and link data keeps true genomic coordinates and strand/orientation.
* Rendering maps candidate/link coordinates into the current track display direction and marks reverse-strand candidates/links clearly.
* Report statistics combine overall context with per-query/per-pair summaries.
* Variant/link statistics in the detail card are based on the currently selected candidate combination.
* Preserve existing single-query report behavior or provide a compatible path for it.
* Current non-multi/two-genome code paths should be aligned to the new shared conventions so that `1 ref + 1 query` and `1 ref + N queries` produce structurally consistent artifacts and reports.
* Pairwise query-query comparisons should expose a user-configurable candidate limit.
  * Default value is `3`.
  * CLI omitted parameter means top `3`.
  * CLI full/unlimited mode is explicit, e.g. `--pairwise-all`.
  * GUI default value is `3`; user can choose full/unlimited mode.
  * Candidate limit applies per query genome before generating pairwise query-query comparisons.
* Pairwise query-query comparisons are precomputed during analysis for the kept candidate set.
* Generated HTML reports are static: candidate switching and track reordering only change which precomputed candidates/links are displayed, and do not rerun BLAST or call the backend.

## Implementation Plan

* Phase 1: migrate the current two-genome comparison into the shared `1 ref + N query` framework with `N=1`.
  * Define shared genome-entry, candidate, pairwise-link, statistics, and report data schema.
  * Keep current Gene ID, Location, and Sequence single-query behavior working through the new schema.
  * Generate the new per-input output structure and `report/index.html + report/data.json` for the `N=1` case.
  * Preserve existing `.coords/.snps/.hl` artifacts under the new pair-oriented directories.
* Phase 2: extend CLI and GUI inputs to accept multiple query genome entries.
  * CLI supports repeated `-qry` entries with preloaded, explicit, and mixed forms.
  * GUI replaces the single query selector with the shared add/remove query list.
* Phase 3: add multi-query candidate collection and query-query pairwise precomputation.
  * Align the ref-derived sequence against every query genome.
  * Rank and keep candidates per query genome.
  * Generate pairwise query-query comparisons among kept candidates.
* Phase 4: implement multi-track report rendering and interactions.
  * Overview supports multi-genome candidate display and selection.
  * Detail supports selected candidate combination, per-track candidate switching, drag/reorder, optional annotations, and reverse-strand rendering.
  * Statistics cards update for `1 ref + N query` and selected detail combinations.

## Acceptance Criteria

* [ ] Gene ID mode accepts one ref and two or more query genomes in GUI and CLI.
* [ ] Location mode accepts one ref and two or more query genomes in GUI and CLI.
* [ ] Sequence mode accepts one selected ref context and two or more comparison genomes in GUI and CLI, with behavior defined before implementation.
* [ ] Multi-input runs keep one subdirectory per input item.
* [ ] Each input item directory clearly separates ref-derived artifacts, per-query candidate artifacts, pairwise artifacts, and report files.
* [ ] Existing `.coords/.snps/.hl` style artifacts are available per ref-query and query-query pair.
* [ ] Each report directory contains `index.html` and `data.json`.
* [ ] Query genomes without annotation still produce tracks and reports.
* [ ] Query genomes with annotation show annotation on their corresponding tracks.
* [ ] Report overview highlights each genome's top match.
* [ ] Clicking a candidate in the overview switches that genome's candidate in the detail view.
* [ ] Report statistics cards show `1 ref + N query` context, per-query candidate counts, and selected-combination variant/link counts.
* [ ] Report detail view defaults to the best candidate per genome.
* [ ] Detail view allows switching candidates per genome track.
* [ ] Detail view allows track order changes.
* [ ] Detail view only renders links/variants for the currently selected candidate combination.
* [ ] Reverse-strand candidates and links render correctly without mutating stored genomic coordinates.
* [ ] Switching candidates in an HTML report does not require rerunning BLAST or making a backend call.
* [ ] Users can control the number of candidate hits per genome used for pairwise query-query comparisons, including an unlimited mode.
* [ ] Existing single-query workflows still pass existing checks.
* [ ] Existing two-genome runs use the same result schema/report conventions as the multi-query implementation with one query.

## Definition of Done

* Tests added or updated for core data model / candidate selection logic.
* GUI and CLI paths produce consistent result data.
* Python compile checks pass.
* Manual or automated smoke test covers single-query and multi-query cases.
* Documentation/README updated for CLI and GUI usage.
* Work committed and pushed after completion.

## Out of Scope

* Replacing BLAST with MUMmer/minimap2.
* Reworking genome database storage beyond what multi-query selection requires.
* Full all-vs-all whole-genome comparison.

## Technical Notes

* Existing two-track visualizer likely needs a new multi-track result model rather than incremental patching of single-query `VisualizationResult`.
* Need a stable schema for:
  * genome entries
  * per-genome candidate hits
  * pairwise links between selected candidates
  * optional per-track GFF features
* Pairwise query-query candidate strategy:
  * A ref-derived sequence is first aligned to each query genome.
  * Hits per query genome are ranked by the same score used for overview/default selection.
  * Candidate ranking is by total alignment length descending, then identity descending, then coverage descending.
  * The configured candidate limit keeps the top N per query genome.
  * Unlimited mode keeps all passing hits.
  * CLI default is top `3`; full mode must be explicitly requested.
  * Pairwise links are generated among the kept candidates at report generation time.
  * Report interactions switch among precomputed candidate/link combinations only.

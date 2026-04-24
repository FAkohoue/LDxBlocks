## LDxBlocks 0.3.1.9000 (development)

### New feature: block_match = "position" in compare_block_effects() and compare_gwas_effects()

Both cross-population comparison functions now support matching LD blocks between
populations by **genomic interval overlap** rather than by `block_id` string equality.

**Why this matters.** LD block boundaries are population-specific: the same causal
QTL region may be carved into a 100 kb block in Population A and a 130 kb block in
Population B, producing different `block_id` strings (`"block_1_10000_85000"` vs
`"block_1_10000_91000"`). With the default `block_match = "id"`, these would not be
compared — the block appears as Pop1-only and the replication signal is lost.

**New parameters:**

- `block_match = c("id", "position")` — `"id"` (default) preserves backward
  compatibility; `"position"` matches by Intersection-over-Union (IoU) in base pairs.
- `overlap_min = 0.50` — minimum IoU for two blocks to be considered the same region.
  Blocks below this threshold are labelled `"pop1_only"`.

**New output column `match_type`** in `$concordance`:

- `"exact"` — same `block_id` string (boundaries identical)
- `"position"` — matched by genomic overlap (boundaries differ but IoU ≥ `overlap_min`)
- `"pop1_only"` — no Pop2 block overlaps this Pop1 block at the threshold
- `NA` — no block tables were supplied

**New internal function `.match_blocks_by_position()`** performs the interval join
using a CHR-filtered IoU search. For each Pop1 block, it finds the best-matching Pop2
block by IoU and records the match type. Accessible via `LDxBlocks:::.match_blocks_by_position()`.

8 new tests in `test-association.R` (122 → 130 total).

## LDxBlocks 0.3.1.9000 (development)

### New function: compare_gwas_effects()

Cross-population effect concordance from **external GWAS results**. Complements
`compare_block_effects()` for users who ran association analysis outside
LDxBlocks (GAPIT, TASSEL, FarmCPU, PLINK, or any other tool).

**Two input paths:**

- **Pre-mapped (recommended):** supply `define_qtl_regions()` output for each
  population. Block assignment is explicit and auditable.
- **Raw GWAS + blocks (convenience):** supply raw GWAS data frames and block
  tables; `define_qtl_regions()` is called internally.

**SE derivation.** When the SE column is absent (common in GAPIT/FarmCPU
output), SE is derived from the z-score: `SE = |BETA| / |Φ⁻¹(P/2)|`. The
`se_derived_pop1` / `se_derived_pop2` output columns flag which population
required this step.

**Key differences from `compare_block_effects()`:**
External GWAS produces one lead SNP per block rather than multiple haplotype
allele effects. Consequently:
- `effect_correlation` is always `NA` (needs ≥ 3 alleles).
- `direction_agreement` is 0 or 1 only.
- `Q_stat`, `Q_p`, `I2` are always `NA` (Cochran Q undefined with df = 0).
- `replicated` uses `meta_p ≤ 0.05` instead of `Q_p > 0.05`.

**GWAS-specific output columns** added to `$concordance`:
`lead_snp_pop1`, `lead_snp_pop2`, `lead_p_pop1`, `lead_p_pop2`,
`se_derived_pop1`, `se_derived_pop2`, `both_pleiotropic`.

**Flexible column naming.** The `beta_col`, `se_col`, and `p_col` arguments
accept any column names (e.g. `"effect"`, `"std_err"`, `"pvalue"`). The
`Marker` column is accepted as an alias for `SNP`.

**Output class.** Returns `LDxBlocks_effect_concordance` — the same class as
`compare_block_effects()`. The existing `print()` method works immediately.

18 new tests in `test-association.R` (104 → 122 total).

## LDxBlocks 0.3.1.9000 (development)

### New function: compare_block_effects()

Cross-population haplotype effect concordance. Takes two `test_block_haplotypes()`
result objects and returns per-block statistics for systematic GWAS replication:

- **IVW meta-analysis**: inverse-variance weighted combined effect and SE per
  block, identical framework to two-sample Mendelian randomisation.
- **Cochran Q heterogeneity**: tests whether effect sizes differ significantly
  between populations. Significant Q (low Q_p) flags G×E interaction or
  population-specific LD structure differences.
- **I² inconsistency**: 0–100% measure of between-population heterogeneity.
  Values > 50% indicate that the two populations are telling a different
  biological story at that block.
- **Direction agreement**: fraction of shared alleles with the same effect sign.
  Controlled by `direction_threshold` (default 0.75).
- **`replicated` flag**: composite criterion — `enough_shared AND
  directionally_concordant AND Q_p > 0.05`.
- **Block boundary diagnostics**: when `blocks_pop1` and `blocks_pop2` are
  supplied, `boundary_overlap_ratio` is **automatically computed** (not
  user-set) as bp(intersection) / bp(union) for every block. This quantifies
  how similarly the two populations carved the region into LD blocks.
  `boundary_overlap_warn` (input parameter, default `0.80`) is the threshold
  below which `boundary_warning = TRUE` is set in the output.
- **`$shared_alleles`**: per-allele detail table with `ivw_effect`, `ivw_SE`,
  `direction_agree`, and raw effects/SEs from both populations.
- 15 new tests in `test-association.R` (89 → 104 total).

### Enhancement: simpleM multiple-testing correction in test_block_haplotypes()

`test_block_haplotypes()` now implements the simpleM procedure
(Gao et al. 2008, 2010, 2011) for LD-aware multiple-testing correction of
correlated haplotype allele tests. simpleM estimates the effective number of
independent tests (Meff) from the eigenspectrum of the haplotype allele dosage
correlation matrix, replacing raw test counting with an LD-aware effective
count that is less conservative than Bonferroni while providing family-wise
error control.

**New parameters:**

- `sig_metric`: which p-value drives `significant` / `significant_omnibus`.
  One of `"p_wald"`, `"p_fdr"`, `"p_simplem"` (Bonferroni-style), or
  `"p_simplem_sidak"` (Šidák-style, recommended). Default `"p_wald"`.
- `meff_scope`: scope for Meff estimation — `"chromosome"` (recommended),
  `"global"`, or `"block"`. Default `"chromosome"`.
- `meff_percent_cut`: variance threshold for simpleM eigendecomposition.
  Default `0.995` (99.5%), following the original simpleM recommendation.
- `meff_max_cols`: chunk size for large eigendecompositions. Default `1000L`.

**New output columns — always present regardless of `sig_metric`:**

- `allele_tests`: `Meff`, `alpha_simplem`, `alpha_simplem_sidak`,
  `p_simplem`, `p_simplem_sidak`, `p_fdr`
- `block_tests`: `Meff`, `alpha_simplem`, `alpha_simplem_sidak`,
  `p_omnibus_fdr`, `p_omnibus_simplem`, `p_omnibus_simplem_sidak`
  (plus `p_omnibus_adj` retained for backward compatibility)

**New return list elements:** `meff_scope`, `meff_percent_cut`, `meff`
(nested list of Meff summaries per trait: `$allele$global`,
`$allele$chromosome`, `$allele$block`, `$block$global`, `$block$chromosome`).

13 new tests in `test-association.R` (76 → 89).

### Bug fix: read_phased_vcf() dot-ID synthesis

When a phased VCF has `.` or empty string in the ID column (column 3),
`read_phased_vcf()` now synthesises `CHR_POS` identifiers for those rows,
matching the behaviour of `read_geno()`. Without this fix, all-dot VCFs
produced duplicate `rownames(hap1)` that caused `extract_haplotypes()` to
fail when slicing by SNP name.

- Handles `NA`, `"."`, and `""` (empty string) row-by-row.
- Real rsIDs in mixed VCFs are preserved; only missing IDs are replaced.
- Verbose message reports the count of synthesised IDs.
- 3 new regression tests in `test-phasing.R` (30 → 33).

Note: this gap only affects direct calls to `read_phased_vcf()` on
externally phased VCFs with dot IDs. The `run_ldx_pipeline()` path was
already safe because `read_geno()` synthesises `CHR_POS` IDs before writing
the cleaned VCF to Beagle.

## LDxBlocks 0.3.1.9000 (development patch)

### Breaking change: SNP ID separator changed from `:` to `_`

When a VCF or GDS file has no rsID (i.e. the ID field is `.` or empty),
`read_geno()` previously generated fallback SNP identifiers in `CHR:POS`
format (e.g. `"1:4106"`). These are now generated as `CHR_POS` format
(e.g. `"1_4106"`).

**Reason:** The colon `:` is a reserved character in many genomic file
formats (BED, VCF INFO field, R data frames used as rownames) and causes
silent parsing errors downstream. The underscore `_` is safe in all contexts.

**Impact:** Any cached bigmemory backing files (`ldxblocks_bm*.rds`,
`ldxblocks_bm*.bin`, `ldxblocks_bm*.desc`) built with the old format must
be deleted before rerunning. The pipeline will rebuild them automatically.

**Files to delete before rerunning:**
```r
results_dir <- "/your/results/dir"
for (f in list.files(results_dir, pattern = "^ldxblocks_bm", full.names = TRUE))
  file.remove(f)
```

## LDxBlocks 0.3.1.9000 (development patch)

### Bug fixes

- **`extract_chr_haplotypes_cpp()` — `retained_idx` field added (critical)**
  The C++ extractor now returns a `retained_idx` integer vector (1-based) giving
  the row index in the input `block_sb`/`block_eb` arrays for each retained block.
  Previously the compact `hap_strings` and `n_snps` arrays (indexed over retained
  blocks only) were paired in R with `chr_blk_srt[b, ]` (indexed over all blocks),
  causing a systematic mismatch: singleton coordinates were paired with multi-SNP
  haplotype strings and inflated `n_snps` values. Both the backend streaming path
  and the matrix path in `extract_haplotypes()` now use `ret_idx <- cpp_res$retained_idx`
  to fetch block metadata, eliminating the mismatch entirely.
  Sanity checks stop immediately with an informative message if `retained_idx` is
  missing or malformed (indicating the old compiled binary is still in use).
  Both extraction paths also pre-filter `chr_blk_srt` with `findInterval` before
  calling C++, providing a first-pass guard independent of compilation.

- **Pipeline QC singleton check** — `run_ldx_pipeline()` now asserts that no
  `block_info` row has `start_bp == end_bp` with `n_snps > 1` immediately after
  `extract_haplotypes()`, stopping with a clear message if the stale binary is
  still loaded.

- **`.validate_hap_output()` added** — called at end of pipeline; reports NA
  count in feature matrix, duplicated `hap_id` values, n_snps range, and
  hap_id / column alignment.

### Improvements

- **`alleles` column in writer** — `write_haplotype_numeric()` now decodes
  nucleotide sequences directly from `hap_info$hap_string` + `snp_info` per row,
  without calling `decode_haplotype_strings()`. Single source of truth: the
  dosage string stored in `hap_info` is always the one used to build the matrix
  column.

- **Backend extraction list accumulation** — the backend streaming path in
  `extract_haplotypes()` now accumulates `block_info` rows in a pre-allocated
  list (`bi_rows`) and binds once at the end with `data.table::rbindlist()`,
  matching the matrix path and eliminating O(n²) rbind growth.

- **`.prefilter_blocks_by_span()` helper** — shared `findInterval`-based
  pre-filter extracted to a single internal function used by both extraction
  paths; eliminates duplicated logic and future drift risk.

- **`extract_chr_haplotypes_phased_cpp()` — new C++ phased extractor (Item 5)**
  A dedicated phased haplotype extractor is now compiled into `ld_core.cpp`.
  It accepts `hap1_chr` and `hap2_chr` (0/1 gamete matrices) and builds
  `"gamete1|gamete2"` strings in one OpenMP pass, counting gamete frequencies
  correctly (each individual contributes two observations). Returns the same
  contract as `extract_chr_haplotypes_cpp()`: `retained_idx`, `retained_start_bp`,
  `retained_end_bp`, `hap_freq`, `freq_dominant`, etc. The matrix path in
  `extract_haplotypes()` now dispatches to this function when `isp = TRUE`,
  eliminating the previous R loop that called `build_hap_strings_cpp()` twice
  per block and concatenated results with `paste0()`.
  The R `for(b)` loop is now identical for phased and unphased: both read
  `cpp_res$hap_strings[[b]]` from their respective C++ extractor.

- **Imputed backend parameter fingerprinting (Item 6)** — `.make_imputed_backend()`
  now accepts `maf_cut`, `min_callrate`, and `impute_method` parameters and
  saves them as `ldxblocks_bm_imputed_params.rds` alongside the other four
  backing files. On reattach the saved fingerprint is compared to the current
  parameters; if they differ (or the fingerprint file is absent, indicating an
  old-format cache), all five files are cleaned and the backend is rebuilt.
  Without this, a run with `maf_cut = 0.10` would silently reuse an imputed
  backend built with `maf_cut = 0.05`. `run_ldx_pipeline()` passes `maf_cut`,
  `min_callrate`, and the `impute` argument through to `.make_imputed_backend()`.

- **C++ as bp coordinate truth — `retained_start_bp` and `retained_end_bp` (Item 8)**
  Both `extract_chr_haplotypes_cpp()` and `extract_chr_haplotypes_phased_cpp()`
  now return `retained_start_bp` and `retained_end_bp` integer vectors: for each
  retained block, the values are direct copies of `block_sb[b]` and `block_eb[b]`
  from the C++ arrays — no R-side index lookup into `chr_blk_srt`. Both R
  extraction loops check `!is.null(cpp_res$retained_start_bp)` and read
  coordinates directly from C++ when available. The R-side fallback
  (`chr_blk_srt[ret_idx[b], ]`) remains for binary-compatibility with old
  compiled objects. This removes the last avenue for a coordinate mismatch
  between the C++ compact array index and the R full-block-table row index.

---

## LDxBlocks 0.3.1 (development)

### Performance: C++ overlap resolution and WGS stall fixes

The `Big_LD()` post-segment pipeline now routes through C++ for all overlap
resolution, eliminating several O(n²)–O(n·p) bottlenecks that caused
chromosome processing to stall after segment detection completed.

**`resolve_overlap_cpp()` — new exported C++ function**

Replaces the two sequential R `.resolve_overlap()` calls with a single C++
pass implementing the identical cumulative-score split rule. Four improvements
over the R version:

- **BLAS DGEMM scoring** — for each overlapping block pair, disputed SNP
  scores are computed as `rowMeans(C_L²) − rowMeans(C_R²)` where
  `C_L = (Z_overlap.t() × Z_left_reps) / (n−1)` via Armadillo matrix
  multiply. This replaces a `vapply()` loop calling `col_r2_cpp()` per SNP
  against all p columns (314k for chr1). Cost per disputed SNP drops from
  O(n × p) to O(n × 2k_rep) — a 15,700× reduction for chr1.
- **Lazy column cache** — each column of `adj_mat` is standardised at most
  once. Columns never accessed as representatives or overlap SNPs are never
  touched. For ~1,000 overlapping pairs using 20 reps each, fewer than 1%
  of chr1's 314k columns are ever standardised.
- **Single pass** — one call replaces two. The second R call was needed
  because the first pass could create new adjacent overlaps; the C++ while
  loop handles this naturally within the same pass.
- **OpenMP over pairs** — overlapping pairs are collected first in an O(n_blocks)
  scan, then resolved in `#pragma omp parallel for` (scores computed in
  parallel; boundaries applied serially in reverse to preserve index stability).

**R-level running counters**

Four additional O(n) scan operations replaced with O(1) running counters:

- `sum(!is.na(LDblocks[,1L]))` called 419 times per chromosome → `ld_count`
- `max(which(!is.na(newLDblocks[,1L])))` in re-merge loop → `remerge_count`
- `min(which(is.na(newLDblocks[,1L])))` in re-merge loop → `remerge_count`
- `done <- LDblocks[!is.na(LDblocks[,1L]),]` → `LDblocks[seq_len(ld_count),]`

**Re-merge skip when modeNum=2**

The optional re-merge loop across forced cut-points is now skipped when all
cut-points are forced (`length(atfcut) >= length(cutpoints) - 2`), which
occurs when `cutsequence.modi` switches to `modeNum=2`. In this mode every
cut is forced by definition, so re-running CLQD across all boundaries
produces no benefit. With 419 forced segments on chr1 this eliminated
potentially hundreds of CLQD calls in a serial R loop.

**`intersect()` → `findInterval()` in re-merge check**

The per-pair forced-cut check `length(intersect(eb[2L]:nb[1L], atfcut)) > 0L`
— O(n_blocks × n_atfcut) total — replaced with
`findInterval(eb[2L], atfcut)` binary search: O(log n_atfcut) per call.

**`score_overlap_cpp()` and `resolve_seam_cpp()` (internal static helpers)**

Two additional C++ helpers in `src/ld_core.cpp` implement the BLAS scoring
kernel and a seam-local resolver for future use at truly massive scale (pangenome
panels where the global O(n_blocks) scan itself becomes a bottleneck). These
are `static` functions not exported to R.

### Tests

- `test-resolve-overlap.R` (307 lines, new file): parallel property tests for
  both `.resolve_overlap()` (R reference) and `resolve_overlap_cpp()` (C++),
  plus cross-validation tests verifying identical boundaries on all six
  biological cases: non-overlapping passthrough, union merge (B inside A),
  all-left path, all-right path, tie-breaking (zero-variance overlap),
  mixed cumulative-score split. Output invariant tests: no new blocks created,
  output is non-overlapping, global span preserved.
- `test-cpp.R`: five additional tests for `resolve_overlap_cpp()` covering
  dimensionality, non-overlapping passthrough, no remaining overlaps,
  start ≤ end invariant, and R vs C++ agreement on random overlapping input.


### Haplotype association testing and breeding decision functions

Four new exported functions complete the statistical inference and breeding
decision layer.

**Haplotype association testing**

- **`test_block_haplotypes(haplotypes, blues, blocks, n_pcs, top_n, min_freq,
  id_col, blue_col, blue_cols, alpha, verbose)`** — Block-level haplotype
  association tests via a unified Q+K mixed linear model (EMMAX/GAPIT3
  formulation): y = mu + alpha·x_hap + sum(beta_k·PC_k) + g + e, where PC_k
  are GRM-derived eigenvectors (fixed effects for population structure) and g
  ~ MVN(0, sigma_g^2 G) is the polygenic kinship random effect. GRM is inverted
  once per trait via `rrBLUP::mixed.solve()` (O(n^3)); per-allele scan across
  all blocks is fully vectorised in a single `crossprod()` call (O(n*p), same
  BLAS trick as marginal SNP screening). Returns a `LDxBlocks_haplotype_assoc`
  object with `$allele_tests` (per-allele Wald tests: effect, SE, t, p_wald,
  p_wald_adj, significant) and `$block_tests` (omnibus F-test per block:
  F_stat, p_omnibus, p_omnibus_adj, var_explained, significant_omnibus).
  `n_pcs = 0L` (default): EMMAX pure GRM; `n_pcs > 0`: Q+K model;
  `n_pcs = NULL`: auto-select via scree elbow.

- **`estimate_diplotype_effects(haplotypes, blues, blocks, min_freq,
  min_n_diplotype, id_col, blue_col, blue_cols, verbose)`** — Estimates
  additive and dominance effects from diplotype class means at each LD block,
  after GRM kinship correction. For each allele pair (A, B): additive effect
  a = (mean_BB - mean_AA) / 2; dominance deviation d = mean_AB - midpoint;
  dominance ratio d/a (0 = additive, +/-1 = complete dominance, |d/a| > 1 =
  overdominance). Returns a `LDxBlocks_diplotype` object with
  `$diplotype_means`, `$dominance_table` (a, d, d_over_a, overdominance per
  allele pair per block per trait), and `$omnibus_tests` (F-test per block).

**Breeding decision tools**

- **`score_favorable_haplotypes(haplotypes, allele_effects, min_freq,
  missing_string, normalize)`** — Scores each individual's genome-wide
  haplotype portfolio against a table of known per-allele effects. Block score
  = sum(allele_effect × dosage) per block. Genome-wide stacking index = sum
  across all scored blocks, normalised to [0,1] when `normalize = TRUE`.
  Returns a ranked data frame (one row per individual) with `stacking_index`,
  `n_blocks_scored`, `mean_block_score`, `rank`, and one `score_<block_id>`
  column per scored block for detailed inspection.

- **`summarize_parent_haplotypes(haplotypes, candidate_ids, allele_effects,
  blocks, min_freq, missing_string)`** — Produces a tidy long-format allele
  inventory (one row per individual × block × allele) for candidate parents.
  Reports allele dosage (0/1 unphased; 0/1/2 phased), population allele
  frequency, optional allele effect, and a `is_rare` flag (freq < 0.10).
  Includes rows with dosage = 0 so all candidates can be compared on the same
  rows. Primary tool for identifying complementary rare alleles across
  candidates and designing haplotype stacking crosses.

### New example dataset


- **`ldx_blues_list`** — named list of two environments (`env1`, `env2`) of
  named numeric BLUEs (120 individuals each), generated from the same
  polygenic architecture as `ldx_blues` with environment-specific offsets.
  Used in examples for `run_haplotype_stability()` and
  `cv_haplotype_prediction()`. Flat-file copy at
  `inst/extdata/example_blues_env.csv`.

### Analysis extension functions (v0.3.1 additions)

Ten new exported functions extend the pipeline without modifying any existing
function. All are in two new files: `R/analysis_extensions.R` and
`R/haplotype_inference.R`.

**Cross-validation and model evaluation**

- **`cv_haplotype_prediction()`** — k-fold cross-validation for the haplotype
  GBLUP model. Masks phenotypes fold-by-fold, predicts via `rrBLUP::kin.blup()`
  using the shared haplotype GRM, and returns predictive ability (Pearson r) and
  RMSE per trait per fold. Supports multiple replications and multiple traits.
  Returns an `LDxBlocks_cv` object with `pa_summary`, `pa_mean`, `k`, `n_rep`.

**Population comparison**

- **`compare_haplotype_populations()`** — computes Weir-Cockerham (1984) FST
  and allele frequency differences per block between two named sample groups.
  Returns `FST`, `max_freq_diff`, dominant allele per group, chi-squared
  p-value (Monte Carlo, B=2000), and a `divergent` flag (FST > 0.1 AND p < 0.05).
  Suitable for breeding cycle monitoring and wild/elite panel comparisons.

**Visualisation**

- **`plot_haplotype_network()`** — draws a minimum-spanning network of haplotype
  alleles within one LD block using `igraph::mst()` with Hamming-distance
  edge weights. Node size proportional to frequency. Optional group colouring
  via a `groups` named vector. Returns the `igraph` MST object invisibly.

**Multi-environment stability**

- **`run_haplotype_stability()`** — Finlay-Wilkinson (1963) regression of
  per-block local GEBV contributions against the environmental index across
  environments. Returns slope b (stability coefficient), SE, R², deviation
  mean square s²d, and a `stable` flag (H0: b=1 not rejected at alpha=0.05).
  Requires at least 2 environments.

**Annotation export**

- **`export_candidate_regions()`** — converts `define_qtl_regions()` output to
  BED (0-based, UCSC/BEDtools-compatible), CSV, or a named list ready for
  `biomaRt::getBM()`. Supports `chr_prefix`, LD-extended windows via
  `use_lead_snp`, and `padding_bp`.

**Effect decomposition**

- **`decompose_block_effects()`** — aggregates per-SNP additive effects (from
  `backsolve_snp_effects()`) into a per-haplotype-allele effect table.
  Effect = sum(SNP_effect × allele_dosage) per allele position. Returns
  `allele_effect`, `effect_rank`, and `frequency` per allele per block.
  Directly links prediction model output to selection index construction.

**Genome-wide diversity scanning**

- **`scan_diversity_windows()`** — sliding-window He / Shannon / n_eff_alleles
  scan across the genome independent of LD block boundaries. Window size and
  step controlled by `window_bp` and `step_bp`. Returns a data frame with one
  row per window including `sweep_flag` (freq_dominant >= 0.90).

**True haplotype inference and harmonisation**

- **`infer_block_haplotypes()`** — converts raw haplotype strings to a
  structured per-individual, per-block diplotype table with explicit `hap1`,
  `hap2`, `diplotype` (canonical sorted string), `heterozygous`,
  `phase_ambiguous`, and `missing` columns. Handles both phased input
  (from `read_phased_vcf()`, where `phase_ambiguous`
  is always `FALSE`) and unphased input (where heterozygous genotypes set
  `phase_ambiguous = TRUE` unless `resolve_unphased = TRUE` triggers a
  maximum-parsimony heuristic).

- **`collapse_haplotypes()`** — merges rare haplotype alleles (below
  `min_freq`) into biologically meaningful groups rather than dropping them.
  Three strategies: `"rare_to_other"` (pool into `<other>`), `"nearest"`
  (merge with most similar common allele by Hamming distance), `"tree_based"`
  (UPGMA dendrogram; merges rare alleles at the coarsest cut that avoids
  merging common alleles with each other). Preserves a `label_map` attribute
  recording every original→collapsed mapping for use by `harmonize_haplotypes()`.

- **`harmonize_haplotypes()`** — makes haplotype allele labels transferable
  across training/validation splits, populations, or environments. Builds a
  reference dictionary from alleles above `min_freq_ref` in the reference
  panel; matches target alleles by exact string first, then nearest Hamming
  neighbour (up to `max_hamming`), then labels unmatched alleles `"<novel>"`.
  Attaches a `harmonization_report` attribute reporting `n_exact`, `n_nearest`,
  `n_novel`, and `mean_hamming_dist` per block.

### New functions

- **`compute_ld_decay()`** -- LD decay analysis per chromosome. Estimates the
  distance at which r² (or rV²) drops below a critical threshold using a
  memory-efficient position-first random sampling strategy: pair indices are
  sampled from SNP positions only, then `read_chunk()` loads only the unique
  SNP columns involved (~2% of a WGS chromosome for 50k random pairs).
  `compute_r2_sparse_cpp()` (Armadillo BLAS + OpenMP) is called once per
  chromosome in place of an R pair-by-pair loop. Two threshold approaches:
  fixed numeric (e.g. 0.1, standard GWAS practice) and parametric (95th
  percentile of r² between unlinked markers on different chromosomes, measuring
  the background kinship-induced LD level). Optional LOESS and Hill-Weir (1988)
  nonlinear decay model fitting. Chromosome-specific decay distances can be
  passed directly to `define_qtl_regions(ld_decay = decay)` to define
  biologically justified candidate gene windows. Censored distances (threshold
  never crossed within `max_dist`) are flagged with `censored = TRUE` in output
  and emitted as warnings. Requires no change to existing code.

- **`plot_ld_decay()`** -- ggplot2 visualisation of r² vs physical distance
  per chromosome, with optional raw points, threshold line, per-chromosome
  decay distance markers, and facet option.

### Algorithm improvements

- **LD-informed overlap resolution** replacing blind union merge in `Big_LD()`.
  Overlapping blocks at sub-segment seams are now resolved by a cumulative-score
  boundary rule: for each disputed SNP, `score = mean_r2(with left core) -
  mean_r2(with right core)`. The cumulative sum is tracked across the overlap
  zone and the split boundary is placed at the last position where the cumulative
  score is >= 0. Representatives are selected from boundary-adjacent SNPs of
  each core (nearest to the overlap zone, not first-k). Three clean cases:
  all-left (block A keeps overlap, block B shrinks), all-right (block A shrinks,
  block B keeps overlap), mixed (split at cumulative boundary). Falls back to
  union merge only when one core is empty (one block fully inside the other).

### C++ updates

- **`compute_r2_sparse_cpp()`** gains an `n_threads` parameter (OpenMP).
  Thread-local vector accumulation prevents contention; results are merged
  after the parallel loop. The R wrapper in `RcppExports.R` is updated to
  expose the new argument with default `n_threads = 1L`.

### Bug fixes and robustness

- **`compute_ld_decay()` pair sampling** -- `snp_info` is now sorted by POS
  within each chromosome before any pair index sampling. Previously, unsorted
  input could produce incorrect pair distances silently.
- **`compute_ld_decay()` column precomputation** -- prepared columns (mean-imputed,
  optionally whitened) are now computed once per chunk, not per pair. Monomorphic
  columns (sd < 1e-6 after preparation) are skipped before any r² computation.
- **`compute_ld_decay()` `both` mode** -- sliding-window and random pairs are
  now kept separate: sliding window feeds the decay curve, random pairs feed
  the parametric threshold only. Previously they were merged into one pool,
  biasing the curve shape estimate.
- **`compute_ld_decay()` fitting failures** -- LOESS and Hill-Weir nonlinear
  failures now emit `warning()` with chromosome name and error message instead
  of being silently swallowed by `tryCatch(..., error = function(e) NULL)`.
- **`compute_ld_decay()` threshold crossing** -- linear interpolation between
  the last-above and first-below bins for sub-bin precision (was: first-below
  bin midpoint). Non-monotone curves fall back to first-below.
- **`r2_threshold` validation** -- numeric thresholds outside [0,1], NA values,
  and unrecognised strings now produce informative errors instead of silent
  wrong results.

### Performance improvements

- **`extract_haplotypes()` C++ optimisation** — haplotype strings now built
  via `build_hap_strings_cpp()`, replacing an R `vapply()` loop that incurred
  one function call per individual per block. For a 3M-SNP panel with 17,078
  blocks and 204 individuals: ~3.5 million R calls -> one C++ call per block.
  Expected speedup: 20-50x (3.5 hours -> ~5-10 minutes for that panel).
- **`read_chunk()` bigmemory NA fix** — char-type big.matrix stores NA as -128;
  these are now correctly restored to `NA_integer_`. Also removed a redundant
  integer->double type conversion.

### WGS-scale acceleration — three new approaches

- **Louvain/Leiden community detection** (`CLQmode = "Louvain"` or `"Leiden"`):
  Polynomial-time O(n log n) community detection replaces Bron-Kerbosch
  (`igraph::max_cliques()`) which has exponential worst-case complexity.
  The original Big-LD run on a 3M-SNP WGS panel found 4.26 million maximal
  cliques in a single 1500-SNP window and ran for > 1 hour without completing.
  Louvain/Leiden finish the same window in < 1 second. Block boundaries are
  equivalent to or better than the Density mode for WGS panels where LD
  structure is dense. Available in `CLQD()`, `Big_LD()`, `run_Big_LD_all_chr()`,
  `tune_LD_params()`, and `run_ldx_pipeline()`. Edge weights are inversely
  proportional to base-pair distance so local LD structure is respected.

- **Sparse LD computation** (`max_bp_distance` parameter):
  When `max_bp_distance > 0`, only SNP pairs within that physical distance
  have their r² computed via `compute_r2_sparse_cpp()`. Pairs beyond the
  threshold are assumed to be in negligible LD and set to zero in the
  adjacency matrix. This reduces O(p²) LD computation to near-O(p) for
  WGS panels — at 500 kb, approximately 70–90% of pairs are skipped.
  Available in `CLQD()`, propagated through all callers. Default `0L`
  (disabled) preserves original behaviour.

- **Memory-mapped genotype store** (`read_geno_bigmemory()`):
  Wraps any genotype source in a `bigmemory::big.matrix` backed by a binary
  file on disk. `read_chunk()` retrieves columns via OS page faults — only
  the requested bytes are loaded into RAM. Peak RAM is proportional to
  `n_samples × subSegmSize × 8 bytes`, not the full genome matrix.
  Backing files persist across R sessions: supply the `.desc` path to
  reattach without re-loading source data. Storage type `"char"` (1 byte per
  cell) saves 8× RAM vs double for 0/1/2 dosage values.
  Requires `bigmemory` (added to Suggests).

### Recommended configuration for 3M-SNP WGS panels

```r
blocks <- run_Big_LD_all_chr(
  be,
  CLQmode         = "Leiden",    # polynomial — guaranteed connected communities
  CLQcut          = 0.70,        # sparser LD graph
  max_bp_distance = 500000L,     # skip pairs > 500 kb
  subSegmSize     = 500L,        # smaller windows for safety
  leng            = 50L,         # narrow boundary scan for WGS density
  checkLargest    = TRUE,        # belt-and-suspenders for Density mode
  n_threads       = n_threads
)
```

---

## LDxBlocks 0.3.0

### Breaking changes

- `read_geno()`: removed `pheno_ids` parameter. Sample subsetting is the
  caller's responsibility -- use `be$sample_ids` to match externally.
- Removed `read_pheno()` and `align_geno_pheno()` entirely. Phenotype handling
  is not part of LD block detection and these functions had no algorithmic role.
- `run_Big_LD_all_chr()`: parameter `rV2method` renamed to `kin_method`.
- `digits` default changed from `6` to `-1` (no rounding) throughout.
- SeqArray backend replaced by SNPRelate throughout. GDS files are now created
  via `SNPRelate::snpgdsVCF2GDS()` and read via `SNPRelate::snpgdsGetGeno()`.
  Existing `.gds` files created by SeqArray are not compatible; delete and
  re-convert with `read_geno(..., verbose = TRUE)`.

### New features -- haplotype module (comprehensive rewrite)

- **Phasing functions**:
  - `read_phased_vcf()`: read pre-phased VCF with `0|1` GT fields; returns
    `hap1`/`hap2` gamete matrices plus combined dosage.
  - `phase_with_beagle()`: call Beagle 5.x for statistical phasing of WGS
    data; returns path to phased VCF.gz.
  - `unphase_to_dosage()`: collapse phased gamete matrices back to 0/1/2
    (internal helper; not exported).
- **Haplotype extraction** (`extract_haplotypes()`):
  - Auto-detects phased vs unphased input from list structure.
  - Unphased mode: diploid string `"012201"` per individual per block.
  - Phased mode: two-gamete string `"011|100"` per individual per block.
  - Blocks are strictly per-chromosome; cross-chromosome blocks are
    architecturally impossible.
- **Enhanced diversity metrics** (`compute_haplotype_diversity()`):
  - `n_eff_alleles`: effective number of alleles = 1/Σpᵢ² (Hill 1973),
    ranging from 1 (monomorphic) to k (equal frequencies).
  - `sweep_flag`: TRUE when `freq_dominant ≥ 0.90`, flagging possible
    selective sweeps or strong founder effects (Difabachew et al. 2023).
  - `He` is now sample-size corrected following Nei (1973):
    He = n/(n-1) × (1 − Σpᵢ²).
  - Output now includes `CHR`, `start_bp`, `end_bp`, `n_snps` columns
    for direct use in genomic region analyses.

- **Haplotype string decoder** (`decode_haplotype_strings()`):
  - Converts raw dosage strings (e.g. `"02110"`) to nucleotide sequences
    (e.g. `"AGTTA"`) using REF/ALT from `snp_info`.
  - Returns a data frame: block_id, CHR, start_bp, end_bp, hap_rank, hap_id,
    dosage_string, nucleotide_sequence, frequency, n_carriers, snp_positions,
    snp_alleles.
- **Feature matrix** (`build_haplotype_feature_matrix()`):
  - New `encoding =` parameter: `"additive_012"` (default) or `"presence_01"`.
  - Phased + `"additive_012"`: true 0/1/2 allele counts per gamete (0=absent,
    1=one gamete, 2=both gametes).
  - Unphased + `"additive_012"`: 0/1/NA — 1=present, 0=absent. The value 2
    is not used because homozygosity cannot be inferred from unphased strings.
  - `"presence_01"`: 0/1 presence/absence for kernel methods or random forests
    (formerly `"presence_02"`; `"presence_02"` accepted as a backward-compat alias).
  - New `min_freq =` parameter to drop rare haplotype allele columns.
  - `top_n` default changed from `5L` to `NULL` (retain all alleles above
    `min_freq`).
- **QTL region definition** (`define_qtl_regions()`):
  - Maps significant GWAS markers onto LD blocks post-GWAS.
  - Flags pleiotropic blocks (significant hits from multiple traits).
  - Implements the haploblock-based QTL cataloguing approach of
    Tong et al. (2024) *Theoretical and Applied Genetics* 137:274.
  - Now accepts optional `BETA` column in `gwas_results`. When supplied,
    output includes: `lead_beta` (effect of lead SNP), `sig_snps`
    (all significant SNP IDs, semicolon-separated), `sig_betas` (their
    marginal effects). See `?define_qtl_regions` for block effect
    estimation approaches.
- **Output writers** (all produce rows = haplotype alleles, cols = individuals,
  with metadata columns hap_id, CHR, start_bp, end_bp, n_snps before individual
  columns):
  - `write_haplotype_numeric()`: tab-delimited dosage matrix (0/1/2/NA).
    Metadata column `alleles` shows the nucleotide sequence of this allele.
    Compatible with rrBLUP, BGLR, ASReml-R.
  - `write_haplotype_character()`: tab-delimited nucleotide matrix. Each cell
    shows the nucleotide sequence if the individual carries that allele, `"-"`
    if absent, `"."` if missing. Heterozygous SNP positions are encoded with
    IUPAC ambiguity codes (R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T, M=A/C),
    keeping the sequence the same length as `n_snps`. `Alleles` column shows
    all REF/ALT at every block SNP (semicolon-separated).
  - `write_haplotype_diversity()`: CSV of per-block diversity metrics with
    optional genome-wide mean summary row.
  - `write_haplotype_hapmap()` is **removed**. The diploid AA/AT/TT HapMap
    encoding was not meaningful for multi-SNP haplotype alleles.
- **Haplotype string decoder** (`decode_haplotype_strings()`): exported function
  mapping dosage strings to nucleotide sequences with full SNP-level metadata.
- **Multi-trait haplotype prediction** (extended `run_haplotype_prediction()`):
  - Extends `run_haplotype_prediction()` to k traits simultaneously.
  - One trait-agnostic GRM (computed once) shared across all traits.
  - GBLUP solver: attempts `sommer::mmer()` multi-trait model first;
    automatically falls back to `rrBLUP::kin.blup()` per-trait loop if
    sommer is not installed or the model fails to converge.
  - Block importance aggregated across traits: `var_scaled_mean`,
    `n_traits_important`, `important_any`, `important_all` per block.
  - `importance_rule` argument controls combined flag:
    `'any'` (≥ 1 trait), `'all'` (all traits), `'mean'` (mean ≥ 0.9).
  - `blues` accepts a wide data frame (one column per trait) or a named
    list of named numeric vectors (different individuals per trait).
  - `sommer` added to `Suggests`; `rrBLUP` is the required fallback.
  - `run_haplotype_prediction_mt()` is removed; its functionality is
    fully absorbed into the unified `run_haplotype_prediction()` API.
  - Cites Covarrubias-Pazaran (2016) for sommer and Endelman (2011) for rrBLUP.

- **Haplotype-based genomic prediction** (Tong et al. 2024, 2025):
  - `compute_haplotype_grm()`: VanRaden (2008) GRM from haplotype feature
    matrix. Mean-imputes NA, clamps frequencies, returns symmetric n×n matrix.
  - `backsolve_snp_effects()`: derives per-SNP additive effects from GEBV
    without refitting the marker model: α̂ = M′G⁻¹ĝ / 2Σp(1−p)
    (Tong et al. 2025 *Theor Appl Genet* 138:267).
  - `compute_local_gebv()`: local haplotype GEBV per block = sum of SNP
    effects within block. Ranks blocks by Var(local GEBV); `important` flag
    when scaled variance ≥ 0.90 (Tong et al. 2024 *Theor Appl Genet* 137:274).
  - `prepare_gblup_inputs()`: aligns haplotype feature matrix with a
    phenotype data frame; computes and returns a bended GRM ready for
    rrBLUP, sommer, ASReml-R, or BGLR.
  - `run_haplotype_prediction()`: end-to-end Tong et al. (2024/2025)
    pipeline from pre-adjusted phenotype values (BLUEs/BLUPs) to block
    importance. Accepts `blues` as a named numeric vector or a data frame
    with `id_col` and `blue_col` arguments. Uses `rrBLUP::kin.blup()` for
    REML-based GBLUP. Returns all standard pipeline outputs plus GEBV,
    per-SNP effects, local haplotype GEBV matrix, and block importance table.
  - `integrate_gwas_haplotypes()`: combines GWAS evidence (`has_gwas_hit`),
    variance evidence (`is_important`, scaled Var(local GEBV) ≥ 0.9), and
    diversity evidence (`is_diverse`, He ≥ threshold) into a `priority_score`
    (0–3) with plain-language `recommendation` per block.
  - `rank_haplotype_blocks()`: unified block ranking across three use cases:
    (1) diversity-only — rank by He; (2) GWAS-only — binary GWAS-hit flag
    then He as tiebreaker, p-value not used for ranking; (3) phenotype —
    rank by scaled Var(local GEBV). Returns all standard pipeline outputs
    plus `ranked_blocks` data frame.

- **End-to-end pipeline** (`run_ldx_pipeline()`):
  - `hap_format = "hapmap"` renamed to `hap_format = "character"`.
  - All Big_LD arguments now exposed: `clstgap`, `split`, `appendrare`,
    `singleton_as_block`, `checkLargest`, `digits`, `kin_method`, `CLQmode`.
    Previously 8 of these were silently using Big_LD defaults.
  - `min_freq` parameter exposed (was previously hardcoded inside
    `build_haplotype_feature_matrix`).
  - Return list now includes `geno_matrix` (individuals × SNPs, MAF-filtered),
    enabling direct use with `tune_LD_params()` and `run_haplotype_prediction()`
    without reloading the genotype file.
  - `tune_LD_params()` default grid now jointly optimises `CLQcut` (4 values)
    and `min_freq` (2 values), giving 8 combinations. Weber et al. (2023)
    show both are hyperparameters requiring dataset-specific tuning.

### New features -- I/O and backend

- **Multi-format I/O**: `read_geno()` supports numeric dosage CSV, HapMap,
  VCF / bgzipped VCF, SNPRelate GDS, PLINK BED/BIM/FAM, and plain R matrices.
- **SNPRelate GDS streaming backend**: peak RAM proportional to one window.
  Replaces the previous SeqArray backend entirely.
- **PLINK BED backend**: `BEDMatrix`-backed memory-mapped access.
- **`read_chunk()`**, **`close_backend()`**: unified accessor interface.
- **`run_Big_LD_all_chr()`** accepts `LDxBlocks_backend` directly.
- **`min_snps_chr`** parameter: skip scaffolds with fewer SNPs (default 10).
- **`clean_malformed`** parameter in `read_geno()`, `run_Big_LD_all_chr()`,
  and `run_ldx_pipeline()`: streams input files to remove malformed lines
  before GDS conversion (for NGSEP and some GATK-origin VCFs).
- **`chr`** parameter in `run_Big_LD_all_chr()` and `run_ldx_pipeline()`:
  restrict processing to specific chromosomes.
- **`ldx_gwas`** example dataset now includes a `trait` column
  (`"TraitA"` / `"TraitB"`) to demonstrate multi-trait pleiotropic block
  detection via `define_qtl_regions(trait_col = "trait")`.

### Performance improvements -- C++ core

- `boundary_scan_cpp()`: C++ boundary scan replaces R inner loop (~20x faster).
- `maf_filter_cpp()`: combined MAF + monomorphic filter, single O(np) pass
  (~10x faster for panels > 100 k SNPs).
- `build_adj_matrix_cpp()`: in-place adjacency construction, eliminates
  intermediate allocation.
- `compute_r2_sparse_cpp()`: sparse r² within bp window, avoids O(p²) cost
  for large sub-segments.

### Performance improvements -- never-full-genome memory model

- **Chunked pre-allocated numeric CSV reader**: two-pass strategy -- header
  scan only in Pass 1, then fixed 50,000-row chunks fill a single
  pre-allocated matrix in Pass 2. Peak RAM = one chunk, not 2x the file.
- **VCF and HapMap mandatory SNPRelate GDS auto-conversion**: transparent
  streaming GDS backend on first call; cache reused on subsequent calls.
- **Explicit per-chromosome gc()** in `run_Big_LD_all_chr()` and
  `extract_haplotypes()`: prevents heap fragmentation across 20-30 chromosome
  passes.

### Tests

- `test-ld-decay.R`: 33 tests for `compute_ld_decay()` and `plot_ld_decay()`
  covering: input validation, `LDxBlocks_decay` object structure, all sampling
  modes, unsorted `snp_info` handling, all threshold types (fixed, parametric,
  both, NULL), `fit_model` options, decay distance correctness, censored flag,
  all three backend types (matrix, `LDxBlocks_backend`, bigmemory), print and
  plot methods, and `define_qtl_regions()` integration with `ld_decay=`.
- `test-basic.R`: smoke tests covering all major functions via the 230-SNP
  example dataset (3 chromosomes, 9 LD blocks).
- `test-algorithm.R`: property-based tests for the Big-LD segmentation
  algorithm: `CLQD`, `Big_LD`, `run_Big_LD_all_chr`, `summarise_blocks`,
  `plot_ld_blocks`.
- `test-haplotypes.R`: comprehensive tests for haplotype extraction (phased
  and unphased), diversity metrics (`n_eff_alleles`, `sweep_flag`), QTL
  region definition (`lead_beta`, `sig_snps`, `sig_betas`), feature matrix
  encoding (`additive_012` / `presence_01`), output writers,
  `rank_haplotype_blocks`, and `integrate_gwas_haplotypes`.

---

## LDxBlocks 0.2.0 (2025-04-07)

### New features

- C++/Armadillo core: `compute_r2_cpp()`, `compute_rV2_cpp()`,
  `maf_filter_cpp()`, `build_adj_matrix_cpp()`, `col_r2_cpp()`,
  `compute_r2_sparse_cpp()`, `boundary_scan_cpp()`.
- OpenMP parallelism via `n_threads` parameter.
- Dual LD metric: `method = "r2"` or `"rV2"`.
- `prepare_geno()`, `compute_ld()`.

### Bug fixes

- `CLQD()`: fixed `re_idx` tracking after dense-core pre-pass.
- `Big_LD()`: `vapply()` used throughout for type-stability.

---

## LDxBlocks 0.1.0 (2025-04-07)

Initial release.

- `Big_LD()`, `CLQD()`, `run_Big_LD_all_chr()`, `tune_LD_params()`.
- `extract_haplotypes()`, `compute_haplotype_diversity()`,
  `build_haplotype_feature_matrix()`.
- `summarise_blocks()`, `plot_ld_blocks()`.

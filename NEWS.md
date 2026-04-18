## LDxBlocks 0.3.1 (development)

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
  (from `read_phased_vcf()` / `phase_with_pedigree()`, where `phase_ambiguous`
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
  - `phase_with_pedigree()`: Mendelian allele transmission phasing within
    parent-offspring trios; exact when parents are homozygous.
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

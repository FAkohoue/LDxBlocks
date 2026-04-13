# Changelog

## LDxBlocks 0.3.1 (development)

### WGS-scale acceleration — three new approaches

- **Louvain/Leiden community detection** (`CLQmode = "Louvain"` or
  `"Leiden"`): Polynomial-time O(n log n) community detection replaces
  Bron-Kerbosch
  ([`igraph::max_cliques()`](https://r.igraph.org/reference/cliques.html))
  which has exponential worst-case complexity. The original Big-LD run
  on a 3M-SNP WGS panel found 4.26 million maximal cliques in a single
  1500-SNP window and ran for \> 1 hour without completing.
  Louvain/Leiden finish the same window in \< 1 second. Block boundaries
  are equivalent to or better than the Density mode for WGS panels where
  LD structure is dense. Available in
  [`CLQD()`](https://FAkohoue.github.io/LDxBlocks/reference/CLQD.md),
  [`Big_LD()`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md),
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
  [`tune_LD_params()`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md),
  and
  [`run_ldx_pipeline()`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md).
  Edge weights are inversely proportional to base-pair distance so local
  LD structure is respected.

- **Sparse LD computation** (`max_bp_distance` parameter): When
  `max_bp_distance > 0`, only SNP pairs within that physical distance
  have their r² computed via `compute_r2_sparse_cpp()`. Pairs beyond the
  threshold are assumed to be in negligible LD and set to zero in the
  adjacency matrix. This reduces O(p²) LD computation to near-O(p) for
  WGS panels — at 500 kb, approximately 70–90% of pairs are skipped.
  Available in
  [`CLQD()`](https://FAkohoue.github.io/LDxBlocks/reference/CLQD.md),
  propagated through all callers. Default `0L` (disabled) preserves
  original behaviour.

- **Memory-mapped genotype store**
  ([`read_geno_bigmemory()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno_bigmemory.md)):
  Wraps any genotype source in a
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  backed by a binary file on disk.
  [`read_chunk()`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md)
  retrieves columns via OS page faults — only the requested bytes are
  loaded into RAM. Peak RAM is proportional to
  `n_samples × subSegmSize × 8 bytes`, not the full genome matrix.
  Backing files persist across R sessions: supply the `.desc` path to
  reattach without re-loading source data. Storage type `"char"` (1 byte
  per cell) saves 8× RAM vs double for 0/1/2 dosage values. Requires
  `bigmemory` (added to Suggests).

### Recommended configuration for 3M-SNP WGS panels

``` r
blocks <- run_Big_LD_all_chr(
  be,
  CLQmode         = "Louvain",   # polynomial — no clique blowup
  CLQcut          = 0.70,        # sparser LD graph
  max_bp_distance = 500000L,     # skip pairs > 500 kb
  subSegmSize     = 500L,        # smaller windows for safety
  leng            = 50L,         # narrow boundary scan for WGS density
  checkLargest    = TRUE,        # belt-and-suspenders for Density mode
  n_threads       = n_threads
)
```

------------------------------------------------------------------------

## LDxBlocks 0.3.0

### Breaking changes

- [`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md):
  removed `pheno_ids` parameter. Sample subsetting is the caller’s
  responsibility – use `be$sample_ids` to match externally.
- Removed `read_pheno()` and `align_geno_pheno()` entirely. Phenotype
  handling is not part of LD block detection and these functions had no
  algorithmic role.
- [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md):
  parameter `rV2method` renamed to `kin_method`.
- `digits` default changed from `6` to `-1` (no rounding) throughout.
- SeqArray backend replaced by SNPRelate throughout. GDS files are now
  created via
  [`SNPRelate::snpgdsVCF2GDS()`](https://rdrr.io/pkg/SNPRelate/man/snpgdsVCF2GDS.html)
  and read via
  [`SNPRelate::snpgdsGetGeno()`](https://rdrr.io/pkg/SNPRelate/man/snpgdsGetGeno.html).
  Existing `.gds` files created by SeqArray are not compatible; delete
  and re-convert with `read_geno(..., verbose = TRUE)`.

### New features – haplotype module (comprehensive rewrite)

- **Phasing functions**:
  - [`read_phased_vcf()`](https://FAkohoue.github.io/LDxBlocks/reference/read_phased_vcf.md):
    read pre-phased VCF with `0|1` GT fields; returns `hap1`/`hap2`
    gamete matrices plus combined dosage.
  - [`phase_with_beagle()`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_beagle.md):
    call Beagle 5.x for statistical phasing of WGS data; returns path to
    phased VCF.gz.
  - [`phase_with_pedigree()`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_pedigree.md):
    Mendelian allele transmission phasing within parent-offspring trios;
    exact when parents are homozygous.
  - [`unphase_to_dosage()`](https://FAkohoue.github.io/LDxBlocks/reference/unphase_to_dosage.md):
    collapse phased gamete matrices back to 0/1/2 (internal helper; not
    exported).
- **Haplotype extraction**
  ([`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)):
  - Auto-detects phased vs unphased input from list structure.
  - Unphased mode: diploid string `"012201"` per individual per block.
  - Phased mode: two-gamete string `"011|100"` per individual per block.
  - Blocks are strictly per-chromosome; cross-chromosome blocks are
    architecturally impossible.
- **Enhanced diversity metrics**
  ([`compute_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md)):
  - `n_eff_alleles`: effective number of alleles = 1/Σpᵢ² (Hill 1973),
    ranging from 1 (monomorphic) to k (equal frequencies).
  - `sweep_flag`: TRUE when `freq_dominant ≥ 0.90`, flagging possible
    selective sweeps or strong founder effects (Difabachew et al. 2023).
  - `He` is now sample-size corrected following Nei (1973): He = n/(n-1)
    × (1 − Σpᵢ²).
  - Output now includes `CHR`, `start_bp`, `end_bp`, `n_snps` columns
    for direct use in genomic region analyses.
- **Haplotype string decoder**
  ([`decode_haplotype_strings()`](https://FAkohoue.github.io/LDxBlocks/reference/decode_haplotype_strings.md)):
  - Converts raw dosage strings (e.g. `"02110"`) to nucleotide sequences
    (e.g. `"AGTTA"`) using REF/ALT from `snp_info`.
  - Returns a data frame: block_id, CHR, start_bp, end_bp, hap_rank,
    hap_id, dosage_string, nucleotide_sequence, frequency, n_carriers,
    snp_positions, snp_alleles.
- **Feature matrix**
  ([`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)):
  - New `encoding =` parameter: `"additive_012"` (default) or
    `"presence_01"`.
  - Phased + `"additive_012"`: true 0/1/2 allele counts per gamete
    (0=absent, 1=one gamete, 2=both gametes).
  - Unphased + `"additive_012"`: 0/1/NA — 1=present, 0=absent. The value
    2 is not used because homozygosity cannot be inferred from unphased
    strings.
  - `"presence_01"`: 0/1 presence/absence for kernel methods or random
    forests (formerly `"presence_02"`; `"presence_02"` accepted as a
    backward-compat alias).
  - New `min_freq =` parameter to drop rare haplotype allele columns.
  - `top_n` default changed from `5L` to `NULL` (retain all alleles
    above `min_freq`).
- **QTL region definition**
  ([`define_qtl_regions()`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)):
  - Maps significant GWAS markers onto LD blocks post-GWAS.
  - Flags pleiotropic blocks (significant hits from multiple traits).
  - Implements the haploblock-based QTL cataloguing approach of Tong et
    al. (2024) *Theoretical and Applied Genetics* 137:274.
  - Now accepts optional `BETA` column in `gwas_results`. When supplied,
    output includes: `lead_beta` (effect of lead SNP), `sig_snps` (all
    significant SNP IDs, semicolon-separated), `sig_betas` (their
    marginal effects). See
    [`?define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
    for block effect estimation approaches.
- **Output writers** (all produce rows = haplotype alleles, cols =
  individuals, with metadata columns hap_id, CHR, start_bp, end_bp,
  n_snps before individual columns):
  - [`write_haplotype_numeric()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_numeric.md):
    tab-delimited dosage matrix (0/1/2/NA). Metadata column `alleles`
    shows the nucleotide sequence of this allele. Compatible with
    rrBLUP, BGLR, ASReml-R.
  - [`write_haplotype_character()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_character.md):
    tab-delimited nucleotide matrix. Each cell shows the nucleotide
    sequence if the individual carries that allele, `"-"` if absent,
    `"."` if missing. Heterozygous SNP positions are encoded with IUPAC
    ambiguity codes (R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T, M=A/C), keeping
    the sequence the same length as `n_snps`. `Alleles` column shows all
    REF/ALT at every block SNP (semicolon-separated).
  - [`write_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_diversity.md):
    CSV of per-block diversity metrics with optional genome-wide mean
    summary row.
  - `write_haplotype_hapmap()` is **removed**. The diploid AA/AT/TT
    HapMap encoding was not meaningful for multi-SNP haplotype alleles.
- **Haplotype string decoder**
  ([`decode_haplotype_strings()`](https://FAkohoue.github.io/LDxBlocks/reference/decode_haplotype_strings.md)):
  exported function mapping dosage strings to nucleotide sequences with
  full SNP-level metadata.
- **Multi-trait haplotype prediction** (extended
  [`run_haplotype_prediction()`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)):
  - Extends
    [`run_haplotype_prediction()`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)
    to k traits simultaneously.
  - One trait-agnostic GRM (computed once) shared across all traits.
  - GBLUP solver: attempts
    [`sommer::mmer()`](https://rdrr.io/pkg/sommer/man/mmer.html)
    multi-trait model first; automatically falls back to
    [`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html)
    per-trait loop if sommer is not installed or the model fails to
    converge.
  - Block importance aggregated across traits: `var_scaled_mean`,
    `n_traits_important`, `important_any`, `important_all` per block.
  - `importance_rule` argument controls combined flag: `'any'` (≥ 1
    trait), `'all'` (all traits), `'mean'` (mean ≥ 0.9).
  - `blues` accepts a wide data frame (one column per trait) or a named
    list of named numeric vectors (different individuals per trait).
  - `sommer` added to `Suggests`; `rrBLUP` is the required fallback.
  - `run_haplotype_prediction_mt()` is removed; its functionality is
    fully absorbed into the unified
    [`run_haplotype_prediction()`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)
    API.
  - Cites Covarrubias-Pazaran (2016) for sommer and Endelman (2011) for
    rrBLUP.
- **Haplotype-based genomic prediction** (Tong et al. 2024, 2025):
  - [`compute_haplotype_grm()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_grm.md):
    VanRaden (2008) GRM from haplotype feature matrix. Mean-imputes NA,
    clamps frequencies, returns symmetric n×n matrix.
  - [`backsolve_snp_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/backsolve_snp_effects.md):
    derives per-SNP additive effects from GEBV without refitting the
    marker model: α̂ = M′G⁻¹ĝ / 2Σp(1−p) (Tong et al. 2025 *Theor Appl
    Genet* 138:267).
  - [`compute_local_gebv()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_local_gebv.md):
    local haplotype GEBV per block = sum of SNP effects within block.
    Ranks blocks by Var(local GEBV); `important` flag when scaled
    variance ≥ 0.90 (Tong et al. 2024 *Theor Appl Genet* 137:274).
  - [`prepare_gblup_inputs()`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_gblup_inputs.md):
    aligns haplotype feature matrix with a phenotype data frame;
    computes and returns a bended GRM ready for rrBLUP, sommer,
    ASReml-R, or BGLR.
  - [`run_haplotype_prediction()`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md):
    end-to-end Tong et al. (2024/2025) pipeline from pre-adjusted
    phenotype values (BLUEs/BLUPs) to block importance. Accepts `blues`
    as a named numeric vector or a data frame with `id_col` and
    `blue_col` arguments. Uses
    [`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html)
    for REML-based GBLUP. Returns all standard pipeline outputs plus
    GEBV, per-SNP effects, local haplotype GEBV matrix, and block
    importance table.
  - [`integrate_gwas_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/integrate_gwas_haplotypes.md):
    combines GWAS evidence (`has_gwas_hit`), variance evidence
    (`is_important`, scaled Var(local GEBV) ≥ 0.9), and diversity
    evidence (`is_diverse`, He ≥ threshold) into a `priority_score`
    (0–3) with plain-language `recommendation` per block.
  - [`rank_haplotype_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/rank_haplotype_blocks.md):
    unified block ranking across three use cases:
    1.  diversity-only — rank by He; (2) GWAS-only — binary GWAS-hit
        flag then He as tiebreaker, p-value not used for ranking; (3)
        phenotype — rank by scaled Var(local GEBV). Returns all standard
        pipeline outputs plus `ranked_blocks` data frame.
- **End-to-end pipeline**
  ([`run_ldx_pipeline()`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md)):
  - `hap_format = "hapmap"` renamed to `hap_format = "character"`.
  - All Big_LD arguments now exposed: `clstgap`, `split`, `appendrare`,
    `singleton_as_block`, `checkLargest`, `digits`, `kin_method`,
    `CLQmode`. Previously 8 of these were silently using Big_LD
    defaults.
  - `min_freq` parameter exposed (was previously hardcoded inside
    `build_haplotype_feature_matrix`).
  - Return list now includes `geno_matrix` (individuals × SNPs,
    MAF-filtered), enabling direct use with
    [`tune_LD_params()`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md)
    and
    [`run_haplotype_prediction()`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)
    without reloading the genotype file.
  - [`tune_LD_params()`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md)
    default grid now jointly optimises `CLQcut` (4 values) and
    `min_freq` (2 values), giving 8 combinations. Weber et al. (2023)
    show both are hyperparameters requiring dataset-specific tuning.

### New features – I/O and backend

- **Multi-format I/O**:
  [`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
  supports numeric dosage CSV, HapMap, VCF / bgzipped VCF, SNPRelate
  GDS, PLINK BED/BIM/FAM, and plain R matrices.
- **SNPRelate GDS streaming backend**: peak RAM proportional to one
  window. Replaces the previous SeqArray backend entirely.
- **PLINK BED backend**: `BEDMatrix`-backed memory-mapped access.
- **[`read_chunk()`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md)**,
  **[`close_backend()`](https://FAkohoue.github.io/LDxBlocks/reference/close_backend.md)**:
  unified accessor interface.
- **[`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)**
  accepts `LDxBlocks_backend` directly.
- **`min_snps_chr`** parameter: skip scaffolds with fewer SNPs (default
  10).
- **`clean_malformed`** parameter in
  [`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md),
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
  and
  [`run_ldx_pipeline()`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md):
  streams input files to remove malformed lines before GDS conversion
  (for NGSEP and some GATK-origin VCFs).
- **`chr`** parameter in
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
  and
  [`run_ldx_pipeline()`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md):
  restrict processing to specific chromosomes.
- **`ldx_gwas`** example dataset now includes a `trait` column
  (`"TraitA"` / `"TraitB"`) to demonstrate multi-trait pleiotropic block
  detection via `define_qtl_regions(trait_col = "trait")`.

### Performance improvements – C++ core

- `boundary_scan_cpp()`: C++ boundary scan replaces R inner loop (~20x
  faster).
- `maf_filter_cpp()`: combined MAF + monomorphic filter, single O(np)
  pass (~10x faster for panels \> 100 k SNPs).
- `build_adj_matrix_cpp()`: in-place adjacency construction, eliminates
  intermediate allocation.
- `compute_r2_sparse_cpp()`: sparse r² within bp window, avoids O(p²)
  cost for large sub-segments.

### Performance improvements – never-full-genome memory model

- **Chunked pre-allocated numeric CSV reader**: two-pass strategy –
  header scan only in Pass 1, then fixed 50,000-row chunks fill a single
  pre-allocated matrix in Pass 2. Peak RAM = one chunk, not 2x the file.
- **VCF and HapMap mandatory SNPRelate GDS auto-conversion**:
  transparent streaming GDS backend on first call; cache reused on
  subsequent calls.
- **Explicit per-chromosome gc()** in
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
  and
  [`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md):
  prevents heap fragmentation across 20-30 chromosome passes.

### Tests

- `test-basic.R`: smoke tests covering all major functions via the
  230-SNP example dataset (3 chromosomes, 9 LD blocks).
- `test-algorithm.R`: property-based tests for the Big-LD segmentation
  algorithm: `CLQD`, `Big_LD`, `run_Big_LD_all_chr`, `summarise_blocks`,
  `plot_ld_blocks`.
- `test-haplotypes.R`: comprehensive tests for haplotype extraction
  (phased and unphased), diversity metrics (`n_eff_alleles`,
  `sweep_flag`), QTL region definition (`lead_beta`, `sig_snps`,
  `sig_betas`), feature matrix encoding (`additive_012` /
  `presence_01`), output writers, `rank_haplotype_blocks`, and
  `integrate_gwas_haplotypes`.

------------------------------------------------------------------------

## LDxBlocks 0.2.0 (2025-04-07)

### New features

- C++/Armadillo core: `compute_r2_cpp()`, `compute_rV2_cpp()`,
  `maf_filter_cpp()`, `build_adj_matrix_cpp()`, `col_r2_cpp()`,
  `compute_r2_sparse_cpp()`, `boundary_scan_cpp()`.
- OpenMP parallelism via `n_threads` parameter.
- Dual LD metric: `method = "r2"` or `"rV2"`.
- [`prepare_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_geno.md),
  [`compute_ld()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_ld.md).

### Bug fixes

- [`CLQD()`](https://FAkohoue.github.io/LDxBlocks/reference/CLQD.md):
  fixed `re_idx` tracking after dense-core pre-pass.
- [`Big_LD()`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md):
  [`vapply()`](https://rdrr.io/r/base/lapply.html) used throughout for
  type-stability.

------------------------------------------------------------------------

## LDxBlocks 0.1.0 (2025-04-07)

Initial release.

- [`Big_LD()`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md),
  [`CLQD()`](https://FAkohoue.github.io/LDxBlocks/reference/CLQD.md),
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
  [`tune_LD_params()`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md).
- [`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
  [`compute_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md),
  [`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md).
- [`summarise_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/summarise_blocks.md),
  [`plot_ld_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/plot_ld_blocks.md).

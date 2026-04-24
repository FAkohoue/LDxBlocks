# Statistical Phasing via Beagle 5.x

Wrapper around Beagle 5.x for statistical haplotype phasing. Calls
`java -jar beagle.jar` and writes the phased VCF.gz to disk. This is the
recommended phasing method in LDxBlocks: Beagle performs
chromosome-level statistical phasing using population-LD across all
markers simultaneously, producing true inferred haplotypes suitable for
block-level frequency estimation and diversity analysis.

## Usage

``` r
phase_with_beagle(
  input_vcf,
  out_prefix,
  beagle_jar = NULL,
  java_path = "java",
  java_mem_gb = NULL,
  nthreads = 1L,
  ref_panel = NULL,
  map_file = NULL,
  chrom = NULL,
  seed = NULL,
  burnin = NULL,
  iterations = NULL,
  window = NULL,
  overlap = NULL,
  beagle_args = "",
  verbose = TRUE
)
```

## Arguments

- input_vcf:

  Path to input VCF or VCF.gz (unphased).

- out_prefix:

  Output path prefix. Beagle appends `.vcf.gz`.

- beagle_jar:

  Path to `beagle.jar`. Default: searched in `dirname(out_prefix)`, then
  standard locations.

- java_path:

  Java executable. Default `"java"`.

- java_mem_gb:

  Java heap size in GB (e.g. `8` sets `-Xmx8g`). Default `NULL` (uses
  JVM default). Increase for large VCFs.

- nthreads:

  Beagle threads. Default `1L`.

- ref_panel:

  Optional path to phased reference VCF. Default `NULL`.

- map_file:

  Optional genetic map file for more accurate phasing in structured
  populations. Default `NULL`.

- chrom:

  Optional chromosome string passed to Beagle (e.g. `"chr1"` or `"1"`).
  Default `NULL` (all chromosomes).

- seed:

  Integer random seed for reproducibility. Default `NULL`.

- burnin:

  Beagle burn-in iterations. Default `NULL` (Beagle default).

- iterations:

  Beagle phasing iterations. Default `NULL` (Beagle default).

- window:

  Beagle window size. Default `NULL` (Beagle default).

- overlap:

  Beagle window overlap. Default `NULL` (Beagle default).

- beagle_args:

  Additional Beagle arguments string, space-separated. Default `""`.

- verbose:

  Logical. Default `TRUE`.

## Value

Invisibly returns the path to the phased VCF.gz. Beagle stdout and
stderr are written to `out_prefix.log`.

## Note

`Sys.which("beagle.jar")` typically fails to find `.jar` files on PATH.
Supply `beagle_jar` explicitly or place `beagle.jar` next to
`out_prefix`. Download:
<https://faculty.washington.edu/browning/beagle/beagle.html>

## References

Browning et al. (2018) Am J Hum Genet 103:338-348.

# Statistical Phasing via Beagle 5.x

Statistical Phasing via Beagle 5.x

## Usage

``` r
phase_with_beagle(
  input_vcf,
  out_prefix,
  beagle_jar = NULL,
  java_path = "java",
  nthreads = 1L,
  ref_panel = NULL,
  beagle_args = "",
  verbose = TRUE
)
```

## Arguments

- input_vcf:

  Unphased VCF path.

- out_prefix:

  Output prefix (Beagle appends .vcf.gz).

- beagle_jar:

  Path to beagle.jar. Default: auto-detected.

- java_path:

  Java executable. Default "java".

- nthreads:

  Threads. Default 1L.

- ref_panel:

  Optional phased reference VCF. Default NULL.

- beagle_args:

  Additional Beagle arguments string.

- verbose:

  Logical. Default TRUE.

## Value

Path to phased VCF.gz. Load with read_phased_vcf().

## References

Browning et al. (2018) Am J Hum Genet 103:338-348.

# Test enrichment across a cluster of motifs using a background set of sequences

Test for enrichment of any motif within a cluster across a set of
sequences using a background set to derive a NULL hypothesis

## Usage

``` r
testClusterEnrich(
  cl,
  stringset,
  bg,
  var = "iteration",
  model = c("poisson", "hypergeometric", "quasipoisson", "glm_poisson", "iteration"),
  sort_by = c("p", "none"),
  mc.cores = 1,
  prior.count = 1,
  seed = 100,
  ...
)
```

## Arguments

- cl:

  A list of Position Weight Matrices, universalmotifs, with each element
  representing clusters of related matrices

- stringset:

  An XStringSet with equal sequence widths

- bg:

  An XStringSet with the same sequence widths as the test XStringset

- var:

  A column in the mcols element of bg, usually denoting an iteration
  number

- model:

  The model used for analysis

- sort_by:

  Column to sort results by

- mc.cores:

  Passed to [mclapply](https://rdrr.io/r/parallel/mclapply.html)

- prior.count:

  Added to all counts to better manage zero counts in background
  sequences. For analysis under QuasiPoisson models prior counts are
  added as Poisson noise using this value as expected counts

- seed:

  Used for reproducibility when adding Poisson noise

- ...:

  Passed to
  [getPwmMatches](https://smped.github.io/motifTestR/reference/getPwmMatches.md)
  or
  [countPwmMatches](https://smped.github.io/motifTestR/reference/countPwmMatches.md)

## Value

See
[testMotifEnrich](https://smped.github.io/motifTestR/reference/testMotifEnrich.md)

## Details

This extends the analytic methods offered by
[testMotifEnrich](https://smped.github.io/motifTestR/reference/testMotifEnrich.md)
using PWMs grouped into a set of clusters. As with all cluster-level
approaches, hits from multiple PWMs which overlap are counted as a
single hit ensuring that duplicated matches are not double-counted, and
that only individual positions within the sequences are.

## See also

[`makeRMRanges()`](https://smped.github.io/motifTestR/reference/makeRMRanges-methods.md),
[`getClusterMatches()`](https://smped.github.io/motifTestR/reference/getClusterMatches.md),
[`countClusterMatches()`](https://smped.github.io/motifTestR/reference/getClusterMatches.md),
[`testMotifEnrich()`](https://smped.github.io/motifTestR/reference/testMotifEnrich.md)

## Examples

``` r
## Load the example peaks & the sequences
data("ar_er_peaks")
data("ar_er_seq")
sq <- seqinfo(ar_er_peaks)
## Now sample size-matched ranges 10 times larger. In real-world analyses,
## this set should be sampled as at least 1000x larger, ensuring features
## are matched to your requirements. This example masks regions with known N
## content, including centromeres & telomeres
data("hg19_mask")
set.seed(305)
bg_ranges <- makeRMRanges(
  ar_er_peaks, GRanges(sq)[1], exclude = hg19_mask, n_iter = 10
)

## Convert ranges to DNAStringSets
library(BSgenome.Hsapiens.UCSC.hg19)
#> Loading required package: BSgenome
#> Loading required package: BiocIO
#> Loading required package: rtracklayer
genome <- BSgenome.Hsapiens.UCSC.hg19
bg_seq <- getSeq(genome, bg_ranges)

## Test for enrichment of clustered motifs
data("ex_pfm")
cl <- list(A = ex_pfm[1], B = ex_pfm[2:3])
testClusterEnrich(cl, ar_er_seq, bg_seq, model = "poisson")
#>   sequences matches expected enrichment        Z            p          fdr
#> B       849    1089    603.1   1.805671 19.78574 1.243312e-70 2.486625e-70
#> A       849     200     54.3   3.683241 19.77241 4.143439e-52 4.143439e-52
#>   est_bg_rate
#> B   0.7103651
#> A   0.0639576

```

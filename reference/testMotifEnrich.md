# Test motif enrichment using a background set of sequences

Test for motif enrichment within a set of sequences using a background
set to derive a NULL hypothesis

## Usage

``` r
testMotifEnrich(
  pwm,
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

- pwm:

  A Position Weight Matrix or list of PWMs

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

A data.frame with columns: `sequences`, `matches`, `expected`,
`enrichment`, and `p`, with additional columns `Z`, `est_bg_rate`
(Poisson), `odds_ratio` (Hypergeometric) or `Z`, `sd_bg`, `n_iter` and
`iter_p` (Iterations). The numbers of sequences and matches refer to the
test set of sequences, whilst expected is the expected number of matches
under the Poisson or iterative null distribution. The ratio of matches
to expected is given as `enrichment`, along with the Z score and
p-value. Whilst the Z score is only representative under the Poisson
model, it is used to directly estimate p-values under the iterative
approach. Under this latter approach, the sd of the matches found in the
background sequences is also given, along with the number of iterations
and the p-values from permutations testing the one-sided hypothesis
hypothesis for enrichment.

It may also be worth noting that if producing background sequences using
[makeRMRanges](https://smped.github.io/motifTestR/reference/makeRMRanges-methods.md)
with `replace = TRUE` and `force_ol = TRUE`, the iterative model
corresponds to a bootstrap, given that the test sequences will overlap
the background sequences and background ranges are able to be sampled
with replacement.

## Details

This function offers four alternative models for assessing the
enrichment of a motif within a set of sequences, in comparison to a
background set of sequences. Selection of the BG sequences plays an
important role and, in conjunction with the question being addressed,
determines the most appropriate model to use for testing, as described
below.

It should also be noted that the larger the BG set of sequences, the
larger the computational burden, and results can take far longer to
return. For many millions of background sequences, this may run beyond
an hour

### Descriptions of Models and Use Cases

#### Hypergeometric Tests

Hypergeometric tests are best suited to the use case where the test set
of sequences represents a subset of a larger set, with a specific
feature or behaviour, whilst the BG set may be the remainder of the set
without that feature. For example, the test set may represent ChIP-Seq
binding sites where signal changes in response to treatment, whilst the
BG set represents the sites where no changed signal was observed.
Testing is one-sided, for enrichment of motifs within the test set.

Due to these relatively smaller sized datasets, setting model =
"hypergeometric", will generally return results quickly

#### Poisson Tests

This approach requires a set of background sequences which should be
much larger than the test set of sequences. The parameters for a Poisson
model are estimated in a per-sequence manner on the set of BG sequences,
and the observed rate of motif-matches within the test set is then
tested using [poisson.test](https://rdrr.io/r/stats/poisson.test.html).
Testing is two-sided.

This approach assumes that all matches follow a Poisson distribution,
which is often true, but data can also be over-dispersed. Given that
this model can also return results relatively quickly, is it primarily
suitable for data exploration, such as quickly checking for expected
behaviours, but not for final results.

#### Quasi-Poisson Test

The quasipoisson model allows for over-dispersion and will return more
conservative results than using the standard Poisson model. Under the
method currently implemented here, BG sequences should be divided into
blocks (i.e. iterations), identical in size to the test set of
sequences. Model parameters are estimated per iteration across the BG
set of sequences, with the rate of matches in the test set being
compared against these blocks. This ensures more conservative results
that if analysing test and bg sequences as collections of individual
sequences.

It is expected that the BG set will matched for the features of interest
and chosen using
[makeRMRanges](https://smped.github.io/motifTestR/reference/makeRMRanges-methods.md)
with a large number of iterations, e.g. n_iter = 1000. Due to this
parameterisation, quasipoisson approaches can be computationally
time-consuming, as this is effectively an iterative approach. Testing is
two-sided.

#### GLM-Poisson Test

This follows the same approach as the Quasi-Poisson model, relying on
fitting iterations using [`glm()`](https://rdrr.io/r/stats/glm.html).
For this model however, no over-dispersions are estimated and the
underlying family is simply the
[`poisson()`](https://rdrr.io/r/stats/family.html) family

#### Iteration

Setting the model as "iteration" performs a non-parametric analysis,
with the exception of returning Z-scores under the Central Limit
Theorem. Mean and SD of matches is found for each iteration, and used to
return Z scores, with p-values returned from both a Z-test and from
comparing observed values directly to sampled values obtained from the
BG sequences. Sampled values are calculated directly and as such, are
limited in precision.

As for the QuasiPoisson model, a very large number of iterations is
expected to be used, to ensure the CLT holds, again making this a
computationally demanding test. Each iteration/block is expected to be
identically-sized to the test set, and matched for any features as
appropriate using
[`makeRMRanges()`](https://smped.github.io/motifTestR/reference/makeRMRanges-methods.md).

## See also

[`makeRMRanges()`](https://smped.github.io/motifTestR/reference/makeRMRanges-methods.md),
[`getPwmMatches()`](https://smped.github.io/motifTestR/reference/getPwmMatches.md),
[`countPwmMatches()`](https://smped.github.io/motifTestR/reference/countPwmMatches.md)

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
genome <- BSgenome.Hsapiens.UCSC.hg19
bg_seq <- getSeq(genome, bg_ranges)
mcols(bg_seq)$iteraton <-bg_ranges$iteration

## Test for enrichment of the ESR1 motif
data("ex_pfm")
esr1 <- ex_pfm$ESR1
testMotifEnrich(esr1, ar_er_seq, bg_seq, model = "poisson")
#>   sequences matches expected enrichment        Z            p          fdr
#> 1       849      23      0.6   38.33333 28.91828 1.719458e-28 1.719458e-28
#>    est_bg_rate
#> 1 0.0007067138

## Test all motifs
testMotifEnrich(ex_pfm, ar_er_seq, bg_seq, model = "poisson")
#>       sequences matches expected enrichment         Z            p          fdr
#> ZN143       849      22      0.1 220.000000 69.253881 8.085297e-44 4.042648e-43
#> FOXA1       849     114     25.6   4.453125 17.471584 1.333620e-37 3.334051e-37
#> ESR1        849      23      0.6  38.333333 28.918276 1.719458e-28 2.865764e-28
#> ANDR        849       9      1.5   6.000000  6.123724 2.773582e-05 3.466977e-05
#> ZN281       849      15      6.5   2.307692  3.333974 4.459196e-03 4.459196e-03
#>        est_bg_rate
#> ZN143 0.0001177856
#> FOXA1 0.0301531213
#> ESR1  0.0007067138
#> ANDR  0.0017667845
#> ZN281 0.0076560660

```

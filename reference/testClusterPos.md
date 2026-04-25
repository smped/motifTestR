# Test positional bias motifs within a cluster

Test positional bias for all motifs within a given cluster

## Usage

``` r
testClusterPos(
  x,
  stringset,
  binwidth = 10,
  abs = FALSE,
  rc = TRUE,
  min_score = "50%",
  break_ties = "all",
  alt = c("greater", "less", "two.sided"),
  sort_by = c("p", "none"),
  mc.cores = 1,
  ...
)
```

## Arguments

- x:

  A Position Weight Matrix, universalmotif object or list thereof.
  Alternatively can be a single DataFrame or list of DataFrames as
  returned by
  [getClusterMatches](https://smped.github.io/motifTestR/reference/getClusterMatches.md)
  with `best_only = TRUE`

- stringset:

  An XStringSet. Not required if matches are supplied as x

- binwidth:

  Width of bins across the range to group data into

- abs:

  Use absolute positions around zero to find symmetrical enrichment

- rc:

  logical(1) Also find matches using the reverse complement of each PWM

- min_score:

  The minimum score to return a match

- break_ties:

  Choose how to resolve matches with tied scores

- alt:

  Alternative hypothesis for the binomial test

- sort_by:

  Column to sort results by

- mc.cores:

  Passed to [mclapply](https://rdrr.io/r/parallel/mclapply.html)

- ...:

  Passed to [matchPWM](https://rdrr.io/pkg/Biostrings/man/matchPWM.html)

## Value

A data.frame with columns `start`, `end`, `centre`, `width`,
`total_matches`, `matches_in_region`, `expected`, `enrichment`,
`prop_total`, `p` and `consensus_motif` The total matches represent the
total number of matches within the set of sequences, whilst the number
observed in the final region are also given, along with the proportion
of the total this represents. Enrichment is simply the ratio of observed
to expected based on the expectation of the null hypothesis

The consensus motif across all matches is returned as a Position
Frequency Matrix (PFM) using
[consensusMatrix](https://rdrr.io/pkg/Biostrings/man/letterFrequency.html).

## Details

This is a reimplementation of
[testMotifPos](https://smped.github.io/motifTestR/reference/testMotifPos.md)
for sets of motifs which have been clustered for similarity. The
positions test the bias of any motifs within the cluster given that
overlapping matches are only counted once, and with the match retained
being the one with the highest relative score.

It should also be noted that some motif clusters will contain PWMs of
varying length. When finding positional bias, the widest motif is taken
as the width for all, and any matches from narrower motifs outside of
the range allowed by wider motifs are discarded. This reduction in
signal will make a small difference in the outer bins, but is not
considered to be problematic for the larger analysis.

## Examples

``` r
## Load the example PWM
data("ex_pfm")
## Load the example sequences
data("ar_er_seq")

## Cluster the motifs
cl <- list(A = ex_pfm[1], B = ex_pfm[2:3])

## Get the best match and use this data
matches <- getClusterMatches(cl, ar_er_seq, best_only = TRUE)
## Test for enrichment in any position
testClusterPos(matches)
#>   start end centre width total_matches matches_in_region  expected enrichment
#> A    -5   5      0    10           175                23  4.521964   5.086286
#> B   -35  25     -5    60           595                92 46.484375   1.979160
#>   prop_total            p          fdr consensus_motif
#> A  0.1314286 9.269662e-09 1.853932e-08    102, 5, ....
#> B  0.1546218 1.968539e-03 1.968539e-03    15, 5, 3....

## Or just pass the clustered matrices
## Here we've set abs = TRUE to test absolute distance from the centre
testClusterPos(cl, ar_er_seq, abs = TRUE, binwidth = 10)
#>   start end centre width total_matches matches_in_region  expected enrichment
#> A     0  10      5    10           175                34  9.067358   3.749714
#> B     0  20     10    20           595               116 61.979167   1.871597
#>   prop_total            p          fdr consensus_motif
#> A  0.1942857 5.538352e-10 1.107670e-09    102, 5, ....
#> B  0.1949580 3.866257e-05 3.866257e-05    15, 5, 3....
```

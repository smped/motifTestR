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
  min_score = "80%",
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
#>   start end centre width total_matches matches_in_region expected enrichment
#> B   -35  75     20   110           114                45 17.81250   2.526316
#> A  -195 135    -30   330            22                22  9.03876   2.433962
#>   prop_total          p        fdr consensus_motif
#> B  0.3947368 0.03553632 0.07107263    0, 0, 0,....
#> A  1.0000000 0.83704449 0.83704449    15, 0, 7....

## Or just pass the clustered matrices
## Here we've set abs = TRUE to test absolute distance from the centre
testClusterPos(cl, ar_er_seq, abs = TRUE, binwidth = 10)
#>   start end centre width total_matches matches_in_region expected enrichment
#> B    10  20     15    10           114                17 5.937500   2.863158
#> A     0 140     70   140            22                16 7.979275   2.005195
#>   prop_total           p         fdr consensus_motif
#> B  0.1491228 0.001779427 0.003558853    0, 0, 0,....
#> A  0.7272727 0.657021567 0.657021567    15, 0, 7....
```

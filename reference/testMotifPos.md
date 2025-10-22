# Test for a Uniform Distribution across a set of best matches

Test for a Uniform Distribution across a set of best matches

## Usage

``` r
testMotifPos(
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
  [getPwmMatches](https://smped.github.io/motifTestR/reference/getPwmMatches.md)
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

This function tests for an even positional spread of motif matches
across a set of sequences, using the assumption (i.e. H~0~) that if
there is no positional bias, matches will be evenly distributed across
all positions within a set of sequences. Conversely, if there is
positional bias, typically but not necessarily near the centre of a
range, this function intends to detect this signal, as a rejection of
the null hypothesis.

Input can be provided as the output from
[getPwmMatches](https://smped.github.io/motifTestR/reference/getPwmMatches.md)
setting `best_only = TRUE` if these matches have already been
identified. If choosing to provide this object to the argument
`matches`, nothing is required for the arguments `pwm`, `stringset`,
`rc`, `min_score` or `break_ties` Otherwise, a Position Weight Matrix
(PWM) and an `XStringSet` are required, along with the relevant
arguments, with best matches identified within the function.

The set of best matches are then grouped into bins along the range, with
the central bin containing zero, and tallied. Setting `abs` to `TRUE`
will set all positions from the centre as *absolute values*, returning
counts purely as bins with distances from zero, marking this as an
inclusive lower bound. Motif alignments are assigned into bins based on
the central position of the match, as provided in the column
`from_centre` when calling
[getPwmMatches](https://smped.github.io/motifTestR/reference/getPwmMatches.md).

The [binom.test](https://rdrr.io/r/stats/binom.test.html) is performed
on each bin using the alternative hypothesis, with the returned p-values
across all bins combined using the Harmonic Mean p-value (HMP) (See
[p.hmp](https://rdrr.io/pkg/harmonicmeanp/man/p.hmp.html)). All bins
with raw p-values below the HMP are identified and the returned values
for start, end, centre, width, matches in region, expected and
enrichment are across this set of bins. The expectation is that where a
positional bias is evident, this will be a narrow range containing a
non-trivial proportion of the total matches.

It should also be noted that
[`binom.test()`](https://rdrr.io/r/stats/binom.test.html) can return
p-values of zero, as beyond machine precision. In these instances, zero
p-values are excluded from calculation of the HMP. This will give a very
slight conservative bias, and assumes that for these extreme cases,
neighbouring bins are highly likely to also return extremely low
p-values and no significance will be lost.

## Examples

``` r
## Load the example PWM
data("ex_pfm")
esr1 <- ex_pfm$ESR1

## Load the example sequences
data("ar_er_seq")

## Get the best match and use this data
matches <- getPwmMatches(esr1, ar_er_seq, best_only = TRUE)
## Test for enrichment in any position
testMotifPos(matches)
#>   start end centre width total_matches matches_in_region expected enrichment
#> 1  -195 135    -30   330            22                22  9.03876   2.433962
#>   prop_total         p       fdr consensus_motif
#> 1          1 0.8370445 0.8370445    15, 0, 7....

## Provide a list of PWMs, testing for distance from zero
testMotifPos(ex_pfm, ar_er_seq, abs = TRUE, binwidth = 10)
#>       start end centre width total_matches matches_in_region expected
#> FOXA1    10  20     15    10           107                16 5.487179
#> ZN143    20  40     30    20            15                 8 1.578947
#> ESR1      0 140     70   140            22                16 7.979275
#> ANDR     10 140     75   130             8                 8 2.916667
#> ZN281    10 180     95   170            13                11 6.062176
#>       enrichment prop_total           p        fdr consensus_motif
#> FOXA1   2.915888  0.1495327 0.002314982 0.01157491    0, 0, 0,....
#> ZN143   5.066667  0.5333333 0.086140631 0.21535158    3, 1, 6,....
#> ESR1    2.005195  0.7272727 0.657021567 0.94718429    15, 0, 7....
#> ANDR    2.742857  1.0000000 0.900392478 0.94718429    0, 0, 0,....
#> ZN281   1.814530  0.8461538 0.947184286 0.94718429    2, 1, 9,....

```

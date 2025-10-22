# motifTestR: Perform Key Analyses on Transcription Factor Binding Motifs

The package `motifTestR` has been designed for two primary analyses of
TFBMs, testing for positional bias and overall enrichment.

## Details

The package `motifTestR` provides two primary functions for testing
TFBMs within a set of sequences

- [`testMotifPos()`](https://smped.github.io/motifTestR/reference/testMotifPos.md)
  for detecting positional bias within a set of test sequences

- [`testMotifEnrich()`](https://smped.github.io/motifTestR/reference/testMotifEnrich.md)
  for testing overall enrichment of a TFBM within a set of test
  sequences

Motifs are also able to be clustered for analysis as a cluster, or for
grouping results. Clusters from external approaches can also be
incorporated.

- [`testClusterPos()`](https://smped.github.io/motifTestR/reference/testClusterPos.md)
  for detecting positional bias for matches to any motif annotated to a
  cluster, within a set of test sequences

- [`testClusterEnrich()`](https://smped.github.io/motifTestR/reference/testClusterEnrich.md)
  for testing overall enrichment of any TFBM annotated to a cluster,
  within a set of test sequences

The main functions rely on lower-level functions such as:

- [`countPwmMatches()`](https://smped.github.io/motifTestR/reference/countPwmMatches.md)
  simply counts the number of matches within an `XStringSet`

- [`getPwmMatches()`](https://smped.github.io/motifTestR/reference/getPwmMatches.md)
  returns the position of matches within an `XStringSet`

- [`countClusterMatches()`](https://smped.github.io/motifTestR/reference/getClusterMatches.md)
  simply counts the number of matches to motifs annotated to a cluster
  within an `XStringSet`

- [`getClusterMatches()`](https://smped.github.io/motifTestR/reference/getClusterMatches.md)
  returns the position of matches to motifs annotated to a cluster
  within an `XStringSet`

- [`makeRMRanges()`](https://smped.github.io/motifTestR/reference/makeRMRanges-methods.md)
  which produces a set of random, matching ranges based on key
  characteristics of the set of test sequences/ranges

A simple utility function is provided to enable visualisation of results

- [`plotMatchPos()`](https://smped.github.io/motifTestR/reference/plotMatchPos.md)
  enables visualisation of the matches within a set of sequences using
  multiple strategies

## See also

Useful links:

- <https://github.com/smped/motifTestR>

- Report bugs at <https://github.com/smped/motifTestR/issues>

## Author

Stevie Pederson

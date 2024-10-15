# motifTestR <img id="motiftestr_logo" alt="MotifTestR Logo" src="man/figures/favicon.png" align="right" width = "125" />

<!-- badges: start -->
[![Build Status](https://github.com/smped/motifTestR/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/smped/motifTestR/actions)
[![Repo Status](https://img.shields.io/badge/repo%20status-Active-green.svg)](https://shields.io/)
[![Codecov test coverage](https://codecov.io/gh/smped/motifTestR/branch/gh-actions/graph/badge.svg)](https://codecov.io/gh/smped/motifTestR?branch=gh-actions)
<!-- badges: end -->


The package `motifTestR` provides a small set of functions for analysing transcription factor binding motifs (TFBMs).
Testing for positional bias is enabled using a novel approach, and testing for enrichment relative to a set of background sequences is enabled using multiple statistical models.

Testing for positional bias is intended to be an R native alternative to [CentriMo](https://meme-suite.org/meme/doc/centrimo.html) from the MEME-suite, and will detect any deviation of note from an even distribution across the width of sequences being tested, not just the centrality of motifs.
Given the conventional statistical approach taken, results are easily interpretable directly through adjusted p-values.
Enrichment testing follows well-worn modelling and iterative strategies, and instead offers a novel approach to selection of control, or background, sequences.

To install the stable version of transmogR from Bioconductor please try the following.

``` r
if (!require("BiocManager")) {
  install.packages("BiocManager")
}
BiocManager::install("motifTestR")
```

Alternatively, the latest build can be installed using

``` r
if (!require("BiocManager")) {
  install.packages("BiocManager")
}
BiocManager::install("smped/motifTestR")
```
